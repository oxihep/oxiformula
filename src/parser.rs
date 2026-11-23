//! A precedence-climbing (Pratt) parser from tokens to an [`Expr`], with
//! ROOT's shortcuts (`gaus`, `expo`, `polN`, …) and ROOT's numbering of
//! parameters.

use std::cmp::Ordering;
use std::ops::Range;

use crate::ParseError;
use crate::expr::{BinaryOp, CONSTANTS, Constant, Expr};
use crate::func::{self, Function};
use crate::lexer::{Tok, Token, lex};

/// The deepest nesting of parentheses, operators and calls a formula may
/// have. Formulas come from files, so parsing (and evaluating and dropping
/// the tree, which recurse over it) must not be able to exhaust the stack.
pub(crate) const MAX_DEPTH: usize = 200;

/// The largest parameter index a formula may use. It bounds the parameter
/// vector a formula from a file can make the reader allocate.
pub(crate) const MAX_PARAMS: usize = 10_000;

/// Parse result: the tree, the parameter names by index, and the number of
/// variables.
pub(crate) struct Parsed {
    pub(crate) expr: Expr,
    pub(crate) params: Vec<String>,
    pub(crate) ndim: usize,
}

pub(crate) fn parse(src: &str, strict: bool) -> Result<Parsed, ParseError> {
    let toks = lex(src)?;
    let mut p = Parser {
        src,
        toks,
        pos: 0,
        depth: 0,
        strict,
        names: Vec::new(),
        fill_to: None,
        max_var: None,
        whole: None,
    };
    if p.toks.is_empty() {
        return Err(ParseError::new("the formula is empty", 0..0));
    }
    let mut expr = p.expr_ternary()?;
    if let Some(t) = p.toks.get(p.pos) {
        return Err(p.trailing(t));
    }
    let (params, map) = root_order(&p.names, p.fill_to);
    remap(&mut expr, &map);
    // As ROOT's `SetParName(0, "Constant")`, …: the first parameters by
    // index, even when an offset (`gaus(3)`) leaves them unused.
    let mut params = params;
    if let Some(names) = p.whole {
        for (slot, name) in params.iter_mut().zip(names) {
            (*name).clone_into(slot);
        }
    }
    let ndim = p.max_var.map_or(1, |v| v + 1);
    Ok(Parsed { expr, params, ndim })
}

struct Parser<'s> {
    src: &'s str,
    toks: Vec<Token>,
    pos: usize,
    /// Current nesting, checked against [`MAX_DEPTH`].
    depth: usize,
    /// Refuse what ROOT does not accept (`ln`, a one-argument `sign`, …).
    strict: bool,
    /// Parameter names by first appearance; numbered ones as `p<N>`.
    names: Vec<String>,
    /// The highest `[N]` written with digits only: ROOT then defines every
    /// `p0` … `pN`.
    fill_to: Option<usize>,
    max_var: Option<usize>,
    /// The parameter names ROOT gives a formula that is a single shortcut
    /// (`gaus` → `Constant`, `Mean`, `Sigma`).
    whole: Option<&'static [&'static str]>,
}

impl Parser<'_> {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos).map(|t| &t.tok)
    }

    fn span(&self) -> Range<usize> {
        self.toks
            .get(self.pos)
            .map_or(self.src.len()..self.src.len(), |t| t.span.clone())
    }

    fn text(&self, span: &Range<usize>) -> &str {
        self.src.get(span.clone()).unwrap_or("")
    }

    fn eat(&mut self, want: &Tok, what: &str) -> Result<Range<usize>, ParseError> {
        let span = self.span();
        if self.peek() == Some(want) {
            self.pos += 1;
            Ok(span)
        } else if self.pos >= self.toks.len() {
            Err(ParseError::new(
                format!("the formula ends where {what} was expected"),
                span,
            ))
        } else {
            Err(ParseError::new(
                format!("expected {what}, found `{}`", self.text(&span)),
                span,
            ))
        }
    }

    /// The error for a token after a complete expression.
    fn trailing(&self, t: &Token) -> ParseError {
        let e = ParseError::new(
            format!(
                "unexpected `{}` after a complete expression",
                self.text(&t.span)
            ),
            t.span.clone(),
        );
        let before = self.pos.checked_sub(1).and_then(|i| self.toks.get(i));
        let juxtaposed = matches!(
            before.map(|b| &b.tok),
            Some(Tok::Num(_) | Tok::RParen | Tok::Ident(_) | Tok::Bracket(_))
        ) && matches!(
            t.tok,
            Tok::Num(_) | Tok::Ident(_) | Tok::LParen | Tok::Bracket(_)
        );
        if juxtaposed {
            e.with_hint("ROOT has no implicit multiplication: write `2*x`, not `2x`")
        } else {
            e
        }
    }

    /// The lowest precedence: the right-associative `cond ? a : b`.
    fn expr_ternary(&mut self) -> Result<Expr, ParseError> {
        let cond = self.expr_bp(0)?;
        if self.peek() == Some(&Tok::Question) {
            self.pos += 1;
            let then = self.expr_ternary()?;
            self.eat(&Tok::Colon, "`:`")?;
            let otherwise = self.expr_ternary()?;
            Ok(Expr::Cond(
                Box::new(cond),
                Box::new(then),
                Box::new(otherwise),
            ))
        } else {
            Ok(cond)
        }
    }

    /// Precedence climbing over the binary operators. Every recursive path of
    /// the parser goes through here, so nesting is limited here.
    fn expr_bp(&mut self, min_bp: u8) -> Result<Expr, ParseError> {
        if self.depth >= MAX_DEPTH {
            return Err(ParseError::new(
                "the formula is nested too deeply",
                self.span(),
            ));
        }
        self.depth += 1;
        let result = self.expr_bp_inner(min_bp);
        self.depth -= 1;
        result
    }

    fn expr_bp_inner(&mut self, min_bp: u8) -> Result<Expr, ParseError> {
        let mut lhs = self.prefix()?;
        while let Some((op, l_bp, r_bp)) = self.peek().and_then(binop) {
            if l_bp < min_bp {
                break;
            }
            self.pos += 1;
            let rhs = self.expr_bp(r_bp)?;
            lhs = Expr::Binary(op, Box::new(lhs), Box::new(rhs));
        }
        Ok(lhs)
    }

    /// A unary `-`, `+` or `!`, or a primary.
    fn prefix(&mut self) -> Result<Expr, ParseError> {
        match self.peek() {
            Some(Tok::Minus) => {
                self.pos += 1;
                Ok(match self.expr_bp(PREFIX_BP)? {
                    Expr::Num(v) => Expr::Num(-v),
                    e => Expr::Neg(Box::new(e)),
                })
            }
            Some(Tok::Plus) => {
                self.pos += 1;
                self.expr_bp(PREFIX_BP)
            }
            Some(Tok::Not) => {
                self.pos += 1;
                Ok(Expr::Not(Box::new(self.expr_bp(PREFIX_BP)?)))
            }
            _ => self.primary(),
        }
    }

    fn primary(&mut self) -> Result<Expr, ParseError> {
        let span = self.span();
        let Some(tok) = self.peek().cloned() else {
            return Err(ParseError::new("the formula ends too early", span));
        };
        self.pos += 1;
        match tok {
            Tok::Num(v) => Ok(Expr::Num(v)),
            Tok::Bracket(name) => self.param(&name, span),
            Tok::LParen => {
                let e = self.expr_ternary()?;
                self.eat(&Tok::RParen, "`)`")?;
                Ok(e)
            }
            Tok::Ident(name) => self.ident(&name, span),
            _ => Err(ParseError::new(
                format!("expected a value, found `{}`", self.text(&span)),
                span,
            )),
        }
    }

    /// A parameter `[…]`: `[3]` (which, as in ROOT, defines `p0` … `p3`),
    /// `[p3]`, or a name.
    fn param(&mut self, name: &str, span: Range<usize>) -> Result<Expr, ParseError> {
        let canonical = if name.bytes().all(|b| b.is_ascii_digit()) {
            let n = index(name, &span)?;
            self.fill_to = Some(self.fill_to.map_or(n, |f| f.max(n)));
            format!("p{n}")
        } else if let Some(digits) = name
            .strip_prefix('p')
            .filter(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()))
        {
            format!("p{}", index(digits, &span)?)
        } else if name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_') {
            name.to_owned()
        } else {
            return Err(
                ParseError::new(format!("`[{name}]` is not a parameter"), span).with_hint(
                    "a parameter is `[0]`, `[p1]` or a name of letters, digits and `_`; \
                     variables in brackets go after a shortcut, as in `gaus[x,y]`",
                ),
            );
        };
        let i = self
            .names
            .iter()
            .position(|n| *n == canonical)
            .unwrap_or_else(|| {
                self.names.push(canonical);
                self.names.len() - 1
            });
        Ok(Expr::Param(i))
    }

    fn var(&mut self, v: usize) -> Expr {
        self.max_var = Some(self.max_var.map_or(v, |m| m.max(v)));
        Expr::Var(v)
    }

    /// A name: a variable, a constant, a shortcut or a function call.
    ///
    /// This and [`call`](Self::call) stay on the stack while a nested
    /// formula is parsed, so they keep their frames small: the error
    /// messages are built in functions of their own.
    fn ident(&mut self, name: &str, span: Range<usize>) -> Result<Expr, ParseError> {
        if let Some(v) = var_index(name) {
            // `x[1]` is ROOT's other spelling of `y`.
            if v == 0 && matches!(self.peek(), Some(Tok::Bracket(_))) {
                return self.indexed_x();
            }
            return Ok(self.var(v));
        }
        if let Some(e) = self.constant(name, &span)? {
            return Ok(e);
        }
        if let Some(e) = self.shortcut(name, &span)? {
            return Ok(e);
        }
        self.call(name, span)
    }

    /// `x[i]`, the bracket next.
    #[inline(never)]
    fn indexed_x(&mut self) -> Result<Expr, ParseError> {
        let bracket = self.span();
        let Some(Tok::Bracket(i)) = self.peek().cloned() else {
            return Err(ParseError::new("expected `[`", bracket));
        };
        self.pos += 1;
        let v = i.parse::<usize>().ok().filter(|v| *v < 4).ok_or_else(|| {
            ParseError::new(format!("`x[{i}]` is not a variable"), bracket)
                .with_hint("the variables are x[0] to x[3], or x, y, z and t")
        })?;
        Ok(self.var(v))
    }

    /// A constant: `pi`, or `TMath::Pi()`.
    fn constant(&mut self, name: &str, span: &Range<usize>) -> Result<Option<Expr>, ParseError> {
        let bare = CONSTANTS.iter().find(|(n, _, _)| *n == name);
        let call = CONSTANTS
            .iter()
            .find(|(_, f, _)| !f.is_empty() && *f == name);
        let Some((_, _, c)) = bare.or(call) else {
            return Ok(None);
        };
        let called = self.peek() == Some(&Tok::LParen);
        if called {
            self.pos += 1;
            self.eat(&Tok::RParen, "`)`: a constant takes no arguments")?;
        }
        if call.is_some() && !called {
            return Err(
                ParseError::new(format!("`{name}` needs `()`"), span.clone())
                    .with_hint(format!("write `{name}()`, or `{}`", c.root_name())),
            );
        }
        if *c == Constant::Infinity && bare.is_some() && self.strict {
            return Err(
                ParseError::new("ROOT cannot compile `infinity`", span.clone())
                    .with_hint("ROOT turns it into an undeclared `inf`: use `TMath::Infinity()`"),
            );
        }
        if bare.is_some() && called && self.strict {
            return Err(
                ParseError::new(format!("ROOT's `{name}` is a constant"), span.clone())
                    .with_hint(format!("write `{name}` without `()`")),
            );
        }
        Ok(Some(Expr::Const(*c)))
    }

    /// A function call `name(args…)`.
    fn call(&mut self, name: &str, span: Range<usize>) -> Result<Expr, ParseError> {
        if self.peek() != Some(&Tok::LParen) {
            return Err(unknown_name(name, span));
        }
        self.pos += 1;
        let mut args = Vec::new();
        if self.peek() != Some(&Tok::RParen) {
            args.push(self.expr_ternary()?);
            while self.peek() == Some(&Tok::Comma) {
                self.pos += 1;
                args.push(self.expr_ternary()?);
            }
        }
        let close = self.eat(&Tok::RParen, "`,` or `)`")?;
        self.resolve(name, span, close.end, args)
    }

    /// The function `name(args…)` calls.
    #[inline(never)]
    fn resolve(
        &self,
        name: &str,
        span: Range<usize>,
        end: usize,
        args: Vec<Expr>,
    ) -> Result<Expr, ParseError> {
        let Some((f, root_ok)) = func::lookup(name, args.len()) else {
            return Err(match func::known_name(name) {
                Some(f) => {
                    let (lo, hi) = f.arity();
                    let count = if lo == hi {
                        format!("{lo}")
                    } else {
                        format!("{lo} to {hi}")
                    };
                    ParseError::new(
                        format!(
                            "`{name}` takes {count} argument{}, not {}",
                            if hi == 1 { "" } else { "s" },
                            args.len()
                        ),
                        span.start..end,
                    )
                }
                None => unknown_name(name, span),
            });
        };
        if self.strict && !root_ok {
            let hint = match f {
                Function::Signum => {
                    "ROOT's `sign(a, b)` takes two arguments: `TMath::Sign(1., x)` is the sign of x"
                        .to_owned()
                }
                _ => format!("use `{}`", f.root_name()),
            };
            return Err(
                ParseError::new(format!("ROOT has no function `{name}`"), span).with_hint(hint),
            );
        }
        Ok(match (f, args.as_slice()) {
            // One tree for one formula: `pow(a, b)` is `a^b`, and ROOT's
            // `TMath::Sign(1., x)` is the sign of x.
            (Function::Pow, [_, _]) => {
                let mut args = args.into_iter();
                let (a, b) = (args.next(), args.next());
                match (a, b) {
                    (Some(a), Some(b)) => Expr::Binary(BinaryOp::Pow, Box::new(a), Box::new(b)),
                    _ => return Err(ParseError::new("`pow` takes 2 arguments", span)),
                }
            }
            (Function::Sign, [Expr::Num(one), _]) if *one == 1.0 => {
                Expr::Call(Function::Signum, args.into_iter().skip(1).collect())
            }
            _ => Expr::Call(f, args),
        })
    }

    /// Expand a ROOT shortcut, if `name` is one, as ROOT does: `gaus`,
    /// `gausn`, `expo`, `landau`, `landaun`, `crystalball`, `breitwigner`,
    /// `bigaus`, `chebN` and `polN`, with variables (`xygaus`, `gaus[x,y]`,
    /// `ypol2`, `pol2(y, 0)`) and a first parameter (`gaus(3)`).
    fn shortcut(&mut self, name: &str, span: &Range<usize>) -> Result<Option<Expr>, ParseError> {
        let Some((kind, prefix)) = Shortcut::from_name(name) else {
            return Ok(None);
        };
        let first = self.pos - 1;
        let mut vars: Vec<String> = prefix
            .map(|p| p.chars().map(String::from).collect())
            .unwrap_or_default();
        let mut offset = 0;
        let mut named: Vec<String> = Vec::new();
        if vars.is_empty()
            && !matches!(kind, Shortcut::Pol(_))
            && let Some(Tok::Bracket(list)) = self.peek().cloned()
        {
            self.pos += 1;
            vars = list.split(',').map(|v| v.trim().to_owned()).collect();
        }
        if self.peek() == Some(&Tok::LParen) {
            let open = self.span();
            self.pos += 1;
            match (self.peek().cloned(), kind) {
                (Some(Tok::Num(v)), _) => {
                    self.pos += 1;
                    offset = whole_index(v, &open)?;
                }
                (Some(Tok::Ident(var)), Shortcut::Pol(_)) if var_index(&var).is_some() => {
                    // ROOT's newer `pol2(y, 3)` and `pol2(x, [a], [b], [c])`.
                    self.pos += 1;
                    vars = vec![var];
                    self.eat(&Tok::Comma, "`,` and the first parameter")?;
                    if let Some(Tok::Num(v)) = self.peek().cloned() {
                        self.pos += 1;
                        offset = whole_index(v, &open)?;
                    } else {
                        loop {
                            let at = self.span();
                            match self.peek().cloned() {
                                Some(Tok::Bracket(p)) => {
                                    self.pos += 1;
                                    named.push(p);
                                }
                                _ => {
                                    return Err(ParseError::new(
                                        "expected a parameter such as `[a]`",
                                        at,
                                    ));
                                }
                            }
                            if self.peek() != Some(&Tok::Comma) {
                                break;
                            }
                            self.pos += 1;
                        }
                    }
                }
                _ => {
                    return Err(ParseError::new(
                        format!("`{name}(…)` takes the index of its first parameter"),
                        open,
                    )
                    .with_hint(format!(
                        "`{name}(3)` starts at `[3]`; for another variable write `{}`",
                        kind.variable_example()
                    )));
                }
            }
            self.eat(&Tok::RParen, "`)`")?;
        }
        let dim = if vars.is_empty() {
            kind.default_dim()
        } else {
            vars.len()
        };
        if vars.is_empty() {
            vars = ["x", "y", "z"][..dim]
                .iter()
                .map(|v| (*v).to_owned())
                .collect();
        }
        for v in &vars {
            if var_index(v).is_none() {
                return Err(ParseError::new(
                    format!("`{v}` is not a variable of `{name}`"),
                    span.clone(),
                )
                .with_hint("the variables are x, y, z and t"));
            }
        }
        let body = kind
            .body(dim, offset, &vars, &named)
            .map_err(|message| ParseError::new(format!("`{name}`: {message}"), span.clone()))?;
        let expr = self.sub_parse(&body, span)?;
        if first == 0 && self.pos == self.toks.len() {
            self.whole = kind.names(dim);
        }
        Ok(Some(expr))
    }

    /// Parse a shortcut's expansion with this parser's parameters.
    fn sub_parse(&mut self, body: &str, span: &Range<usize>) -> Result<Expr, ParseError> {
        let internal = |e: ParseError| {
            ParseError::new(
                format!("internal error expanding a shortcut: {}", e.message()),
                span.clone(),
            )
        };
        let toks = lex(body).map_err(internal)?;
        let saved = (std::mem::replace(&mut self.toks, toks), self.pos, self.src);
        self.pos = 0;
        // Spans inside the body point into it, not the formula; errors are
        // reported at the shortcut.
        let result = self.expr_ternary();
        let complete = self.pos == self.toks.len();
        (self.toks, self.pos) = (saved.0, saved.1);
        self.src = saved.2;
        let expr = result.map_err(internal)?;
        if complete {
            Ok(expr)
        } else {
            Err(ParseError::new(
                "internal error expanding a shortcut",
                span.clone(),
            ))
        }
    }
}

/// Binding power of the unary `-`, `+` and `!`: tighter than `*`, looser
/// than `^` (so `-x^2` is `-(x^2)`).
const PREFIX_BP: u8 = 13;

/// `(operator, left_bp, right_bp)` for a binary token. Right-associative `^`
/// has `right_bp < left_bp`. As in C++, `==`/`!=` bind looser than `<`/`>`.
fn binop(t: &Tok) -> Option<(BinaryOp, u8, u8)> {
    Some(match t {
        Tok::Or => (BinaryOp::Or, 1, 2),
        Tok::And => (BinaryOp::And, 3, 4),
        Tok::EqEq => (BinaryOp::Eq, 5, 6),
        Tok::Ne => (BinaryOp::Ne, 5, 6),
        Tok::Lt => (BinaryOp::Lt, 7, 8),
        Tok::Gt => (BinaryOp::Gt, 7, 8),
        Tok::Le => (BinaryOp::Le, 7, 8),
        Tok::Ge => (BinaryOp::Ge, 7, 8),
        Tok::Plus => (BinaryOp::Add, 9, 10),
        Tok::Minus => (BinaryOp::Sub, 9, 10),
        Tok::Star => (BinaryOp::Mul, 11, 12),
        Tok::Slash => (BinaryOp::Div, 11, 12),
        Tok::Caret => (BinaryOp::Pow, 16, 15),
        _ => return None,
    })
}

/// The variable index of `x`, `y`, `z` and `t`.
pub(crate) fn var_index(name: &str) -> Option<usize> {
    match name {
        "x" => Some(0),
        "y" => Some(1),
        "z" => Some(2),
        "t" => Some(3),
        _ => None,
    }
}

fn index(digits: &str, span: &Range<usize>) -> Result<usize, ParseError> {
    digits
        .parse::<usize>()
        .ok()
        .filter(|n| *n < MAX_PARAMS)
        .ok_or_else(|| {
            ParseError::new(
                format!("parameter index {digits} is too large"),
                span.clone(),
            )
            .with_hint(format!("a formula has at most {MAX_PARAMS} parameters"))
        })
}

fn whole_index(v: f64, span: &Range<usize>) -> Result<usize, ParseError> {
    if v >= 0.0 && v.fract() == 0.0 && v < MAX_PARAMS as f64 {
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a whole number in 0..MAX_PARAMS"
        )]
        let n = v as usize;
        Ok(n)
    } else {
        Err(ParseError::new(
            "a shortcut's first parameter is a whole number",
            span.clone(),
        ))
    }
}

fn unknown_name(name: &str, span: Range<usize>) -> ParseError {
    let e = ParseError::new(format!("unknown name `{name}`"), span);
    if let Some(f) = func::known_name(name) {
        return e.with_hint(format!(
            "`{name}` is a function: write `{}(…)`",
            f.root_name()
        ));
    }
    if let Some(close) = func::closest(name) {
        return e.with_hint(format!("did you mean `{close}`?"));
    }
    let constant = CONSTANTS
        .iter()
        .map(|(n, _, _)| (func::edit_distance(name, n), *n))
        .filter(|(d, n)| *d <= 1 && n.len() > 2)
        .min_by_key(|(d, _)| *d);
    if let Some((_, c)) = constant {
        return e.with_hint(format!("did you mean `{c}`?"));
    }
    e.with_hint("the variables are x, y, z and t; parameters are written `[0]` or `[name]`")
}

/// ROOT's order of parameters (`TFormulaParamOrder`): names of the form
/// `p<N>` first, by number, then the others by byte order. Returns the names
/// in index order (with every `p0` … `p<fill_to>` ROOT defines) and, for each
/// name in `names`, its index.
fn root_order(names: &[String], fill_to: Option<usize>) -> (Vec<String>, Vec<usize>) {
    let mut all: Vec<String> = names.to_vec();
    if let Some(n) = fill_to {
        for i in 0..=n {
            let p = format!("p{i}");
            if !all.contains(&p) {
                all.push(p);
            }
        }
    }
    let number = |s: &str| s.strip_prefix('p').and_then(|d| d.parse::<usize>().ok());
    all.sort_by(|a, b| match (number(a), number(b)) {
        (Some(x), Some(y)) => x.cmp(&y),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => a.as_bytes().cmp(b.as_bytes()),
    });
    let map = names
        .iter()
        .map(|n| all.iter().position(|a| a == n).unwrap_or(0))
        .collect();
    (all, map)
}

/// Renumber the parameters of `e` from first-appearance order to `map`.
fn remap(e: &mut Expr, map: &[usize]) {
    match e {
        Expr::Param(i) => *i = map.get(*i).copied().unwrap_or(*i),
        Expr::Num(_) | Expr::Const(_) | Expr::Var(_) => {}
        Expr::Neg(a) | Expr::Not(a) => remap(a, map),
        Expr::Binary(_, a, b) => {
            remap(a, map);
            remap(b, map);
        }
        Expr::Call(_, args) => args.iter_mut().for_each(|a| remap(a, map)),
        Expr::Cond(c, t, f) => {
            remap(c, map);
            remap(t, map);
            remap(f, map);
        }
    }
}

/// A ROOT shortcut.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Shortcut {
    Gaus { norm: bool },
    Expo,
    Landau { norm: bool },
    CrystalBall,
    BreitWigner,
    BiGaus,
    Cheb(u8),
    Pol(usize),
}

impl Shortcut {
    /// The shortcut `name` names, and the variables its prefix gives
    /// (`xy` of `xygaus`).
    fn from_name(name: &str) -> Option<(Shortcut, Option<&'static str>)> {
        let base = |n: &str| -> Option<Shortcut> {
            Some(match n {
                "gaus" => Shortcut::Gaus { norm: false },
                "gausn" => Shortcut::Gaus { norm: true },
                "expo" => Shortcut::Expo,
                "landau" => Shortcut::Landau { norm: false },
                "landaun" => Shortcut::Landau { norm: true },
                "crystalball" => Shortcut::CrystalBall,
                "breitwigner" => Shortcut::BreitWigner,
                "bigaus" => Shortcut::BiGaus,
                _ => {
                    let digits = |d: &str| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit());
                    if let Some(d) = n.strip_prefix("cheb").filter(|d| digits(d)) {
                        Shortcut::Cheb(d.parse().ok().filter(|d| *d <= 10)?)
                    } else {
                        let d = n.strip_prefix("pol").filter(|d| digits(d))?;
                        Shortcut::Pol(d.parse().ok()?)
                    }
                }
            })
        };
        if let Some(s) = base(name) {
            return Some((s, None));
        }
        for prefix in ["xyz", "xy", "x", "y", "z", "t"] {
            let Some(rest) = name.strip_prefix(prefix) else {
                continue;
            };
            let s = base(rest)?;
            let allowed = match s {
                Shortcut::Pol(_) => prefix.len() == 1,
                Shortcut::Gaus { .. } | Shortcut::Expo | Shortcut::Landau { .. } => prefix != "t",
                _ => false,
            };
            return allowed.then_some((s, Some(prefix)));
        }
        None
    }

    fn default_dim(self) -> usize {
        if self == Shortcut::BiGaus { 2 } else { 1 }
    }

    fn variable_example(self) -> &'static str {
        match self {
            Shortcut::Pol(_) => "`pol2(y, 0)` or `ypol2`",
            _ => "`gaus[y]` or `ygaus`",
        }
    }

    /// ROOT's expansion (`TFormula::FillParametrizedFunctions`), with the
    /// variables put in and the parameters moved up by `offset`.
    fn body(
        self,
        dim: usize,
        offset: usize,
        vars: &[String],
        named: &[String],
    ) -> Result<String, String> {
        let template: String = match (self, dim) {
            (Shortcut::Gaus { norm: false }, 1) => {
                "[0]*exp(-0.5*(({V0}-[1])/[2])*(({V0}-[1])/[2]))".into()
            }
            (Shortcut::Gaus { norm: true }, 1) => {
                "[0]*exp(-0.5*(({V0}-[1])/[2])*(({V0}-[1])/[2]))/(sqrt(2*pi)*[2])".into()
            }
            (Shortcut::Gaus { norm: false }, 2) => {
                "[0]*exp(-0.5*(({V0}-[1])/[2])^2 - 0.5*(({V1}-[3])/[4])^2)".into()
            }
            (Shortcut::Gaus { norm: false }, 3) => {
                "[0]*exp(-0.5*(({V0}-[1])/[2])^2 - 0.5*(({V1}-[3])/[4])^2 - 0.5*(({V2}-[5])/[6])^2)"
                    .into()
            }
            (Shortcut::Expo, 1) => "exp([0]+[1]*{V0})".into(),
            (Shortcut::Expo, 2) => "exp([0]+[1]*{V0}+[2]*{V1})".into(),
            (Shortcut::Landau { norm }, 1) => {
                format!("[0]*TMath::Landau({{V0}},[1],[2],{norm})")
            }
            (Shortcut::Landau { norm: false }, 2) => {
                "[0]*TMath::Landau({V0},[1],[2],false)*TMath::Landau({V1},[3],[4],false)".into()
            }
            (Shortcut::Landau { norm: true }, 2) => {
                "TMath::Landau({V0},[0],[1],true)*TMath::Landau({V1},[2],[3],true)".into()
            }
            (Shortcut::CrystalBall, 1) => {
                "[0]*ROOT::Math::crystalball_function({V0},[3],[4],[2],[1])".into()
            }
            (Shortcut::BreitWigner, 1) => "[0]*ROOT::Math::breitwigner_pdf({V0},[2],[1])".into(),
            (Shortcut::BiGaus, 2) => {
                "[0]*ROOT::Math::bigaussian_pdf({V0},{V1},[2],[4],[5],[1],[3])".into()
            }
            (Shortcut::Cheb(n), 1) => {
                let coeffs: Vec<String> = (0..=n).map(|k| format!("[{k}]")).collect();
                format!("ROOT::Math::Chebyshev{n}({{V0}},{})", coeffs.join(","))
            }
            (Shortcut::Pol(n), 1) => {
                if n + offset >= MAX_PARAMS {
                    return Err(format!(
                        "a polynomial of degree {n} has too many parameters"
                    ));
                }
                let v = &vars[0];
                let param = |k: usize| -> Result<String, String> {
                    if named.is_empty() {
                        Ok(format!("[{}]", k + offset))
                    } else {
                        named
                            .get(k)
                            .map(|p| format!("[{p}]"))
                            .ok_or_else(|| format!("needs {} parameter names", n + 1))
                    }
                };
                let mut s = format!("({}", param(0)?);
                for k in 1..=n {
                    s.push('+');
                    s.push_str(&param(k)?);
                    s.push('*');
                    s.push_str(v);
                    if k > 1 {
                        s.push_str(&format!("^{k}"));
                    }
                }
                s.push(')');
                return Ok(s);
            }
            (Shortcut::Gaus { norm: true }, _) => {
                return Err("has no normalized form in more than one dimension".into());
            }
            _ => return Err(format!("has no {dim}-dimensional form")),
        };
        let mut out = String::with_capacity(template.len() + 8);
        let mut rest = template.as_str();
        while let Some(i) = rest.find(['{', '[']) {
            out.push_str(&rest[..i]);
            let close = if rest.as_bytes()[i] == b'{' { '}' } else { ']' };
            let end = i + rest[i..].find(close).unwrap_or(0);
            let inner = &rest[i + 1..end];
            if close == '}' {
                let k: usize = inner[1..].parse().unwrap_or(0);
                out.push_str(&vars[k]);
            } else {
                let k: usize = inner.parse().unwrap_or(0);
                out.push_str(&format!("[{}]", k + offset));
            }
            rest = &rest[end + 1..];
        }
        out.push_str(rest);
        Ok(out)
    }

    /// The parameter names ROOT gives a formula that is this shortcut alone
    /// (`TFormula::SetPredefinedParamNames`).
    fn names(self, dim: usize) -> Option<&'static [&'static str]> {
        Some(match (self, dim) {
            (Shortcut::Gaus { .. }, 1) => &["Constant", "Mean", "Sigma"],
            (Shortcut::Gaus { .. }, 2) => &["Constant", "MeanX", "SigmaX", "MeanY", "SigmaY"],
            (Shortcut::Gaus { .. }, 3) => &[
                "Constant", "MeanX", "SigmaX", "MeanY", "SigmaY", "MeanZ", "SigmaZ",
            ],
            (Shortcut::BiGaus, 2) => &["Constant", "MeanX", "SigmaX", "MeanY", "SigmaY", "Rho"],
            (Shortcut::Expo, 1) => &["Constant", "Slope"],
            (Shortcut::Landau { .. }, 1) => &["Constant", "MPV", "Sigma"],
            (Shortcut::CrystalBall, 1) => &["Constant", "Mean", "Sigma", "Alpha", "N"],
            _ => return None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn root_order_is_roots() {
        let names = |v: &[&str]| v.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>();
        let (all, map) = root_order(&names(&["b", "a"]), None);
        assert_eq!((all, map), (names(&["a", "b"]), vec![1, 0]));
        let (all, _) = root_order(&names(&["mean", "p1"]), None);
        assert_eq!(all, names(&["p1", "mean"]));
        let (all, _) = root_order(&names(&["p10", "p2"]), None);
        assert_eq!(all, names(&["p2", "p10"]));
        let (all, _) = root_order(&names(&["B", "a", "_c"]), None);
        assert_eq!(all, names(&["B", "_c", "a"]));
        let (all, _) = root_order(&names(&["p1", "mean"]), Some(1));
        assert_eq!(all, names(&["p0", "p1", "mean"]));
    }
}
