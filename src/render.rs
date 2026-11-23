//! Writing a formula out: as ROOT text (what a file stores and ROOT
//! compiles), as readable text (`Display`), and as LaTeX.

use std::fmt::Write;

use crate::expr::{BinaryOp, Constant, Expr};
use crate::func::Function;

/// How to write a tree as text.
#[derive(Clone, Copy, PartialEq)]
pub(crate) enum Style {
    /// What ROOT's `TFormula` compiles: parameters by index (`[2]`),
    /// `pow(a, b)`, and ROOT's function names.
    Root,
    /// As `Root`, with parameters as ROOT itself writes them (`[p2]`): for a
    /// formula that uses every parameter. (ROOT's `[2]` also defines `p0` and
    /// `p1`; `[p2]` does not.)
    RootP,
    /// As `Root`, with every parameter by its name (`[Mean]`, `[p2]`): what
    /// a ROOT file stores as `fFormula`, beside the map from each name to
    /// its index.
    RootNamed,
    /// For people: parameters by name (`[Mean]`, or `[2]` for `p2`) and `^`.
    Display,
}

/// Binding strength of a node when printed, higher binds tighter.
fn prec(e: &Expr) -> u8 {
    match e {
        Expr::Cond(..) => 1,
        Expr::Binary(op, ..) => op_prec(*op),
        Expr::Neg(_) | Expr::Not(_) => 8,
        Expr::Num(v) if v.is_sign_negative() => 8,
        _ => 10,
    }
}

fn op_prec(op: BinaryOp) -> u8 {
    match op {
        BinaryOp::Or => 2,
        BinaryOp::And => 3,
        BinaryOp::Eq | BinaryOp::Ne => 4,
        BinaryOp::Lt | BinaryOp::Gt | BinaryOp::Le | BinaryOp::Ge => 5,
        BinaryOp::Add | BinaryOp::Sub => 6,
        BinaryOp::Mul | BinaryOp::Div => 7,
        BinaryOp::Pow => 9,
        // `#[non_exhaustive]` is for users: every operator is above.
    }
}

fn op_text(op: BinaryOp) -> &'static str {
    match op {
        BinaryOp::Add => "+",
        BinaryOp::Sub => "-",
        BinaryOp::Mul => "*",
        BinaryOp::Div => "/",
        BinaryOp::Pow => "^",
        BinaryOp::Lt => "<",
        BinaryOp::Gt => ">",
        BinaryOp::Le => "<=",
        BinaryOp::Ge => ">=",
        BinaryOp::Eq => "==",
        BinaryOp::Ne => "!=",
        BinaryOp::And => "&&",
        BinaryOp::Or => "||",
    }
}

/// Whether `e` is written starting with a minus sign.
fn starts_negative(e: &Expr) -> bool {
    match e {
        Expr::Neg(_) => true,
        Expr::Num(v) => v.is_sign_negative(),
        Expr::Binary(_, a, _) => starts_negative(a),
        _ => false,
    }
}

/// A number as C++ reads it back exactly: the shortest digits that
/// round-trip, a whole number as an integer, or (`in_call`) as a double
/// (`1.`), so an overloaded or templated ROOT function sees a `double`.
fn number(out: &mut String, v: f64, in_call: bool) {
    if v.is_nan() {
        out.push_str("TMath::QuietNaN()");
    } else if v.is_infinite() {
        out.push_str(if v > 0.0 {
            "TMath::Infinity()"
        } else {
            "-TMath::Infinity()"
        });
    } else if v.fract() == 0.0 && v.abs() < 1e15 {
        let _ = write!(out, "{v:.0}");
        if in_call {
            out.push('.');
        }
    } else {
        let _ = write!(out, "{v:?}");
    }
}

pub(crate) fn text(e: &Expr, params: &[String], style: Style) -> String {
    let mut out = String::new();
    write_text(&mut out, e, params, style, false);
    out
}

fn write_text(out: &mut String, e: &Expr, params: &[String], style: Style, in_call: bool) {
    let child = |out: &mut String, c: &Expr, parens: bool| {
        if parens {
            out.push('(');
        }
        write_text(out, c, params, style, false);
        if parens {
            out.push(')');
        }
    };
    match e {
        Expr::Num(v) => number(out, *v, in_call),
        // ROOT's `infinity` becomes an undeclared `inf` in the C++ it compiles.
        Expr::Const(Constant::Infinity) if style != Style::Display => {
            out.push_str("TMath::Infinity()")
        }
        Expr::Const(c) => out.push_str(c.root_name()),
        Expr::Var(i) => out.push_str(["x", "y", "z", "t"].get(*i).copied().unwrap_or("x")),
        Expr::Param(i) => {
            let name = params.get(*i).map_or("", String::as_str);
            let positional = name.strip_prefix('p') == Some(&i.to_string());
            if style == Style::RootP {
                let _ = write!(out, "[p{i}]");
            } else if style == Style::RootNamed && !name.is_empty() {
                let _ = write!(out, "[{name}]");
            } else if style == Style::Root || positional || name.is_empty() {
                let _ = write!(out, "[{i}]");
            } else {
                let _ = write!(out, "[{name}]");
            }
        }
        Expr::Neg(a) | Expr::Not(a) => {
            let minus = matches!(e, Expr::Neg(_));
            out.push(if minus { '-' } else { '!' });
            // `- -x`, not `--x` (C++'s decrement) or `-(-x)` (one more level
            // of nesting for every minus).
            if minus && starts_negative(a) && prec(a) >= 8 {
                out.push(' ');
            }
            child(out, a, prec(a) < 8);
        }
        Expr::Binary(BinaryOp::Pow, a, b) if style != Style::Display => {
            out.push_str("pow(");
            write_text(out, a, params, style, true);
            out.push(',');
            write_text(out, b, params, style, true);
            out.push(')');
        }
        Expr::Binary(op, a, b) => {
            let p = op_prec(*op);
            let right_assoc = *op == BinaryOp::Pow;
            let left_parens = if right_assoc {
                prec(a) <= p
            } else {
                prec(a) < p
            };
            let right_parens = if right_assoc {
                prec(b) < p
            } else {
                prec(b) <= p
            };
            child(out, a, left_parens);
            out.push_str(op_text(*op));
            if !right_parens && starts_negative(b) {
                out.push(' ');
            }
            child(out, b, right_parens);
        }
        Expr::Call(Function::Signum, args) if style != Style::Display => {
            out.push_str("TMath::Sign(1.,");
            write_args(out, args, params, style);
            out.push(')');
        }
        Expr::Call(f, args) => {
            let name = if style == Style::Display && *f == Function::Signum {
                "sign".into()
            } else {
                f.root_name()
            };
            out.push_str(&name);
            out.push('(');
            write_args(out, args, params, style);
            out.push(')');
        }
        Expr::Cond(c, t, f) => {
            child(out, c, prec(c) <= 1);
            out.push('?');
            child(out, t, prec(t) <= 1);
            out.push(':');
            child(out, f, false);
        }
    }
}

fn write_args(out: &mut String, args: &[Expr], params: &[String], style: Style) {
    for (i, a) in args.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        write_text(out, a, params, style, style != Style::Display);
    }
}

// --- LaTeX ------------------------------------------------------------------

/// `prec`, with a fraction as tight as a name: `\frac` needs no parentheses.
fn lprec(e: &Expr) -> u8 {
    if matches!(e, Expr::Binary(BinaryOp::Div, ..)) {
        10
    } else {
        prec(e)
    }
}

pub(crate) fn latex(e: &Expr, params: &[String]) -> String {
    let mut out = String::new();
    write_latex(&mut out, e, params);
    out
}

fn latex_param(out: &mut String, name: &str, i: usize) {
    if let Some(n) = name
        .strip_prefix('p')
        .filter(|d| d.parse::<usize>().is_ok())
    {
        let _ = write!(out, "p_{{{n}}}");
    } else if name.is_empty() {
        let _ = write!(out, "p_{{{i}}}");
    } else if name.chars().count() == 1 {
        out.push_str(name);
    } else {
        let escaped = name.replace('_', r"\_");
        let _ = write!(out, r"\mathrm{{{escaped}}}");
    }
}

fn latex_const(c: Constant) -> &'static str {
    match c {
        Constant::Pi => r"\pi",
        Constant::E => "e",
        Constant::Sqrt2 => r"\sqrt{2}",
        Constant::Infinity => r"\infty",
        Constant::Ln10 => r"\ln 10",
        Constant::LogE => r"\log_{10} e",
        Constant::C => "c",
        Constant::G => "G",
        Constant::H => "h",
        Constant::K => "k",
        Constant::Sigma => r"\sigma",
        Constant::R => "R",
        Constant::EulerGamma => r"\gamma",
        Constant::True => "1",
        Constant::False => "0",
    }
}

fn write_latex(out: &mut String, e: &Expr, params: &[String]) {
    let paren = |out: &mut String, c: &Expr| {
        out.push_str(r"\left(");
        write_latex(out, c, params);
        out.push_str(r"\right)");
    };
    let child = |out: &mut String, c: &Expr, parens: bool| {
        if parens {
            paren(out, c);
        } else {
            write_latex(out, c, params);
        }
    };
    match e {
        Expr::Num(v) => {
            if v.fract() == 0.0 && v.abs() < 1e15 {
                let _ = write!(out, "{v:.0}");
            } else {
                let _ = write!(out, "{v}");
            }
        }
        Expr::Const(c) => out.push_str(latex_const(*c)),
        Expr::Var(i) => out.push_str(["x", "y", "z", "t"].get(*i).copied().unwrap_or("x")),
        Expr::Param(i) => latex_param(out, params.get(*i).map_or("", String::as_str), *i),
        Expr::Neg(a) => {
            out.push('-');
            child(out, a, lprec(a) < 8 || starts_negative(a));
        }
        Expr::Not(a) => {
            out.push_str(r"\lnot ");
            child(out, a, lprec(a) < 8);
        }
        Expr::Binary(BinaryOp::Div, a, b) => {
            out.push_str(r"\frac{");
            write_latex(out, a, params);
            out.push_str("}{");
            write_latex(out, b, params);
            out.push('}');
        }
        Expr::Binary(BinaryOp::Pow, a, b) => {
            out.push('{');
            child(out, a, lprec(a) <= 9 || matches!(**a, Expr::Call(..)));
            out.push_str("}^{");
            write_latex(out, b, params);
            out.push('}');
        }
        Expr::Binary(op, a, b) => {
            let p = op_prec(*op);
            let symbol = match op {
                BinaryOp::Add => " + ",
                BinaryOp::Sub => " - ",
                BinaryOp::Mul => r" \cdot ",
                BinaryOp::Lt => " < ",
                BinaryOp::Gt => " > ",
                BinaryOp::Le => r" \leq ",
                BinaryOp::Ge => r" \geq ",
                BinaryOp::Eq => " = ",
                BinaryOp::Ne => r" \neq ",
                BinaryOp::And => r" \land ",
                BinaryOp::Or => r" \lor ",
                _ => " ? ",
            };
            child(out, a, lprec(a) < p);
            out.push_str(symbol);
            child(out, b, lprec(b) <= p || starts_negative(b));
        }
        Expr::Call(f, args) => latex_call(out, *f, args, params),
        Expr::Cond(c, t, f) => {
            out.push_str(r"\left(");
            write_latex(out, c, params);
            out.push_str(r" \;?\; ");
            write_latex(out, t, params);
            out.push_str(r" : ");
            write_latex(out, f, params);
            out.push_str(r"\right)");
        }
    }
}

fn latex_call(out: &mut String, f: Function, args: &[Expr], params: &[String]) {
    let arg = |out: &mut String, i: usize| {
        if let Some(a) = args.get(i) {
            write_latex(out, a, params);
        }
    };
    let wrapped = |out: &mut String, open: &str, close: &str| {
        out.push_str(open);
        arg(out, 0);
        out.push_str(close);
    };
    let named = |out: &mut String, name: &str| {
        out.push_str(name);
        out.push_str(r"\left(");
        for i in 0..args.len() {
            if i > 0 {
                out.push_str(", ");
            }
            arg(out, i);
        }
        out.push_str(r"\right)");
    };
    match f {
        Function::Sqrt => wrapped(out, r"\sqrt{", "}"),
        Function::Cbrt => wrapped(out, r"\sqrt[3]{", "}"),
        Function::Abs => wrapped(out, r"\left|", r"\right|"),
        Function::Floor => wrapped(out, r"\left\lfloor ", r"\right\rfloor"),
        Function::Ceil => wrapped(out, r"\left\lceil ", r"\right\rceil"),
        Function::Exp => wrapped(out, "e^{", "}"),
        Function::Sq => {
            out.push('{');
            if let Some(a) = args.first() {
                if lprec(a) <= 9 {
                    out.push_str(r"\left(");
                    write_latex(out, a, params);
                    out.push_str(r"\right)");
                } else {
                    write_latex(out, a, params);
                }
            }
            out.push_str("}^{2}");
        }
        Function::Pow => {
            out.push_str(r"{\left(");
            arg(out, 0);
            out.push_str(r"\right)}^{");
            arg(out, 1);
            out.push('}');
        }
        Function::Sin => named(out, r"\sin"),
        Function::Cos => named(out, r"\cos"),
        Function::Tan => named(out, r"\tan"),
        Function::ASin => named(out, r"\arcsin"),
        Function::ACos => named(out, r"\arccos"),
        Function::ATan => named(out, r"\arctan"),
        Function::SinH => named(out, r"\sinh"),
        Function::CosH => named(out, r"\cosh"),
        Function::TanH => named(out, r"\tanh"),
        Function::Log => named(out, r"\ln"),
        Function::Log10 => named(out, r"\log_{10}"),
        Function::Log2 => named(out, r"\log_{2}"),
        Function::LGamma => named(out, r"\ln\Gamma"),
        Function::Erf => named(out, r"\operatorname{erf}"),
        Function::Erfc => named(out, r"\operatorname{erfc}"),
        Function::Min => named(out, r"\min"),
        Function::Max => named(out, r"\max"),
        Function::Signum => named(out, r"\operatorname{sgn}"),
        Function::Gaus | Function::NormalPdf => named(out, r"\operatorname{Gauss}"),
        Function::Landau => named(out, r"\operatorname{Landau}"),
        Function::BreitWigner | Function::BreitWignerPdf => named(out, r"\operatorname{BW}"),
        Function::CrystalBall => named(out, r"\operatorname{CB}"),
        Function::Chebyshev(n) => named(out, &format!(r"\operatorname{{T}}_{{{n}}}")),
        other => {
            let name = other.root_name();
            let bare = name.rsplit("::").next().unwrap_or(&name).to_owned();
            named(
                out,
                &format!(r"\operatorname{{{}}}", bare.replace('_', r"\_")),
            );
        }
    }
}
