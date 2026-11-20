//! The expression AST and its evaluator.

/// A binary operator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum BinOp {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
    Lt,
    Gt,
    Le,
    Ge,
    Eq,
    Ne,
    And,
    Or,
}

/// A supported function, resolved from its name at parse time (so evaluation is
/// a cheap `match`, never a string lookup). ROOT's `TMath::` prefix and case are
/// normalized away by [`Func::resolve`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Func {
    Sin,
    Cos,
    Tan,
    ASin,
    ACos,
    ATan,
    SinH,
    CosH,
    TanH,
    ASinH,
    ACosH,
    ATanH,
    Exp,
    Expm1,
    Log,
    Log10,
    Log2,
    Sqrt,
    Cbrt,
    Abs,
    Floor,
    Ceil,
    Sign,
    Erf,
    Erfc,
    LGamma,
    ATan2,
    Pow,
    Hypot,
    Min,
    Max,
    Pi,
}

impl Func {
    /// Resolve a function name (after stripping a `TMath::` prefix and
    /// lowercasing) to a `(function, arity)`, or `None` if unknown.
    pub(crate) fn resolve(name: &str) -> Option<(Func, usize)> {
        let n = name
            .strip_prefix("TMath::")
            .unwrap_or(name)
            .to_ascii_lowercase();
        use Func::*;
        Some(match n.as_str() {
            "sin" => (Sin, 1),
            "cos" => (Cos, 1),
            "tan" => (Tan, 1),
            "asin" => (ASin, 1),
            "acos" => (ACos, 1),
            "atan" => (ATan, 1),
            "sinh" => (SinH, 1),
            "cosh" => (CosH, 1),
            "tanh" => (TanH, 1),
            "asinh" => (ASinH, 1),
            "acosh" => (ACosH, 1),
            "atanh" => (ATanH, 1),
            "exp" => (Exp, 1),
            "expm1" => (Expm1, 1),
            "log" | "ln" => (Log, 1),
            "log10" => (Log10, 1),
            "log2" => (Log2, 1),
            "sqrt" => (Sqrt, 1),
            "cbrt" => (Cbrt, 1),
            "abs" | "fabs" => (Abs, 1),
            "floor" => (Floor, 1),
            "ceil" => (Ceil, 1),
            "sign" => (Sign, 1),
            "erf" => (Erf, 1),
            "erfc" => (Erfc, 1),
            "lgamma" | "lngamma" => (LGamma, 1),
            "atan2" => (ATan2, 2),
            "pow" | "power" => (Pow, 2),
            "hypot" => (Hypot, 2),
            "min" => (Min, 2),
            "max" => (Max, 2),
            "pi" => (Pi, 0),
            _ => return None,
        })
    }
}

/// An expression node. `Var(0/1/2)` is `x`/`y`/`z`; `Par(i)` is parameter `[i]`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Expr {
    Num(f64),
    Var(usize),
    Par(usize),
    Neg(Box<Expr>),
    Bin(BinOp, Box<Expr>, Box<Expr>),
    Call(Func, Vec<Expr>),
    /// `cond ? then : otherwise`.
    Cond(Box<Expr>, Box<Expr>, Box<Expr>),
}

#[inline]
fn truthy(v: f64) -> bool {
    v != 0.0 && !v.is_nan()
}

impl Expr {
    /// Evaluate the node. `vars` is `[x]`/`[x,y]`/`[x,y,z]`; `params` holds the
    /// parameter values. Out-of-range variable/parameter indices read as 0.
    pub(crate) fn eval(&self, vars: &[f64], params: &[f64]) -> f64 {
        match self {
            Expr::Num(v) => *v,
            Expr::Var(i) => vars.get(*i).copied().unwrap_or(0.0),
            Expr::Par(i) => params.get(*i).copied().unwrap_or(0.0),
            Expr::Neg(e) => -e.eval(vars, params),
            Expr::Bin(op, a, b) => {
                let (l, r) = (a.eval(vars, params), b.eval(vars, params));
                match op {
                    BinOp::Add => l + r,
                    BinOp::Sub => l - r,
                    BinOp::Mul => l * r,
                    BinOp::Div => l / r,
                    BinOp::Pow => l.powf(r),
                    BinOp::Lt => f64::from(u8::from(l < r)),
                    BinOp::Gt => f64::from(u8::from(l > r)),
                    BinOp::Le => f64::from(u8::from(l <= r)),
                    BinOp::Ge => f64::from(u8::from(l >= r)),
                    BinOp::Eq => f64::from(u8::from(l == r)),
                    BinOp::Ne => f64::from(u8::from(l != r)),
                    BinOp::And => f64::from(u8::from(truthy(l) && truthy(r))),
                    BinOp::Or => f64::from(u8::from(truthy(l) || truthy(r))),
                }
            }
            Expr::Call(f, args) => {
                let a = |i: usize| args[i].eval(vars, params);
                use Func::*;
                match f {
                    Sin => a(0).sin(),
                    Cos => a(0).cos(),
                    Tan => a(0).tan(),
                    ASin => a(0).asin(),
                    ACos => a(0).acos(),
                    ATan => a(0).atan(),
                    SinH => a(0).sinh(),
                    CosH => a(0).cosh(),
                    TanH => a(0).tanh(),
                    ASinH => a(0).asinh(),
                    ACosH => a(0).acosh(),
                    ATanH => a(0).atanh(),
                    Exp => a(0).exp(),
                    Expm1 => a(0).exp_m1(),
                    Log => a(0).ln(),
                    Log10 => a(0).log10(),
                    Log2 => a(0).log2(),
                    Sqrt => a(0).sqrt(),
                    Cbrt => a(0).cbrt(),
                    Abs => a(0).abs(),
                    Floor => a(0).floor(),
                    Ceil => a(0).ceil(),
                    Sign => a(0).signum(),
                    Erf => oxistats::erf(a(0)),
                    Erfc => oxistats::erfc(a(0)),
                    LGamma => oxistats::gammaln(a(0)),
                    ATan2 => a(0).atan2(a(1)),
                    Pow => a(0).powf(a(1)),
                    Hypot => a(0).hypot(a(1)),
                    Min => a(0).min(a(1)),
                    Max => a(0).max(a(1)),
                    Pi => std::f64::consts::PI,
                }
            }
            Expr::Cond(c, t, f) => {
                if truthy(c.eval(vars, params)) {
                    t.eval(vars, params)
                } else {
                    f.eval(vars, params)
                }
            }
        }
    }
}

/// Expression builders that fold the zeros and ones a derivative produces, so
/// derivative trees stay small.
fn num(v: f64) -> Expr {
    Expr::Num(v)
}

fn is_num(e: &Expr, v: f64) -> bool {
    matches!(e, Expr::Num(x) if *x == v)
}

fn add(a: Expr, b: Expr) -> Expr {
    if is_num(&a, 0.0) {
        b
    } else if is_num(&b, 0.0) {
        a
    } else {
        Expr::Bin(BinOp::Add, Box::new(a), Box::new(b))
    }
}

fn sub(a: Expr, b: Expr) -> Expr {
    if is_num(&b, 0.0) {
        a
    } else if is_num(&a, 0.0) {
        neg(b)
    } else {
        Expr::Bin(BinOp::Sub, Box::new(a), Box::new(b))
    }
}

fn mul(a: Expr, b: Expr) -> Expr {
    if is_num(&a, 0.0) || is_num(&b, 0.0) {
        num(0.0)
    } else if is_num(&a, 1.0) {
        b
    } else if is_num(&b, 1.0) {
        a
    } else {
        Expr::Bin(BinOp::Mul, Box::new(a), Box::new(b))
    }
}

fn div(a: Expr, b: Expr) -> Expr {
    if is_num(&a, 0.0) {
        num(0.0)
    } else if is_num(&b, 1.0) {
        a
    } else {
        Expr::Bin(BinOp::Div, Box::new(a), Box::new(b))
    }
}

fn neg(a: Expr) -> Expr {
    match a {
        Expr::Num(v) => num(-v),
        other => Expr::Neg(Box::new(other)),
    }
}

fn call(f: Func, args: Vec<Expr>) -> Expr {
    Expr::Call(f, args)
}

fn pow(a: Expr, b: Expr) -> Expr {
    Expr::Bin(BinOp::Pow, Box::new(a), Box::new(b))
}

impl Expr {
    /// Whether the expression depends on parameter `k`.
    fn has_par(&self, k: usize) -> bool {
        match self {
            Expr::Num(_) | Expr::Var(_) => false,
            Expr::Par(i) => *i == k,
            Expr::Neg(e) => e.has_par(k),
            Expr::Bin(_, a, b) => a.has_par(k) || b.has_par(k),
            Expr::Call(_, args) => args.iter().any(|a| a.has_par(k)),
            Expr::Cond(c, t, f) => c.has_par(k) || t.has_par(k) || f.has_par(k),
        }
    }

    /// The derivative with respect to parameter `k`, symbolically; `None` for
    /// a function with no closed-form derivative here (`lgamma`). Steps
    /// (`floor`, `ceil`, `sign`, comparisons) have derivative zero, as almost
    /// everywhere; a condition differentiates each branch.
    pub(crate) fn diff(&self, k: usize) -> Option<Expr> {
        if !self.has_par(k) {
            return Some(num(0.0));
        }
        Some(match self {
            Expr::Num(_) | Expr::Var(_) => num(0.0),
            Expr::Par(i) => num(if *i == k { 1.0 } else { 0.0 }),
            Expr::Neg(e) => neg(e.diff(k)?),
            Expr::Bin(op, a, b) => {
                let (a, b) = (a.as_ref().clone(), b.as_ref().clone());
                let (da, db) = (a.diff(k)?, b.diff(k)?);
                match op {
                    BinOp::Add => add(da, db),
                    BinOp::Sub => sub(da, db),
                    BinOp::Mul => add(mul(da, b), mul(a, db)),
                    BinOp::Div => sub(div(da, b.clone()), div(mul(a, db), mul(b.clone(), b))),
                    BinOp::Pow => {
                        if is_num(&db, 0.0) {
                            // d(aᵇ) = b·aᵇ⁻¹·a' for a constant exponent (any sign of a).
                            mul(mul(b.clone(), pow(a, sub(b, num(1.0)))), da)
                        } else {
                            // d(aᵇ) = aᵇ·(b'·ln a + b·a'/a).
                            mul(
                                pow(a.clone(), b.clone()),
                                add(
                                    mul(db, call(Func::Log, vec![a.clone()])),
                                    div(mul(b, da), a),
                                ),
                            )
                        }
                    }
                    _ => num(0.0), // comparisons and logic: steps
                }
            }
            Expr::Call(f, args) => {
                let a = args.first().cloned().unwrap_or(num(0.0));
                let da = a.diff(k)?;
                use Func::*;
                let unary = |g: Expr| mul(g, da.clone());
                match f {
                    Sin => unary(call(Cos, vec![a])),
                    Cos => unary(neg(call(Sin, vec![a]))),
                    Tan => unary(div(num(1.0), pow(call(Cos, vec![a]), num(2.0)))),
                    ASin => unary(div(
                        num(1.0),
                        call(Sqrt, vec![sub(num(1.0), pow(a, num(2.0)))]),
                    )),
                    ACos => unary(neg(div(
                        num(1.0),
                        call(Sqrt, vec![sub(num(1.0), pow(a, num(2.0)))]),
                    ))),
                    ATan => unary(div(num(1.0), add(num(1.0), pow(a, num(2.0))))),
                    SinH => unary(call(CosH, vec![a])),
                    CosH => unary(call(SinH, vec![a])),
                    TanH => unary(sub(num(1.0), pow(call(TanH, vec![a]), num(2.0)))),
                    ASinH => unary(div(
                        num(1.0),
                        call(Sqrt, vec![add(pow(a, num(2.0)), num(1.0))]),
                    )),
                    ACosH => unary(div(
                        num(1.0),
                        call(Sqrt, vec![sub(pow(a, num(2.0)), num(1.0))]),
                    )),
                    ATanH => unary(div(num(1.0), sub(num(1.0), pow(a, num(2.0))))),
                    Exp | Expm1 => unary(call(Exp, vec![a])),
                    Log => unary(div(num(1.0), a)),
                    Log10 => unary(div(num(1.0 / std::f64::consts::LN_10), a)),
                    Log2 => unary(div(num(1.0 / std::f64::consts::LN_2), a)),
                    Sqrt => unary(div(num(0.5), call(Sqrt, vec![a]))),
                    Cbrt => unary(div(num(1.0 / 3.0), pow(call(Cbrt, vec![a]), num(2.0)))),
                    Abs => unary(call(Sign, vec![a])),
                    Floor | Ceil | Sign | Pi => num(0.0),
                    Erf => unary(mul(
                        num(2.0 / std::f64::consts::PI.sqrt()),
                        call(Exp, vec![neg(pow(a, num(2.0)))]),
                    )),
                    Erfc => unary(mul(
                        num(-2.0 / std::f64::consts::PI.sqrt()),
                        call(Exp, vec![neg(pow(a, num(2.0)))]),
                    )),
                    LGamma => return None,
                    ATan2 => {
                        // atan2(y, x)' = (x·y' − y·x')/(x² + y²).
                        let (y, x) = (a, args[1].clone());
                        let dx = x.diff(k)?;
                        div(
                            sub(mul(x.clone(), da), mul(y.clone(), dx)),
                            add(pow(x, num(2.0)), pow(y, num(2.0))),
                        )
                    }
                    Pow => Expr::Bin(BinOp::Pow, Box::new(a), Box::new(args[1].clone())).diff(k)?,
                    Hypot => {
                        let b = args[1].clone();
                        let db = b.diff(k)?;
                        div(
                            add(mul(a.clone(), da), mul(b.clone(), db)),
                            call(Hypot, vec![a, b]),
                        )
                    }
                    Min | Max => {
                        let b = args[1].clone();
                        let db = b.diff(k)?;
                        let op = if *f == Min { BinOp::Le } else { BinOp::Ge };
                        Expr::Cond(
                            Box::new(Expr::Bin(op, Box::new(a), Box::new(b))),
                            Box::new(da),
                            Box::new(db),
                        )
                    }
                }
            }
            Expr::Cond(c, t, f) => {
                Expr::Cond(c.clone(), Box::new(t.diff(k)?), Box::new(f.diff(k)?))
            }
        })
    }
}
