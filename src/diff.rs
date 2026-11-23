//! Symbolic derivatives of an expression tree, with respect to a parameter
//! or a variable.

use crate::expr::{BinaryOp, Expr};
use crate::func::Function;

/// What a derivative is taken with respect to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Wrt {
    Param(usize),
    Var(usize),
}

// Builders that fold the zeros and ones a derivative produces, so derivative
// trees stay small.

fn num(v: f64) -> Expr {
    Expr::Num(v)
}

fn is_num(e: &Expr, v: f64) -> bool {
    matches!(e, Expr::Num(x) if *x == v)
}

fn bin(op: BinaryOp, a: Expr, b: Expr) -> Expr {
    Expr::Binary(op, Box::new(a), Box::new(b))
}

fn add(a: Expr, b: Expr) -> Expr {
    if is_num(&a, 0.0) {
        b
    } else if is_num(&b, 0.0) {
        a
    } else {
        bin(BinaryOp::Add, a, b)
    }
}

fn sub(a: Expr, b: Expr) -> Expr {
    if is_num(&b, 0.0) {
        a
    } else if is_num(&a, 0.0) {
        neg(b)
    } else {
        bin(BinaryOp::Sub, a, b)
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
        bin(BinaryOp::Mul, a, b)
    }
}

fn div(a: Expr, b: Expr) -> Expr {
    if is_num(&a, 0.0) {
        num(0.0)
    } else if is_num(&b, 1.0) {
        a
    } else {
        bin(BinaryOp::Div, a, b)
    }
}

fn neg(a: Expr) -> Expr {
    match a {
        Expr::Num(v) => num(-v),
        Expr::Neg(inner) => *inner,
        other => Expr::Neg(Box::new(other)),
    }
}

fn pow(a: Expr, b: Expr) -> Expr {
    bin(BinaryOp::Pow, a, b)
}

fn call(f: Function, args: Vec<Expr>) -> Expr {
    Expr::Call(f, args)
}

fn sq(a: Expr) -> Expr {
    pow(a, num(2.0))
}

impl Expr {
    /// Whether the expression depends on `wrt`.
    pub(crate) fn depends(&self, wrt: Wrt) -> bool {
        match self {
            Expr::Num(_) | Expr::Const(_) => false,
            Expr::Param(i) => wrt == Wrt::Param(*i),
            Expr::Var(i) => wrt == Wrt::Var(*i),
            Expr::Neg(e) | Expr::Not(e) => e.depends(wrt),
            Expr::Binary(_, a, b) => a.depends(wrt) || b.depends(wrt),
            Expr::Call(_, args) => args.iter().any(|a| a.depends(wrt)),
            Expr::Cond(c, t, f) => c.depends(wrt) || t.depends(wrt) || f.depends(wrt),
        }
    }

    /// The derivative with respect to `wrt`, symbolically; `None` where a
    /// function it depends on through has no closed-form derivative here
    /// (`lgamma`, the special densities). Steps (`floor`, `ceil`, the sign,
    /// comparisons) have derivative zero, as almost everywhere; a condition
    /// differentiates each branch.
    pub(crate) fn diff(&self, wrt: Wrt) -> Option<Expr> {
        if !self.depends(wrt) {
            return Some(num(0.0));
        }
        Some(match self {
            Expr::Num(_) | Expr::Const(_) | Expr::Not(_) => num(0.0),
            Expr::Param(_) | Expr::Var(_) => num(1.0),
            Expr::Neg(e) => neg(e.diff(wrt)?),
            Expr::Binary(op, a, b) => {
                let (a, b) = (a.as_ref(), b.as_ref());
                match op {
                    BinaryOp::Add => add(a.diff(wrt)?, b.diff(wrt)?),
                    BinaryOp::Sub => sub(a.diff(wrt)?, b.diff(wrt)?),
                    BinaryOp::Mul => {
                        add(mul(a.diff(wrt)?, b.clone()), mul(a.clone(), b.diff(wrt)?))
                    }
                    BinaryOp::Div => sub(
                        div(a.diff(wrt)?, b.clone()),
                        div(mul(a.clone(), b.diff(wrt)?), sq(b.clone())),
                    ),
                    BinaryOp::Pow => power(a, b, wrt)?,
                    _ => num(0.0), // comparisons and logic: steps
                }
            }
            Expr::Call(f, args) => match f {
                Function::Min | Function::Max => {
                    let (a, b) = (&args[0], &args[1]);
                    let op = if *f == Function::Min {
                        BinaryOp::Le
                    } else {
                        BinaryOp::Ge
                    };
                    Expr::Cond(
                        Box::new(bin(op, a.clone(), b.clone())),
                        Box::new(a.diff(wrt)?),
                        Box::new(b.diff(wrt)?),
                    )
                }
                Function::Pow => power(&args[0], &args[1], wrt)?,
                _ => {
                    let mut total = num(0.0);
                    for (i, a) in args.iter().enumerate() {
                        if a.depends(wrt) {
                            total = add(total, mul(partial(*f, args, i)?, a.diff(wrt)?));
                        }
                    }
                    total
                }
            },
            Expr::Cond(c, t, f) => {
                Expr::Cond(c.clone(), Box::new(t.diff(wrt)?), Box::new(f.diff(wrt)?))
            }
        })
    }
}

/// `d(aᵇ)`.
fn power(a: &Expr, b: &Expr, wrt: Wrt) -> Option<Expr> {
    let da = a.diff(wrt)?;
    if !b.depends(wrt) {
        // b·aᵇ⁻¹·a' for a constant exponent (any sign of a).
        return Some(mul(
            mul(b.clone(), pow(a.clone(), sub(b.clone(), num(1.0)))),
            da,
        ));
    }
    // aᵇ·(b'·ln a + b·a'/a).
    Some(mul(
        pow(a.clone(), b.clone()),
        add(
            mul(b.diff(wrt)?, call(Function::Log, vec![a.clone()])),
            div(mul(b.clone(), da), a.clone()),
        ),
    ))
}

/// `∂f/∂argᵢ` at `args`, or `None` where there is no closed form here.
fn partial(f: Function, args: &[Expr], i: usize) -> Option<Expr> {
    use Function::*;
    let a = args[0].clone();
    let one = || num(1.0);
    Some(match (f, i) {
        (Sin, 0) => call(Cos, vec![a]),
        (Cos, 0) => neg(call(Sin, vec![a])),
        (Tan, 0) => div(one(), sq(call(Cos, vec![a]))),
        (ASin, 0) => div(one(), call(Sqrt, vec![sub(one(), sq(a))])),
        (ACos, 0) => neg(div(one(), call(Sqrt, vec![sub(one(), sq(a))]))),
        (ATan, 0) => div(one(), add(one(), sq(a))),
        (ATan2, _) => {
            // atan2(y, x): ∂y = x/(x² + y²), ∂x = −y/(x² + y²).
            let (y, x) = (a, args[1].clone());
            let r2 = add(sq(x.clone()), sq(y.clone()));
            if i == 0 { div(x, r2) } else { neg(div(y, r2)) }
        }
        (SinH, 0) => call(CosH, vec![a]),
        (CosH, 0) => call(SinH, vec![a]),
        (TanH, 0) => sub(one(), sq(call(TanH, vec![a]))),
        (ASinH, 0) => div(one(), call(Sqrt, vec![add(sq(a), one())])),
        (ACosH, 0) => div(one(), call(Sqrt, vec![sub(sq(a), one())])),
        (ATanH, 0) => div(one(), sub(one(), sq(a))),
        (Exp | Expm1, 0) => call(Exp, vec![a]),
        (Log, 0) => div(one(), a),
        (Log1p, 0) => div(one(), add(one(), a)),
        (Log10, 0) => div(num(1.0 / std::f64::consts::LN_10), a),
        (Log2, 0) => div(num(1.0 / std::f64::consts::LN_2), a),
        (Sqrt, 0) => div(num(0.5), call(Sqrt, vec![a])),
        (Cbrt, 0) => div(num(1.0 / 3.0), sq(call(Cbrt, vec![a]))),
        (Sq, 0) => mul(num(2.0), a),
        (Abs, 0) => call(Signum, vec![a]),
        (Floor | Ceil | Signum, _) => num(0.0),
        (Erf, 0) => mul(
            num(2.0 / std::f64::consts::PI.sqrt()),
            call(Exp, vec![neg(sq(a))]),
        ),
        (Erfc, 0) => mul(
            num(-2.0 / std::f64::consts::PI.sqrt()),
            call(Exp, vec![neg(sq(a))]),
        ),
        (Hypot, _) => div(args[i].clone(), call(Hypot, args.to_vec())),
        (Fmod, 0) => one(),
        (Fmod, _) => {
            // fmod(a, b) = a − trunc(a/b)·b: ∂b = −trunc(a/b).
            let q = div(a, args[1].clone());
            neg(Expr::Cond(
                Box::new(bin(BinaryOp::Ge, q.clone(), num(0.0))),
                Box::new(call(Floor, vec![q.clone()])),
                Box::new(call(Ceil, vec![q])),
            ))
        }
        (Sign, 0) => mul(
            call(Signum, vec![a]),
            call(Sign, vec![one(), args[1].clone()]),
        ),
        (Sign, _) => num(0.0),
        (Chebyshev(_), k) if k > 0 => {
            // Linear in its coefficients: ∂/∂cₖ is Tₖ(x).
            let n = k - 1;
            let mut coeffs = vec![num(0.0); n];
            coeffs.push(one());
            let mut cargs = vec![a];
            cargs.extend(coeffs);
            call(Chebyshev(u8::try_from(n).ok()?), cargs)
        }
        _ => return None,
    })
}
