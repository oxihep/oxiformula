//! The expression tree a formula parses to, and its evaluator.

use crate::func::{Function, MAX_ARGS};

/// A formula's expression tree. Variables are numbered `x` = 0, `y` = 1,
/// `z` = 2, `t` = 3; parameters by their index in
/// [`Formula::param_names`](crate::Formula::param_names). ROOT's shortcuts
/// (`gaus`, `pol2`, …) are expanded, as ROOT expands them.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Expr {
    /// A number.
    Num(f64),
    /// A named constant (`pi`, `e`, `TMath::C()`, …).
    Const(Constant),
    /// A variable: 0 is `x`, 1 `y`, 2 `z`, 3 `t`.
    Var(usize),
    /// A parameter, by index.
    Param(usize),
    /// `-a`.
    Neg(Box<Expr>),
    /// `!a`: 1 where `a` is 0, else 0.
    Not(Box<Expr>),
    /// `a op b`.
    Binary(BinaryOp, Box<Expr>, Box<Expr>),
    /// A function call.
    Call(Function, Vec<Expr>),
    /// `condition ? then : otherwise`.
    Cond(Box<Expr>, Box<Expr>, Box<Expr>),
}

/// A binary operator. Comparisons and logic give `1` or `0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum BinaryOp {
    /// `+`.
    Add,
    /// `-`.
    Sub,
    /// `*`.
    Mul,
    /// `/`.
    Div,
    /// `^` (or `**`), right-associative and tighter than a unary minus.
    Pow,
    /// `<`.
    Lt,
    /// `>`.
    Gt,
    /// `<=`.
    Le,
    /// `>=`.
    Ge,
    /// `==`.
    Eq,
    /// `!=`.
    Ne,
    /// `&&`.
    And,
    /// `||`.
    Or,
}

/// A named constant, with ROOT's value (`TMath::C()`, …).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Constant {
    /// `pi`, `TMath::Pi()`.
    Pi,
    /// `e`, `TMath::E()`.
    E,
    /// `sqrt2`, `TMath::Sqrt2()`.
    Sqrt2,
    /// `infinity`, `TMath::Infinity()`.
    Infinity,
    /// `ln10`, `TMath::Ln10()`.
    Ln10,
    /// `loge`, `TMath::LogE()`: log₁₀ e.
    LogE,
    /// `c`, `TMath::C()`: the speed of light, m/s.
    C,
    /// `g`, `TMath::G()`: the gravitational constant, m³/(kg·s²).
    G,
    /// `h`, `TMath::H()`: Planck's constant, J·s.
    H,
    /// `k`, `TMath::K()`: Boltzmann's constant, J/K.
    K,
    /// `sigma`, `TMath::Sigma()`: the Stefan–Boltzmann constant, W/(m²·K⁴).
    Sigma,
    /// `r`, `TMath::R()`: the gas constant, J/(K·mol).
    R,
    /// `eg`, `TMath::EulerGamma()`.
    EulerGamma,
    /// `true` (1).
    True,
    /// `false` (0).
    False,
}

/// Every constant: `(bare name, TMath function, constant)`.
pub(crate) const CONSTANTS: &[(&str, &str, Constant)] = &[
    ("pi", "TMath::Pi", Constant::Pi),
    ("e", "TMath::E", Constant::E),
    ("sqrt2", "TMath::Sqrt2", Constant::Sqrt2),
    ("infinity", "TMath::Infinity", Constant::Infinity),
    ("ln10", "TMath::Ln10", Constant::Ln10),
    ("loge", "TMath::LogE", Constant::LogE),
    ("c", "TMath::C", Constant::C),
    ("g", "TMath::G", Constant::G),
    ("h", "TMath::H", Constant::H),
    ("k", "TMath::K", Constant::K),
    ("sigma", "TMath::Sigma", Constant::Sigma),
    ("r", "TMath::R", Constant::R),
    ("eg", "TMath::EulerGamma", Constant::EulerGamma),
    ("true", "", Constant::True),
    ("false", "", Constant::False),
];

impl Constant {
    /// Its value, as ROOT 6.40's `TMath` gives it.
    #[must_use]
    pub fn value(self) -> f64 {
        match self {
            Constant::Pi => std::f64::consts::PI,
            Constant::E => std::f64::consts::E,
            Constant::Sqrt2 => std::f64::consts::SQRT_2,
            Constant::Infinity => f64::INFINITY,
            Constant::Ln10 => std::f64::consts::LN_10,
            Constant::LogE => std::f64::consts::LOG10_E,
            Constant::C => 299_792_458.0,
            Constant::G => 6.6743e-11,
            Constant::H => 6.626_070_15e-34,
            Constant::K => 1.380_649e-23,
            Constant::Sigma => 5.670_373e-8,
            Constant::R => 8.314_462_618_153_24,
            Constant::EulerGamma => std::f64::consts::EULER_GAMMA,
            Constant::True => 1.0,
            Constant::False => 0.0,
        }
    }

    /// The name it is written with in a ROOT formula (`pi`, `e`, …).
    #[must_use]
    pub fn root_name(self) -> &'static str {
        CONSTANTS
            .iter()
            .find(|(_, _, c)| *c == self)
            .map_or("?", |(n, _, _)| n)
    }
}

#[inline]
fn truthy(v: f64) -> bool {
    v != 0.0 && !v.is_nan()
}

#[inline]
fn bool_value(b: bool) -> f64 {
    f64::from(u8::from(b))
}

/// `l ^ r`, with an integer exponent done by repeated multiplication (exact,
/// and faster than `powf`).
#[inline]
pub(crate) fn pow(l: f64, r: f64) -> f64 {
    if r.fract() == 0.0 && r.abs() <= 64.0 {
        #[allow(
            clippy::cast_possible_truncation,
            reason = "a whole number of magnitude at most 64"
        )]
        let n = r as i32;
        l.powi(n)
    } else {
        l.powf(r)
    }
}

impl BinaryOp {
    #[inline]
    pub(crate) fn apply(self, l: f64, r: f64) -> f64 {
        match self {
            BinaryOp::Add => l + r,
            BinaryOp::Sub => l - r,
            BinaryOp::Mul => l * r,
            BinaryOp::Div => l / r,
            BinaryOp::Pow => pow(l, r),
            BinaryOp::Lt => bool_value(l < r),
            BinaryOp::Gt => bool_value(l > r),
            BinaryOp::Le => bool_value(l <= r),
            BinaryOp::Ge => bool_value(l >= r),
            BinaryOp::Eq => bool_value(l == r),
            BinaryOp::Ne => bool_value(l != r),
            BinaryOp::And => bool_value(truthy(l) && truthy(r)),
            BinaryOp::Or => bool_value(truthy(l) || truthy(r)),
        }
    }
}

impl Expr {
    /// Evaluate the tree. The caller checks that `vars` and `params` are long
    /// enough; an index past either reads as `NaN`.
    pub(crate) fn eval(&self, vars: &[f64], params: &[f64]) -> f64 {
        match self {
            Expr::Num(v) => *v,
            Expr::Const(c) => c.value(),
            Expr::Var(i) => vars.get(*i).copied().unwrap_or(f64::NAN),
            Expr::Param(i) => params.get(*i).copied().unwrap_or(f64::NAN),
            Expr::Neg(e) => -e.eval(vars, params),
            Expr::Not(e) => bool_value(!truthy(e.eval(vars, params))),
            Expr::Binary(op, a, b) => match op {
                // `&&` and `||` short-circuit, as in C++.
                BinaryOp::And => {
                    bool_value(truthy(a.eval(vars, params)) && truthy(b.eval(vars, params)))
                }
                BinaryOp::Or => {
                    bool_value(truthy(a.eval(vars, params)) || truthy(b.eval(vars, params)))
                }
                _ => op.apply(a.eval(vars, params), b.eval(vars, params)),
            },
            Expr::Call(f, args) => {
                let mut buf = [0.0; MAX_ARGS];
                for (slot, arg) in buf.iter_mut().zip(args) {
                    *slot = arg.eval(vars, params);
                }
                f.eval(&buf[..args.len().min(MAX_ARGS)])
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

    /// Call `visit` on this node and every node below it, parents first.
    pub(crate) fn walk(&self, visit: &mut impl FnMut(&Expr)) {
        visit(self);
        match self {
            Expr::Num(_) | Expr::Const(_) | Expr::Var(_) | Expr::Param(_) => {}
            Expr::Neg(e) | Expr::Not(e) => e.walk(visit),
            Expr::Binary(_, a, b) => {
                a.walk(visit);
                b.walk(visit);
            }
            Expr::Call(_, args) => args.iter().for_each(|a| a.walk(visit)),
            Expr::Cond(c, t, f) => {
                c.walk(visit);
                t.walk(visit);
                f.walk(visit);
            }
        }
    }
}
