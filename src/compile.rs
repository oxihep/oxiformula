//! A formula compiled to nested closures: no per-node `match` at evaluation
//! time, several times faster than walking the tree.

use crate::expr::{BinaryOp, Expr, pow};
use crate::func::{Function, MAX_ARGS};

type Node = Box<dyn Fn(&[f64], &[f64]) -> f64 + Send + Sync>;

/// A [`Formula`](crate::Formula) compiled for evaluation: what to keep when
/// one formula is evaluated many times, as in a fit. Made by
/// [`Formula::compile`](crate::Formula::compile).
///
/// ```
/// # use oxifit_formula::Formula;
/// let f = Formula::parse("gaus")?.compile();
/// let p = [2.0, 0.0, 1.0];
/// assert_eq!(f.eval(&[0.0], &p), 2.0);
/// let mut out = [0.0; 3];
/// f.eval_many(&[-1.0, 0.0, 1.0], &p, &mut out);
/// assert_eq!(out[1], 2.0);
/// # Ok::<(), oxifit_formula::ParseError>(())
/// ```
pub struct Compiled {
    node: Node,
    ndim: usize,
    npar: usize,
}

impl std::fmt::Debug for Compiled {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Compiled")
            .field("ndim", &self.ndim)
            .field("npar", &self.npar)
            .finish_non_exhaustive()
    }
}

impl Compiled {
    pub(crate) fn new(expr: &Expr, ndim: usize, npar: usize) -> Compiled {
        Compiled {
            node: compile(expr),
            ndim,
            npar,
        }
    }

    /// Evaluate at `vars` (`x`, or `x, y`, …) with `params`, as
    /// [`Formula::eval`](crate::Formula::eval) does: `NaN` if either holds
    /// fewer values than the formula uses.
    #[must_use]
    #[inline]
    pub fn eval(&self, vars: &[f64], params: &[f64]) -> f64 {
        if vars.len() < self.ndim || params.len() < self.npar {
            return f64::NAN;
        }
        (self.node)(vars, params)
    }

    /// Evaluate a one-variable formula at each of `xs`, into `out` (as many
    /// values as `xs`; extra entries of `out` are left alone). A formula of
    /// more variables, or too few `params`, gives `NaN` everywhere.
    pub fn eval_many(&self, xs: &[f64], params: &[f64], out: &mut [f64]) {
        let ok = self.ndim == 1 && params.len() >= self.npar;
        for (o, x) in out.iter_mut().zip(xs) {
            *o = if ok {
                (self.node)(std::slice::from_ref(x), params)
            } else {
                f64::NAN
            };
        }
    }

    /// The number of variables the formula uses.
    #[must_use]
    pub fn ndim(&self) -> usize {
        self.ndim
    }

    /// The number of parameters the formula takes.
    #[must_use]
    pub fn npar(&self) -> usize {
        self.npar
    }
}

fn truthy(v: f64) -> bool {
    v != 0.0 && !v.is_nan()
}

fn compile(e: &Expr) -> Node {
    match e {
        Expr::Num(v) => {
            let v = *v;
            Box::new(move |_, _| v)
        }
        Expr::Const(c) => {
            let v = c.value();
            Box::new(move |_, _| v)
        }
        Expr::Var(i) => {
            let i = *i;
            Box::new(move |v, _| v[i])
        }
        Expr::Param(i) => {
            let i = *i;
            Box::new(move |_, p| p[i])
        }
        Expr::Neg(a) => {
            let a = compile(a);
            Box::new(move |v, p| -a(v, p))
        }
        Expr::Not(a) => {
            let a = compile(a);
            Box::new(move |v, p| f64::from(u8::from(!truthy(a(v, p)))))
        }
        Expr::Binary(op, a, b) => binary(*op, a, b),
        Expr::Call(f, args) => call(*f, args),
        Expr::Cond(c, t, f) => {
            let (c, t, f) = (compile(c), compile(t), compile(f));
            Box::new(move |v, p| if truthy(c(v, p)) { t(v, p) } else { f(v, p) })
        }
    }
}

fn binary(op: BinaryOp, a: &Expr, b: &Expr) -> Node {
    // The common powers, without `powi`'s call.
    if op == BinaryOp::Pow
        && let Expr::Num(n) = b
    {
        let (a, n) = (compile(a), *n);
        return if n == 2.0 {
            Box::new(move |v, p| {
                let t = a(v, p);
                t * t
            })
        } else if n == 3.0 {
            Box::new(move |v, p| {
                let t = a(v, p);
                t * t * t
            })
        } else {
            Box::new(move |v, p| pow(a(v, p), n))
        };
    }
    // `t*t`, as ROOT writes a square (its `gaus` is `[0]*exp(-0.5*u*u)`):
    // `t` once.
    if op == BinaryOp::Mul && a == b {
        let a = compile(a);
        return Box::new(move |v, p| {
            let t = a(v, p);
            t * t
        });
    }
    if matches!(
        op,
        BinaryOp::Add | BinaryOp::Sub | BinaryOp::Mul | BinaryOp::Div
    ) {
        let (a, b) = (operand(a), operand(b));
        return match op {
            BinaryOp::Add => Box::new(move |v, p| a.get(v, p) + b.get(v, p)),
            BinaryOp::Sub => Box::new(move |v, p| a.get(v, p) - b.get(v, p)),
            BinaryOp::Mul => Box::new(move |v, p| a.get(v, p) * b.get(v, p)),
            _ => Box::new(move |v, p| a.get(v, p) / b.get(v, p)),
        };
    }
    let (a, b) = (compile(a), compile(b));
    match op {
        BinaryOp::And => {
            Box::new(move |v, p| f64::from(u8::from(truthy(a(v, p)) && truthy(b(v, p)))))
        }
        BinaryOp::Or => {
            Box::new(move |v, p| f64::from(u8::from(truthy(a(v, p)) || truthy(b(v, p)))))
        }
        op => Box::new(move |v, p| op.apply(a(v, p), b(v, p))),
    }
}

/// An operand of an arithmetic operator: a number, variable or parameter
/// read in place, or a compiled subtree. Reading a leaf in place saves the
/// call through a closure for most operands.
enum Operand {
    Num(f64),
    Var(usize),
    Param(usize),
    Node(Node),
}

impl Operand {
    #[inline(always)]
    fn get(&self, v: &[f64], p: &[f64]) -> f64 {
        match self {
            Operand::Num(x) => *x,
            Operand::Var(i) => v[*i],
            Operand::Param(i) => p[*i],
            Operand::Node(n) => n(v, p),
        }
    }
}

fn operand(e: &Expr) -> Operand {
    match e {
        Expr::Num(x) => Operand::Num(*x),
        Expr::Const(c) => Operand::Num(c.value()),
        Expr::Var(i) => Operand::Var(*i),
        Expr::Param(i) => Operand::Param(*i),
        e => Operand::Node(compile(e)),
    }
}

fn call(f: Function, args: &[Expr]) -> Node {
    let mut nodes: Vec<Node> = args.iter().map(compile).collect();
    if nodes.len() == 1 {
        let a = nodes.pop().unwrap_or_else(|| Box::new(|_, _| f64::NAN));
        return match f {
            Function::Exp => Box::new(move |v, p| a(v, p).exp()),
            Function::Sin => Box::new(move |v, p| a(v, p).sin()),
            Function::Cos => Box::new(move |v, p| a(v, p).cos()),
            Function::Sqrt => Box::new(move |v, p| a(v, p).sqrt()),
            Function::Log => Box::new(move |v, p| a(v, p).ln()),
            Function::Abs => Box::new(move |v, p| a(v, p).abs()),
            Function::Sq => Box::new(move |v, p| {
                let t = a(v, p);
                t * t
            }),
            f => Box::new(move |v, p| f.eval(&[a(v, p)])),
        };
    }
    if nodes.len() == 2 {
        let b = nodes.pop().unwrap_or_else(|| Box::new(|_, _| f64::NAN));
        let a = nodes.pop().unwrap_or_else(|| Box::new(|_, _| f64::NAN));
        return Box::new(move |v, p| f.eval(&[a(v, p), b(v, p)]));
    }
    Box::new(move |v, p| {
        let mut buf = [0.0; MAX_ARGS];
        for (slot, n) in buf.iter_mut().zip(&nodes) {
            *slot = n(v, p);
        }
        f.eval(&buf[..nodes.len().min(MAX_ARGS)])
    })
}
