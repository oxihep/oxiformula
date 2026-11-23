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
    /// The formula, with each repeated subexpression read from a slot.
    node: Node,
    /// The repeated subexpressions, in the order they are computed: each into
    /// variable slot `SLOT_BASE + k`, before `node`.
    slots: Vec<(usize, Node)>,
    ndim: usize,
    npar: usize,
}

/// Where the slots of repeated subexpressions start among the variables:
/// after `x`, `y`, `z` and `t`.
const SLOT_BASE: usize = 4;
/// The most repeated subexpressions computed once.
const MAX_SLOTS: usize = 8;
/// Above this many nodes, a formula is compiled as it is (finding repeats
/// costs the square of its size).
const MAX_CSE_NODES: usize = 400;

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
        let (root, defs, order) = hoist(expr);
        Compiled {
            node: compile(&root),
            slots: order
                .into_iter()
                .map(|k| (SLOT_BASE + k, compile(&defs[k])))
                .collect(),
            ndim,
            npar,
        }
    }

    #[inline]
    fn run(&self, vars: &[f64], params: &[f64]) -> f64 {
        if self.slots.is_empty() {
            return (self.node)(vars, params);
        }
        let mut buf = [0.0; SLOT_BASE + MAX_SLOTS];
        let n = vars.len().min(SLOT_BASE);
        buf[..n].copy_from_slice(&vars[..n]);
        for (at, slot) in &self.slots {
            buf[*at] = slot(&buf, params);
        }
        (self.node)(&buf, params)
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
        self.run(vars, params)
    }

    /// Evaluate a one-variable formula at each of `xs`, into `out` (as many
    /// values as `xs`; extra entries of `out` are left alone). A formula of
    /// more variables, or too few `params`, gives `NaN` everywhere.
    pub fn eval_many(&self, xs: &[f64], params: &[f64], out: &mut [f64]) {
        let ok = self.ndim == 1 && params.len() >= self.npar;
        for (o, x) in out.iter_mut().zip(xs) {
            *o = if ok {
                self.run(std::slice::from_ref(x), params)
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

/// Common subexpressions: `expr` with each repeated subtree of four or more
/// nodes replaced by a variable slot (`SLOT_BASE + k`), the subtrees by `k`,
/// and an order to compute them in (each after the slots it reads). The
/// largest repeat is taken first; at most `MAX_SLOTS`.
fn hoist(expr: &Expr) -> (Expr, Vec<Expr>, Vec<usize>) {
    let mut root = expr.clone();
    let mut defs: Vec<Expr> = Vec::new();
    if size(expr) > MAX_CSE_NODES {
        return (root, defs, Vec::new());
    }
    while defs.len() < MAX_SLOTS {
        let mut seen: Vec<(Expr, usize)> = Vec::new();
        count(&root, &mut seen);
        for d in &defs {
            count_children(d, &mut seen);
        }
        let Some(best) = seen
            .into_iter()
            .filter(|(e, n)| *n >= 2 && size(e) >= 4)
            .max_by_key(|(e, _)| size(e))
            .map(|(e, _)| e)
        else {
            break;
        };
        let slot = Expr::Var(SLOT_BASE + defs.len());
        replace(&mut root, &best, &slot);
        for d in &mut defs {
            replace_children(d, &best, &slot);
        }
        defs.push(best);
    }
    // Each slot after the slots it reads.
    let mut order = Vec::with_capacity(defs.len());
    let mut done = vec![false; defs.len()];
    fn visit(k: usize, defs: &[Expr], done: &mut [bool], order: &mut Vec<usize>) {
        if done[k] {
            return;
        }
        done[k] = true;
        defs[k].walk(&mut |e| {
            if let Expr::Var(v) = e
                && *v >= SLOT_BASE
                && *v - SLOT_BASE != k
            {
                visit(*v - SLOT_BASE, defs, done, order);
            }
        });
        order.push(k);
    }
    for k in 0..defs.len() {
        visit(k, &defs, &mut done, &mut order);
    }
    (root, defs, order)
}

fn size(e: &Expr) -> usize {
    let mut n = 0;
    e.walk(&mut |_| n += 1);
    n
}

fn is_leaf(e: &Expr) -> bool {
    matches!(
        e,
        Expr::Num(_) | Expr::Const(_) | Expr::Var(_) | Expr::Param(_)
    )
}

/// Count every subtree of `e` that is not a leaf.
fn count(e: &Expr, seen: &mut Vec<(Expr, usize)>) {
    if is_leaf(e) {
        return;
    }
    match seen.iter_mut().find(|(s, _)| s == e) {
        Some((_, n)) => *n += 1,
        None => seen.push((e.clone(), 1)),
    }
    count_children(e, seen);
}

fn count_children(e: &Expr, seen: &mut Vec<(Expr, usize)>) {
    match e {
        Expr::Neg(a) | Expr::Not(a) => count(a, seen),
        Expr::Binary(_, a, b) => {
            count(a, seen);
            count(b, seen);
        }
        Expr::Call(_, args) => args.iter().for_each(|a| count(a, seen)),
        Expr::Cond(c, t, f) => {
            count(c, seen);
            count(t, seen);
            count(f, seen);
        }
        _ => {}
    }
}

/// Replace every occurrence of `what` in `e` by `with`.
fn replace(e: &mut Expr, what: &Expr, with: &Expr) {
    if e == what {
        *e = with.clone();
    } else {
        replace_children(e, what, with);
    }
}

fn replace_children(e: &mut Expr, what: &Expr, with: &Expr) {
    match e {
        Expr::Neg(a) | Expr::Not(a) => replace(a, what, with),
        Expr::Binary(_, a, b) => {
            replace(a, what, with);
            replace(b, what, with);
        }
        Expr::Call(_, args) => args.iter_mut().for_each(|a| replace(a, what, with)),
        Expr::Cond(c, t, f) => {
            replace(c, what, with);
            replace(t, what, with);
            replace(f, what, with);
        }
        _ => {}
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Formula;

    #[test]
    fn repeated_subexpressions_are_computed_once() {
        // ROOT's `gaus` is `[0]*exp(-0.5*u*u)` with `u = (x-[1])/[2]` twice.
        let g = Formula::parse("gaus").unwrap();
        let (_, defs, order) = hoist(g.expr());
        assert_eq!((defs.len(), order), (1, vec![0]));
        // A repeat inside a larger repeat: the inner slot comes first.
        let f = Formula::parse("exp(sin([0]*x+1)*2) + exp(sin([0]*x+1)*2) + sin([0]*x+1)").unwrap();
        let (_, defs, order) = hoist(f.expr());
        assert_eq!(defs.len(), 2);
        assert_eq!(order, vec![1, 0]);
        let c = f.compile();
        for x in [-1.3, 0.0, 0.4, 2.5] {
            assert_eq!(
                c.eval(&[x], &[0.7]).to_bits(),
                f.eval(&[x], &[0.7]).to_bits()
            );
        }
        // A leaf, or one operation on two leaves, is cheaper to recompute.
        assert!(
            hoist(Formula::parse("x*[0]+x*[0]").unwrap().expr())
                .1
                .is_empty()
        );
    }
}
