//! A ROOT [`TFormula`](https://root.cern/doc/master/classTFormula.html)
//! engine: parse a formula, evaluate it, differentiate it, and write it back
//! as ROOT compiles it — plus numerical integration and differentiation.
//!
//! It powers oxiroot's `Func1D`/`Func2D`/`Func3D` and oxifit's
//! `Model::from_formula`, so a ROOT formula can be evaluated, fitted, and
//! stored in a file ROOT reads.
//!
//! ```
//! use oxifit_formula::Formula;
//! // Parameters are [0], [1], … or names; variables are x, y, z and t.
//! let f = Formula::parse("[0]*sin([1]*x) + [2]")?;
//! assert_eq!((f.npar(), f.ndim()), (3, 1));
//! let y = f.eval(&[1.0], &[2.0, 1.5, 0.5]); // 2·sin(1.5) + 0.5
//! assert!((y - 2.494_990).abs() < 1e-6);
//!
//! // ROOT's shortcuts, with ROOT's parameter names.
//! let g = Formula::parse("gaus")?;
//! assert_eq!(g.param_names(), ["Constant", "Mean", "Sigma"]);
//!
//! // Written back in ROOT's spelling, whatever was typed.
//! let h = Formula::parse("[0]*ln(x) + sign(x)")?;
//! assert_eq!(h.to_root(), "[p0]*log(x)+TMath::Sign(1.,x)");
//! # Ok::<(), oxifit_formula::ParseError>(())
//! ```
//!
//! # What a formula can hold
//!
//! - **Operators**: `+ - * /`, `^` and `**` (right-associative, tighter
//!   than a unary minus: `-x^2` is `-(x^2)`), comparisons, `&&`, `||`, `!`,
//!   and `cond ? a : b`, with C++ precedence. A comparison gives 1 or 0.
//! - **Variables**: `x`, `y`, `z`, `t`, or ROOT's `x[0]` … `x[3]`.
//! - **Parameters**: `[0]`, `[p1]` or a name, `[mean]`, numbered as ROOT
//!   numbers them: `p<N>` first by number, then the names in byte order, and
//!   `[N]` defines every `p0` … `pN`. So `[b]*x+[a]` has `a` first.
//! - **Constants**: ROOT's `pi`, `e`, `sqrt2`, `infinity`, `ln10`, `loge`,
//!   `c`, `g`, `h`, `k`, `sigma`, `r`, `eg`, `true`, `false`, and the
//!   `TMath::Pi()` forms.
//! - **Functions**: the C library's and `TMath`'s (`sin`, `TMath::Sin`,
//!   `std::sin`, …, `TMath::Sq`, `TMath::Min`, `TMath::Sign`), and
//!   `TMath::Gaus`, `TMath::Landau`, `TMath::BreitWigner`,
//!   `ROOT::Math::normal_pdf`, `breitwigner_pdf`, `crystalball_function`,
//!   `bigaussian_pdf` and `Chebyshev0` … `Chebyshev10` — see [`Function`].
//! - **Shortcuts**, expanded as ROOT expands them: `gaus`, `gausn`, `expo`,
//!   `landau`, `landaun`, `crystalball`, `breitwigner`, `bigaus`, `chebN`,
//!   `polN`; in other variables (`ygaus`, `xygaus`, `gaus[x,y]`, `ypol2`,
//!   `pol2(y, 0)`), from a later parameter (`gaus(3)`), or with named
//!   parameters (`pol2(x, [a], [b], [c])`). A formula that is one shortcut
//!   alone gets ROOT's names for its parameters (`Constant`, `Mean`, …).
//!
//! [`Formula::parse`] also takes a few spellings ROOT does not (`ln`,
//! `sign(x)`, `min`, `power`, any case of a function name), and writes them
//! back in ROOT's: [`to_root`](Formula::to_root) is always something ROOT
//! compiles. [`Formula::parse_root`] refuses them instead, with the ROOT
//! spelling as a hint.

mod compile;
mod diff;
mod error;
mod expr;
mod func;
mod lexer;
mod numeric;
mod parser;
mod render;

pub use compile::Compiled;
pub use error::ParseError;
pub use expr::{BinaryOp, Constant, Expr};
pub use func::Function;
pub use numeric::{derivative, integrate, integrate_adaptive};

use diff::Wrt;
use render::Style;

/// A parsed formula: its expression tree, its parameters' names (in ROOT's
/// order, which is their index), and the number of variables it uses.
///
/// Two formulas are equal when their trees, parameter names and dimensions
/// are: `[0]*x` equals `[p0]*x`, but not `[1]*x+[0]` written the other way
/// round.
#[derive(Debug, Clone)]
pub struct Formula {
    expr: Expr,
    params: Vec<String>,
    ndim: usize,
    source: String,
}

impl Formula {
    /// Parse a formula in ROOT's syntax. A few spellings ROOT itself does not
    /// accept (`ln`, a one-argument `sign`, `min`/`max`, `power`, `pi()`, any
    /// case of a function name) are accepted too, and written back in ROOT's
    /// by [`to_root`](Self::to_root).
    ///
    /// # Errors
    /// A [`ParseError`] that points at what is wrong, if `source` is not a
    /// formula.
    pub fn parse(source: &str) -> Result<Formula, ParseError> {
        Formula::parse_with(source, false)
    }

    /// Parse a formula exactly as ROOT's `TFormula` would accept it: the
    /// spellings [`parse`](Self::parse) takes beyond ROOT's are errors, with
    /// ROOT's spelling as a hint.
    ///
    /// ```
    /// # use oxifit_formula::Formula;
    /// let e = Formula::parse_root("ln(x)").unwrap_err();
    /// assert_eq!(e.hint(), Some("use `log`"));
    /// ```
    ///
    /// # Errors
    /// A [`ParseError`] if `source` is not a formula ROOT accepts.
    pub fn parse_root(source: &str) -> Result<Formula, ParseError> {
        Formula::parse_with(source, true)
    }

    fn parse_with(source: &str, strict: bool) -> Result<Formula, ParseError> {
        let parsed = parser::parse(source, strict).map_err(|e| e.in_source(source))?;
        Ok(Formula {
            expr: parsed.expr,
            params: parsed.params,
            ndim: parsed.ndim,
            source: source.to_owned(),
        })
    }

    /// Evaluate at `vars` (`x`, or `x, y`, …) with `params`, by index.
    /// `NaN` if `vars` has fewer values than the formula has variables, or
    /// `params` fewer than it has parameters.
    #[must_use]
    pub fn eval(&self, vars: &[f64], params: &[f64]) -> f64 {
        if vars.len() < self.ndim || params.len() < self.params.len() {
            return f64::NAN;
        }
        self.expr.eval(vars, params)
    }

    /// Evaluate a formula of `x` alone at `x`: `eval(&[x], params)`, so
    /// `NaN` for a formula that also uses `y`.
    #[must_use]
    pub fn eval1(&self, x: f64, params: &[f64]) -> f64 {
        self.eval(&[x], params)
    }

    /// The formula compiled for repeated evaluation (several times faster
    /// than [`eval`](Self::eval)).
    #[must_use]
    pub fn compile(&self) -> Compiled {
        Compiled::new(&self.expr, self.ndim, self.params.len())
    }

    /// The derivative with respect to parameter `k`, symbolically: exact
    /// where a finite difference is not, and as cheap as the formula. `None`
    /// if `k` is not a parameter, or the formula depends on it through a
    /// function with no closed-form derivative here (`lgamma`, the special
    /// densities).
    ///
    /// ```
    /// # use oxifit_formula::Formula;
    /// let f = Formula::parse("[0]*exp(-0.5*((x-[1])/[2])^2)")?;
    /// let d = f.derivative(1).unwrap(); // ∂/∂mean
    /// let (x, p) = (1.3, [2.0, 1.0, 0.5]);
    /// let u: f64 = (x - p[1]) / p[2];
    /// let exact = p[0] * (-0.5 * u * u).exp() * u / p[2];
    /// assert!((d.eval1(x, &p) - exact).abs() < 1e-12);
    /// # Ok::<(), oxifit_formula::ParseError>(())
    /// ```
    #[must_use]
    pub fn derivative(&self, k: usize) -> Option<Formula> {
        if k >= self.params.len() {
            return None;
        }
        let name = self.params[k].clone();
        self.derived(Wrt::Param(k), &format!("[{name}]"))
    }

    /// The derivative with respect to variable `v` (0 is `x`, 1 `y`, …),
    /// symbolically; `None` as for [`derivative`](Self::derivative).
    ///
    /// ```
    /// # use oxifit_formula::Formula;
    /// let f = Formula::parse("[0]*x^3 + sin(x)")?;
    /// let d = f.variable_derivative(0).unwrap();
    /// assert!((d.eval1(2.0, &[1.5]) - (1.5 * 12.0 + 2f64.cos())).abs() < 1e-12);
    /// # Ok::<(), oxifit_formula::ParseError>(())
    /// ```
    #[must_use]
    pub fn variable_derivative(&self, v: usize) -> Option<Formula> {
        let name = ["x", "y", "z", "t"].get(v)?;
        self.derived(Wrt::Var(v), name)
    }

    fn derived(&self, wrt: Wrt, name: &str) -> Option<Formula> {
        let expr = self.expr.diff(wrt)?;
        Some(Formula {
            expr,
            params: self.params.clone(),
            ndim: self.ndim,
            source: format!("d({})/d{name}", self.source),
        })
    }

    /// The expression tree.
    #[must_use]
    pub fn expr(&self) -> &Expr {
        &self.expr
    }

    /// The number of parameters.
    #[must_use]
    pub fn npar(&self) -> usize {
        self.params.len()
    }

    /// The number of variables: 1 for a formula in `x` (or in none), 2 in
    /// `x, y`, 3 in `x, y, z`, 4 with `t`.
    #[must_use]
    pub fn ndim(&self) -> usize {
        self.ndim
    }

    /// The parameter names, in index order: `p0`, `p1`, … for numbered
    /// parameters, the names of named ones, and ROOT's names for a shortcut
    /// alone (`Constant`, `Mean`, `Sigma` for `gaus`).
    #[must_use]
    pub fn param_names(&self) -> &[String] {
        &self.params
    }

    /// The index of the parameter named `name` (`p2` and `2` both find
    /// parameter `p2`).
    #[must_use]
    pub fn param_index(&self, name: &str) -> Option<usize> {
        let numbered = (!name.is_empty() && name.bytes().all(|b| b.is_ascii_digit())).then(|| {
            format!(
                "p{}",
                name.trim_start_matches('0').parse::<usize>().unwrap_or(0)
            )
        });
        let name = numbered.as_deref().unwrap_or(name);
        self.params.iter().position(|p| p == name)
    }

    /// The same formula with its parameters renamed, in index order: for a
    /// formula read from a file, whose names the file stores apart from it.
    ///
    /// # Errors
    /// A [`ParseError`] if there is not one name per parameter, or a name is
    /// empty, repeated, or not letters, digits and `_`.
    pub fn with_param_names<S: AsRef<str>>(mut self, names: &[S]) -> Result<Formula, ParseError> {
        let names: Vec<String> = names.iter().map(|n| n.as_ref().to_owned()).collect();
        if names.len() != self.params.len() {
            return Err(ParseError::new(
                format!(
                    "{} parameter names for a formula of {} parameters",
                    names.len(),
                    self.params.len()
                ),
                0..0,
            ));
        }
        for (i, n) in names.iter().enumerate() {
            if n.is_empty() || !n.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_') {
                return Err(ParseError::new(
                    format!("`{n}` is not a parameter name"),
                    0..0,
                ));
            }
            if names[..i].contains(n) {
                return Err(ParseError::new(
                    format!("parameter `{n}` is named twice"),
                    0..0,
                ));
            }
        }
        self.params = names;
        Ok(self)
    }

    /// The variables the formula uses, by name (`["x", "y"]`).
    #[must_use]
    pub fn variables(&self) -> Vec<&'static str> {
        let mut used = [false; 4];
        self.expr.walk(&mut |e| {
            if let Expr::Var(i) = e
                && let Some(u) = used.get_mut(*i)
            {
                *u = true;
            }
        });
        ["x", "y", "z", "t"]
            .into_iter()
            .zip(used)
            .filter_map(|(n, u)| u.then_some(n))
            .collect()
    }

    /// The functions the formula calls, once each, in order of appearance.
    #[must_use]
    pub fn functions(&self) -> Vec<Function> {
        let mut out = Vec::new();
        self.expr.walk(&mut |e| {
            if let Expr::Call(f, _) = e
                && !out.contains(f)
            {
                out.push(*f);
            }
        });
        out
    }

    /// The formula as it was passed to [`parse`](Self::parse) — what ROOT
    /// stores as a function's title (`fTitle`).
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    /// The formula as ROOT's `TFormula` compiles it: parameters by index, as
    /// ROOT writes them (`[p2]`, or `[2]` when an earlier parameter is
    /// unused), `pow(a,b)`, ROOT's names for functions and constants, and
    /// whole numbers passed to a function as doubles (`1.`). Parsing it again
    /// gives the same formula, up to parameter names; it is what a file
    /// stores as `TFormula::fFormula`, with the names beside it.
    ///
    /// ```
    /// # use oxifit_formula::Formula;
    /// let f = Formula::parse("[N]*exp(-0.5*((x-[mu])/[sigma])^2)")?;
    /// assert_eq!(f.param_names(), ["N", "mu", "sigma"]);
    /// assert_eq!(f.to_root(), "[p0]*exp(-0.5*pow((x-[p1])/[p2],2.))");
    /// # Ok::<(), oxifit_formula::ParseError>(())
    /// ```
    #[must_use]
    pub fn to_root(&self) -> String {
        // ROOT's own spelling, `[p2]`, when every parameter appears; `[2]`
        // otherwise, which defines the ones before it too.
        let mut used = vec![false; self.params.len()];
        self.expr.walk(&mut |e| {
            if let Expr::Param(i) = e
                && let Some(u) = used.get_mut(*i)
            {
                *u = true;
            }
        });
        let style = if used.iter().all(|u| *u) {
            Style::RootP
        } else {
            Style::Root
        };
        render::text(&self.expr, &self.params, style)
    }

    /// The formula as a ROOT file stores it (`TFormula::fFormula`): ROOT's
    /// spelling, as [`to_root`](Self::to_root), with every parameter by its
    /// name (`[Constant]`, `[p2]`). It reads back with the same parameter
    /// order only beside its names and indices
    /// ([`param_names`](Self::param_names)), as a file keeps them: on its own,
    /// names are numbered in ROOT's order.
    ///
    /// ```
    /// # use oxifit_formula::Formula;
    /// let f = Formula::parse("crystalball")?;
    /// assert_eq!(
    ///     f.to_root_named(),
    ///     "[Constant]*ROOT::Math::crystalball_function(x,[Alpha],[N],[Sigma],[Mean])"
    /// );
    /// # Ok::<(), oxifit_formula::ParseError>(())
    /// ```
    #[must_use]
    pub fn to_root_named(&self) -> String {
        render::text(&self.expr, &self.params, Style::RootNamed)
    }

    /// The formula in LaTeX math, for a plot label: `p_{0}` for numbered
    /// parameters, `\mathrm{Mean}` for named ones.
    ///
    /// ```
    /// # use oxifit_formula::Formula;
    /// let f = Formula::parse("[0]*exp(-x/[tau])")?;
    /// assert_eq!(f.to_latex(), r"p_{0} \cdot e^{\frac{-x}{\mathrm{tau}}}");
    /// # Ok::<(), oxifit_formula::ParseError>(())
    /// ```
    #[must_use]
    pub fn to_latex(&self) -> String {
        render::latex(&self.expr, &self.params)
    }
}

/// The formula for people to read: named parameters by name (`[Mean]`),
/// numbered ones by index (`[2]`), `^` for powers.
impl std::fmt::Display for Formula {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&render::text(&self.expr, &self.params, Style::Display))
    }
}

impl PartialEq for Formula {
    fn eq(&self, other: &Self) -> bool {
        self.expr == other.expr && self.params == other.params && self.ndim == other.ndim
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrate_stays_bounded_on_nan_and_noise() {
        use std::cell::Cell;
        let calls = Cell::new(0_u64);
        // NaN on part of the range: the error estimate is NaN there, which
        // once made the depth-first bisection run 40 levels deep (2⁴⁰ rules).
        let r = crate::integrate(
            |x| {
                calls.set(calls.get() + 1);
                if x > 0.3 { f64::NAN } else { 1.0 }
            },
            -1.0,
            1.0,
        );
        assert!(
            r.is_nan() && calls.get() < 1000,
            "{r} after {} calls",
            calls.get()
        );
        // A noisy integrand never meets the tolerance: capped at 2000 intervals.
        calls.set(0);
        let r = crate::integrate(
            |x: f64| {
                calls.set(calls.get() + 1);
                (x.to_bits().wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 11) as f64 / (1u64 << 53) as f64
            },
            0.0,
            1.0,
        );
        assert!(
            (0.0..=1.0).contains(&r) && calls.get() <= 15 * (2 * 2000 - 1),
            "{} calls",
            calls.get()
        );
        // A narrow Gaussian — which hung the old depth-first bisection (its
        // halved tolerance fell below the rule's error floor) — to 1e-10.
        let s = 0.01;
        let g =
            |x: f64| (-0.5 * ((x - 0.37) / s).powi(2)).exp() / (s * std::f64::consts::TAU.sqrt());
        assert!((crate::integrate(g, -1.0, 1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn integrate_adaptive_is_bounded_and_finds_a_spike() {
        use std::cell::Cell;
        // A normalized Gaussian of width `s` at an off-centre spot of [-50, 50].
        let calls = Cell::new(0_usize);
        let peak = |s: f64| {
            let calls = &calls;
            move |x: f64| {
                calls.set(calls.get() + 1);
                (-0.5 * ((x - 3.7) / s).powi(2)).exp() / (s * std::f64::consts::TAU.sqrt())
            }
        };
        // Narrow next to the range (a single rule over it would miss it).
        let (r, e) = crate::integrate_adaptive(peak(0.3), -50.0, 50.0, 1e-10, 500);
        assert!((r - 1.0).abs() < 1e-8, "{r} ± {e}");
        assert!(calls.get() <= 15 * (2 * 500 - 1), "{} calls", calls.get());
        // An extreme spike, and a cap far too small: still bounded, finite.
        for (s, cap) in [(1e-6, 500), (0.3, 4)] {
            calls.set(0);
            let (r, e) = crate::integrate_adaptive(peak(s), -50.0, 50.0, 1e-10, cap);
            assert!(calls.get() <= 15 * (2 * cap - 1), "{} calls", calls.get());
            assert!(r.is_finite() && e.is_finite());
        }
        // Reversed limits negate.
        let (a, _) = crate::integrate_adaptive(|x| x * x, 0.0, 3.0, 1e-12, 50);
        let (b, _) = crate::integrate_adaptive(|x| x * x, 3.0, 0.0, 1e-12, 50);
        assert!((a - 9.0).abs() < 1e-12 && (a + b).abs() < 1e-15);
    }

    #[test]
    fn integrate_and_derivative_exact_cases() {
        const EPS: f64 = 1e-9;
        // ∫₀¹ x² dx = 1/3 ; d/dx x³ at 2 = 12
        assert!((integrate(|x| x * x, 0.0, 1.0) - 1.0 / 3.0).abs() < EPS);
        assert!((integrate(|x| x.sin(), 0.0, std::f64::consts::PI) - 2.0).abs() < EPS);
        assert!((derivative(|x| x * x * x, 2.0) - 12.0).abs() < 1e-6);
        // reversed limits negate.
        assert!((integrate(|x| x, 1.0, 0.0) + 0.5).abs() < EPS);
    }
}
