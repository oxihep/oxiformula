//! The functions a formula can call: their ROOT spellings, their arity and
//! their values, as ROOT computes them.
//!
//! Every function is one [`Function`] and one row of [`NAMES`], so adding one
//! is local: a variant, its spellings, its arity, its value (and, in
//! `diff.rs`, its derivative).

use std::borrow::Cow;

/// A function a formula calls, resolved from its name when the formula is
/// parsed, so evaluating it is a `match`, never a string lookup.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Function {
    /// `sin`.
    Sin,
    /// `cos`.
    Cos,
    /// `tan`.
    Tan,
    /// `asin`.
    ASin,
    /// `acos`.
    ACos,
    /// `atan`.
    ATan,
    /// `atan2(y, x)`.
    ATan2,
    /// `sinh`.
    SinH,
    /// `cosh`.
    CosH,
    /// `tanh`.
    TanH,
    /// `asinh`.
    ASinH,
    /// `acosh`.
    ACosH,
    /// `atanh`.
    ATanH,
    /// `exp`.
    Exp,
    /// `expm1`: `exp(x) − 1`, exact near 0.
    Expm1,
    /// `log`, the natural logarithm.
    Log,
    /// `log1p`: `log(1 + x)`, exact near 0.
    Log1p,
    /// `log10`.
    Log10,
    /// `log2`.
    Log2,
    /// `sqrt`.
    Sqrt,
    /// `cbrt`.
    Cbrt,
    /// `TMath::Sq`: the square.
    Sq,
    /// `TMath::Abs`.
    Abs,
    /// `floor`.
    Floor,
    /// `ceil`.
    Ceil,
    /// `erf`.
    Erf,
    /// `erfc`.
    Erfc,
    /// `lgamma` (`TMath::LnGamma`).
    LGamma,
    /// `pow(a, b)`.
    Pow,
    /// `hypot(a, b)`.
    Hypot,
    /// `fmod(a, b)`: the remainder of `a / b`, with the sign of `a`.
    Fmod,
    /// `TMath::Min(a, b)`.
    Min,
    /// `TMath::Max(a, b)`.
    Max,
    /// `TMath::Sign(a, b)`: `|a|` with the sign of `b` (`b ≥ 0` is positive).
    Sign,
    /// The sign of one value, `+1` for `x ≥ 0` and `−1` otherwise: what
    /// `sign(x)` means here. ROOT has no one-argument form; it is written to
    /// ROOT as `TMath::Sign(1., x)`.
    Signum,
    /// `TMath::Gaus(x, mean = 0, sigma = 1, norm = false)`.
    Gaus,
    /// `TMath::Landau(x, mpv = 0, sigma = 1, norm = false)`.
    Landau,
    /// `TMath::BreitWigner(x, mean = 0, gamma = 1)`: a unit-area Lorentzian.
    BreitWigner,
    /// `ROOT::Math::breitwigner_pdf(x, gamma, x0 = 0)`.
    BreitWignerPdf,
    /// `ROOT::Math::normal_pdf(x, sigma = 1, x0 = 0)`.
    NormalPdf,
    /// `ROOT::Math::crystalball_function(x, alpha, n, sigma, mean)`: unit
    /// peak, the tail below `mean − alpha·sigma` (above it for `alpha < 0`).
    CrystalBall,
    /// `ROOT::Math::bigaussian_pdf(x, y, sigmax = 1, sigmay = 1, rho = 0,
    /// x0 = 0, y0 = 0)`.
    BiGaussianPdf,
    /// `ROOT::Math::ChebyshevN(x, c0, …, cN)`: `Σ cₖ·Tₖ(x)`, for N from 0 to
    /// 10.
    Chebyshev(u8),
    /// `ROOT::Math::crystalball_pdf(x, alpha, n, sigma, mean = 0)`: the
    /// crystal ball of unit area (`NaN` for `n ≤ 1`).
    CrystalBallPdf,
    /// `ROOT::Math::landau_pdf(x, xi = 1, x0 = 0)`.
    LandauPdf,
    /// `ROOT::Math::beta_pdf(x, a, b)`.
    BetaPdf,
    /// `ROOT::Math::cauchy_pdf(x, b = 1, x0 = 0)`.
    CauchyPdf,
    /// `ROOT::Math::chisquared_pdf(x, r, x0 = 0)`.
    ChiSquaredPdf,
    /// `ROOT::Math::exponential_pdf(x, lambda, x0 = 0)`.
    ExponentialPdf,
    /// `ROOT::Math::fdistribution_pdf(x, n, m, x0 = 0)`.
    FDistributionPdf,
    /// `ROOT::Math::gamma_pdf(x, alpha, theta, x0 = 0)`.
    GammaPdf,
    /// `ROOT::Math::lognormal_pdf(x, m, s, x0 = 0)`.
    LognormalPdf,
    /// `ROOT::Math::poisson_pdf(n, mu)`, `n` truncated to a whole number as
    /// C++ does (`NaN` for a negative `n`).
    PoissonPdf,
    /// `ROOT::Math::tdistribution_pdf(x, r, x0 = 0)`.
    TDistributionPdf,
    /// `ROOT::Math::uniform_pdf(x, a, b, x0 = 0)`.
    UniformPdf,
    /// `TMath::Poisson(x, mu)`: the Poisson density continued to real `x`.
    Poisson,
}

use Function::*;

/// Every name a function is known by: `(name, function, ROOT accepts it)`.
/// The first ROOT name of a function is the one it is written with, except
/// where [`Function::root_name`] says otherwise. Names ROOT rejects are
/// accepted by [`Formula::parse`](crate::Formula::parse) and refused, with
/// the ROOT spelling as a hint, by
/// [`Formula::parse_root`](crate::Formula::parse_root).
pub(crate) const NAMES: &[(&str, Function, bool)] = &[
    ("sin", Sin, true),
    ("TMath::Sin", Sin, true),
    ("std::sin", Sin, true),
    ("cos", Cos, true),
    ("TMath::Cos", Cos, true),
    ("std::cos", Cos, true),
    ("tan", Tan, true),
    ("TMath::Tan", Tan, true),
    ("std::tan", Tan, true),
    ("asin", ASin, true),
    ("TMath::ASin", ASin, true),
    ("std::asin", ASin, true),
    ("acos", ACos, true),
    ("TMath::ACos", ACos, true),
    ("std::acos", ACos, true),
    ("atan", ATan, true),
    ("TMath::ATan", ATan, true),
    ("std::atan", ATan, true),
    ("atan2", ATan2, true),
    ("TMath::ATan2", ATan2, true),
    ("std::atan2", ATan2, true),
    ("sinh", SinH, true),
    ("TMath::SinH", SinH, true),
    ("std::sinh", SinH, true),
    ("cosh", CosH, true),
    ("TMath::CosH", CosH, true),
    ("std::cosh", CosH, true),
    ("tanh", TanH, true),
    ("TMath::TanH", TanH, true),
    ("std::tanh", TanH, true),
    ("asinh", ASinH, true),
    ("TMath::ASinH", ASinH, true),
    ("std::asinh", ASinH, true),
    ("acosh", ACosH, true),
    ("TMath::ACosH", ACosH, true),
    ("std::acosh", ACosH, true),
    ("atanh", ATanH, true),
    ("TMath::ATanH", ATanH, true),
    ("std::atanh", ATanH, true),
    ("exp", Exp, true),
    ("TMath::Exp", Exp, true),
    ("std::exp", Exp, true),
    ("expm1", Expm1, true),
    ("std::expm1", Expm1, true),
    ("log", Log, true),
    ("TMath::Log", Log, true),
    ("std::log", Log, true),
    ("ln", Log, false),
    ("log1p", Log1p, true),
    ("std::log1p", Log1p, true),
    ("log10", Log10, true),
    ("TMath::Log10", Log10, true),
    ("std::log10", Log10, true),
    ("log2", Log2, true),
    ("TMath::Log2", Log2, true),
    ("std::log2", Log2, true),
    ("sqrt", Sqrt, true),
    ("TMath::Sqrt", Sqrt, true),
    ("std::sqrt", Sqrt, true),
    ("cbrt", Cbrt, true),
    ("std::cbrt", Cbrt, true),
    ("TMath::Sq", Sq, true),
    ("sq", Sq, true),
    ("TMath::Abs", Abs, true),
    ("abs", Abs, true),
    ("fabs", Abs, true),
    ("std::abs", Abs, true),
    ("std::fabs", Abs, true),
    ("floor", Floor, true),
    ("TMath::Floor", Floor, true),
    ("std::floor", Floor, true),
    ("ceil", Ceil, true),
    ("TMath::Ceil", Ceil, true),
    ("std::ceil", Ceil, true),
    ("erf", Erf, true),
    ("TMath::Erf", Erf, true),
    ("std::erf", Erf, true),
    ("erfc", Erfc, true),
    ("TMath::Erfc", Erfc, true),
    ("std::erfc", Erfc, true),
    ("lgamma", LGamma, true),
    ("TMath::LnGamma", LGamma, true),
    ("std::lgamma", LGamma, true),
    ("lngamma", LGamma, false),
    ("pow", Pow, true),
    ("TMath::Power", Pow, true),
    ("std::pow", Pow, true),
    ("power", Pow, false),
    ("hypot", Hypot, true),
    ("TMath::Hypot", Hypot, true),
    ("std::hypot", Hypot, true),
    ("fmod", Fmod, true),
    ("std::fmod", Fmod, true),
    ("TMath::Min", Min, true),
    ("std::min", Min, true),
    ("fmin", Min, true),
    ("min", Min, false),
    ("TMath::Max", Max, true),
    ("std::max", Max, true),
    ("fmax", Max, true),
    ("max", Max, false),
    ("TMath::Sign", Sign, true),
    ("sign", Sign, true),
    ("sign", Signum, false),
    ("TMath::Gaus", Gaus, true),
    ("TMath::Landau", Landau, true),
    ("TMath::BreitWigner", BreitWigner, true),
    ("ROOT::Math::breitwigner_pdf", BreitWignerPdf, true),
    ("ROOT::Math::normal_pdf", NormalPdf, true),
    ("ROOT::Math::gaussian_pdf", NormalPdf, true),
    ("ROOT::Math::crystalball_function", CrystalBall, true),
    ("ROOT::Math::bigaussian_pdf", BiGaussianPdf, true),
    ("ROOT::Math::crystalball_pdf", CrystalBallPdf, true),
    ("ROOT::Math::landau_pdf", LandauPdf, true),
    ("ROOT::Math::beta_pdf", BetaPdf, true),
    ("ROOT::Math::cauchy_pdf", CauchyPdf, true),
    ("ROOT::Math::chisquared_pdf", ChiSquaredPdf, true),
    ("ROOT::Math::exponential_pdf", ExponentialPdf, true),
    ("ROOT::Math::fdistribution_pdf", FDistributionPdf, true),
    ("ROOT::Math::gamma_pdf", GammaPdf, true),
    ("ROOT::Math::lognormal_pdf", LognormalPdf, true),
    ("ROOT::Math::poisson_pdf", PoissonPdf, true),
    ("ROOT::Math::tdistribution_pdf", TDistributionPdf, true),
    ("ROOT::Math::uniform_pdf", UniformPdf, true),
    ("TMath::Poisson", Poisson, true),
    ("ROOT::Math::Chebyshev0", Chebyshev(0), true),
    ("ROOT::Math::Chebyshev1", Chebyshev(1), true),
    ("ROOT::Math::Chebyshev2", Chebyshev(2), true),
    ("ROOT::Math::Chebyshev3", Chebyshev(3), true),
    ("ROOT::Math::Chebyshev4", Chebyshev(4), true),
    ("ROOT::Math::Chebyshev5", Chebyshev(5), true),
    ("ROOT::Math::Chebyshev6", Chebyshev(6), true),
    ("ROOT::Math::Chebyshev7", Chebyshev(7), true),
    ("ROOT::Math::Chebyshev8", Chebyshev(8), true),
    ("ROOT::Math::Chebyshev9", Chebyshev(9), true),
    ("ROOT::Math::Chebyshev10", Chebyshev(10), true),
];

/// The most arguments a function takes (`Chebyshev10`: `x` and 11
/// coefficients).
pub(crate) const MAX_ARGS: usize = 12;

impl Function {
    /// The name this function is written with in a ROOT formula.
    #[must_use]
    pub fn root_name(self) -> Cow<'static, str> {
        if let Chebyshev(n) = self {
            return Cow::Owned(format!("ROOT::Math::Chebyshev{n}"));
        }
        if self == Signum {
            return Cow::Borrowed("TMath::Sign");
        }
        let name = NAMES
            .iter()
            .find(|(_, f, root)| *f == self && *root)
            .map_or("?", |(n, _, _)| n);
        Cow::Borrowed(name)
    }

    /// How many arguments the function takes: `(fewest, most)`. Fewer than
    /// the most leaves ROOT's defaults for the rest.
    #[must_use]
    pub fn arity(self) -> (usize, usize) {
        match self {
            ATan2 | Pow | Hypot | Fmod | Min | Max | Sign => (2, 2),
            Gaus | Landau => (1, 4),
            BreitWigner | NormalPdf => (1, 3),
            BreitWignerPdf => (2, 3),
            CrystalBall | CrystalBallPdf => (4, 5),
            LandauPdf | CauchyPdf => (1, 3),
            ChiSquaredPdf | ExponentialPdf | TDistributionPdf => (2, 3),
            BetaPdf => (3, 3),
            FDistributionPdf | GammaPdf | LognormalPdf | UniformPdf => (3, 4),
            PoissonPdf | Poisson => (2, 2),
            BiGaussianPdf => (2, 7),
            Chebyshev(n) => (usize::from(n) + 2, usize::from(n) + 2),
            _ => (1, 1),
        }
    }

    /// The value at `args` (as many as [`arity`](Self::arity) allows).
    #[must_use]
    pub fn eval(self, args: &[f64]) -> f64 {
        let a = |i: usize, default: f64| args.get(i).copied().unwrap_or(default);
        let x = a(0, f64::NAN);
        match self {
            Sin => x.sin(),
            Cos => x.cos(),
            Tan => x.tan(),
            ASin => x.asin(),
            ACos => x.acos(),
            ATan => x.atan(),
            ATan2 => x.atan2(a(1, f64::NAN)),
            SinH => x.sinh(),
            CosH => x.cosh(),
            TanH => x.tanh(),
            ASinH => x.asinh(),
            ACosH => x.acosh(),
            ATanH => x.atanh(),
            Exp => x.exp(),
            Expm1 => x.exp_m1(),
            Log => x.ln(),
            Log1p => x.ln_1p(),
            Log10 => x.log10(),
            Log2 => x.log2(),
            Sqrt => x.sqrt(),
            Cbrt => x.cbrt(),
            Sq => x * x,
            Abs => x.abs(),
            Floor => x.floor(),
            Ceil => x.ceil(),
            Erf => oxistats::erf(x),
            Erfc => oxistats::erfc(x),
            LGamma => oxistats::gammaln(x),
            Pow => x.powf(a(1, f64::NAN)),
            Hypot => x.hypot(a(1, f64::NAN)),
            Fmod => x % a(1, f64::NAN),
            Min => {
                let y = a(1, f64::NAN);
                if x <= y { x } else { y }
            }
            Max => {
                let y = a(1, f64::NAN);
                if x >= y { x } else { y }
            }
            Sign => {
                if a(1, f64::NAN) >= 0.0 {
                    x.abs()
                } else {
                    -x.abs()
                }
            }
            Signum => {
                if x >= 0.0 {
                    1.0
                } else {
                    -1.0
                }
            }
            Gaus => gaus(x, a(1, 0.0), a(2, 1.0), a(3, 0.0) != 0.0),
            Landau => landau(x, a(1, 0.0), a(2, 1.0), a(3, 0.0) != 0.0),
            BreitWigner => oxistats::breit_wigner(x, a(1, 0.0), a(2, 1.0)),
            BreitWignerPdf => oxistats::breit_wigner(x, a(2, 0.0), a(1, f64::NAN)),
            NormalPdf => {
                let sigma = a(1, 1.0);
                let z = (x - a(2, 0.0)) / sigma;
                (-0.5 * z * z).exp() / ((2.0 * std::f64::consts::PI).sqrt() * sigma.abs())
            }
            CrystalBall => {
                crystalball(x, a(1, f64::NAN), a(2, f64::NAN), a(3, f64::NAN), a(4, 0.0))
            }
            CrystalBallPdf => {
                crystalball_pdf(x, a(1, f64::NAN), a(2, f64::NAN), a(3, f64::NAN), a(4, 0.0))
            }
            LandauPdf => {
                let xi = a(1, 1.0);
                if xi > 0.0 {
                    oxistats::landau(x, a(2, 0.0), xi)
                } else {
                    0.0
                }
            }
            BetaPdf => beta_pdf(x, a(1, f64::NAN), a(2, f64::NAN)),
            CauchyPdf => {
                let (b, d) = (a(1, 1.0), x - a(2, 0.0));
                b / (std::f64::consts::PI * (d * d + b * b))
            }
            ChiSquaredPdf => chisquared_pdf(x - a(2, 0.0), a(1, f64::NAN)),
            ExponentialPdf => {
                let (lambda, d) = (a(1, f64::NAN), x - a(2, 0.0));
                if d < 0.0 {
                    0.0
                } else {
                    lambda * (-lambda * d).exp()
                }
            }
            FDistributionPdf => fdistribution_pdf(x - a(3, 0.0), a(1, f64::NAN), a(2, f64::NAN)),
            GammaPdf => gamma_pdf(x - a(3, 0.0), a(1, f64::NAN), a(2, f64::NAN)),
            LognormalPdf => {
                let (m, sd, d) = (a(1, f64::NAN), a(2, f64::NAN), x - a(3, 0.0));
                if d <= 0.0 {
                    0.0
                } else {
                    let t = (d.ln() - m) / sd;
                    (-(t * t) / 2.0).exp() / (d * sd.abs() * (2.0 * std::f64::consts::PI).sqrt())
                }
            }
            PoissonPdf => poisson_pdf(x, a(1, f64::NAN)),
            TDistributionPdf => {
                let (r, d) = (a(1, f64::NAN), x - a(2, 0.0));
                (oxistats::gammaln((r + 1.0) / 2.0) - oxistats::gammaln(r / 2.0)).exp()
                    / (std::f64::consts::PI * r).sqrt()
                    * (1.0 + d * d / r).powf(-(r + 1.0) / 2.0)
            }
            UniformPdf => {
                let (lo, hi, d) = (a(1, f64::NAN), a(2, f64::NAN), x - a(3, 0.0));
                if d < hi && d >= lo {
                    1.0 / (hi - lo)
                } else {
                    0.0
                }
            }
            Poisson => {
                let mu = a(1, f64::NAN);
                if x < 0.0 {
                    0.0
                } else if x == 0.0 {
                    1.0 / mu.exp()
                } else {
                    (x * mu.ln() - mu - oxistats::gammaln(x + 1.0)).exp()
                }
            }
            BiGaussianPdf => bigaussian(
                x,
                a(1, f64::NAN),
                a(2, 1.0),
                a(3, 1.0),
                a(4, 0.0),
                a(5, 0.0),
                a(6, 0.0),
            ),
            Chebyshev(_) => chebyshev(x, args.get(1..).unwrap_or(&[])),
        }
    }
}

/// `TMath::Gaus`, with its edge cases: a zero width gives `1e30`, and a point
/// more than 39 widths away gives 0.
fn gaus(x: f64, mean: f64, sigma: f64, norm: bool) -> f64 {
    if sigma == 0.0 {
        return 1e30;
    }
    let arg = (x - mean) / sigma;
    if !(-39.0..=39.0).contains(&arg) {
        return 0.0;
    }
    let res = (-0.5 * arg * arg).exp();
    if norm {
        res / (2.506_628_274_631_000_2 * sigma)
    } else {
        res
    }
}

/// `TMath::Landau`: the Landau density of `(x − mpv)/sigma`, divided by
/// `sigma` when `norm`; 0 for a width that is not positive.
fn landau(x: f64, mpv: f64, sigma: f64, norm: bool) -> f64 {
    if sigma.is_nan() || sigma <= 0.0 {
        return 0.0;
    }
    // oxistats' Landau is the unit-area density, already divided by sigma.
    let density = oxistats::landau(x, mpv, sigma);
    if norm { density } else { density * sigma }
}

/// `ROOT::Math::crystalball_function`: unit peak, a power-law tail below
/// `mean − |alpha|·sigma` (above it for `alpha < 0`); 0 for `sigma < 0`.
fn crystalball(x: f64, alpha: f64, n: f64, sigma: f64, mean: f64) -> f64 {
    if sigma < 0.0 {
        return 0.0;
    }
    oxistats::crystal_ball(x, mean, sigma, alpha, n)
}

/// `ROOT::Math::crystalball_pdf`: the crystal ball over its area.
fn crystalball_pdf(x: f64, alpha: f64, n: f64, sigma: f64, mean: f64) -> f64 {
    if sigma < 0.0 {
        return 0.0;
    }
    if n <= 1.0 {
        return f64::NAN;
    }
    let abs_alpha = alpha.abs();
    let c = n / abs_alpha / (n - 1.0) * (-alpha * alpha / 2.0).exp();
    let d = (std::f64::consts::PI / 2.0).sqrt()
        * (1.0 + oxistats::erf(abs_alpha / std::f64::consts::SQRT_2));
    crystalball(x, alpha, n, sigma, mean) / (sigma * (c + d))
}

/// `ROOT::Math::beta_pdf`, with its edge cases at 0 and 1.
fn beta_pdf(x: f64, a: f64, b: f64) -> f64 {
    if !(0.0..=1.0).contains(&x) {
        return 0.0;
    }
    let edge = |p: f64, other: f64| {
        if p < 1.0 {
            f64::INFINITY
        } else if p > 1.0 {
            0.0
        } else {
            other
        }
    };
    if x == 0.0 {
        return edge(a, b);
    }
    if x == 1.0 {
        return edge(b, a);
    }
    (oxistats::gammaln(a + b) - oxistats::gammaln(a) - oxistats::gammaln(b)
        + x.ln() * (a - 1.0)
        + (-x).ln_1p() * (b - 1.0))
        .exp()
}

/// `ROOT::Math::chisquared_pdf` at `d = x − x0`.
fn chisquared_pdf(d: f64, r: f64) -> f64 {
    if d < 0.0 {
        return 0.0;
    }
    if d == 0.0 && r / 2.0 - 1.0 == 0.0 {
        return 0.5;
    }
    ((r / 2.0 - 1.0) * (d / 2.0).ln() - d / 2.0 - oxistats::gammaln(r / 2.0)).exp() / 2.0
}

/// `ROOT::Math::fdistribution_pdf` at `d = x − x0`.
fn fdistribution_pdf(d: f64, n: f64, m: f64) -> f64 {
    if n < 0.0 || m < 0.0 {
        return f64::NAN;
    }
    if d < 0.0 {
        return 0.0;
    }
    ((n / 2.0) * n.ln() + (m / 2.0) * m.ln() + oxistats::gammaln((n + m) / 2.0)
        - oxistats::gammaln(n / 2.0)
        - oxistats::gammaln(m / 2.0)
        + (n / 2.0 - 1.0) * d.ln()
        - ((n + m) / 2.0) * (m + n * d).ln())
    .exp()
}

/// `ROOT::Math::gamma_pdf` at `d = x − x0`.
fn gamma_pdf(d: f64, alpha: f64, theta: f64) -> f64 {
    if d < 0.0 {
        0.0
    } else if d == 0.0 {
        if alpha == 1.0 { 1.0 / theta } else { 0.0 }
    } else if alpha == 1.0 {
        (-d / theta).exp() / theta
    } else {
        ((alpha - 1.0) * (d / theta).ln() - d / theta - oxistats::gammaln(alpha)).exp() / theta
    }
}

/// `ROOT::Math::poisson_pdf(n, mu)`, `n` truncated as C++ converts a double
/// to `unsigned int`.
fn poisson_pdf(x: f64, mu: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    let n = x.trunc();
    if n > 0.0 && mu >= 0.0 {
        (n * mu.ln() - oxistats::gammaln(n + 1.0) - mu).exp()
    } else if mu >= 0.0 {
        (-mu).exp()
    } else {
        f64::NAN
    }
}

/// `ROOT::Math::bigaussian_pdf`.
fn bigaussian(x: f64, y: f64, sx: f64, sy: f64, rho: f64, x0: f64, y0: f64) -> f64 {
    let u = (x - x0) / sx;
    let v = (y - y0) / sy;
    let c = 1.0 - rho * rho;
    let z = u * u - 2.0 * rho * u * v + v * v;
    (-z / (2.0 * c)).exp() / (2.0 * std::f64::consts::PI * sx * sy * c.sqrt())
}

/// `Σ cₖ·Tₖ(x)` by Clenshaw's recurrence, as `ROOT::Math::ChebyshevN`.
fn chebyshev(x: f64, c: &[f64]) -> f64 {
    let (mut b1, mut b2) = (0.0, 0.0);
    for &ck in c.iter().skip(1).rev() {
        let b0 = 2.0 * x * b1 - b2 + ck;
        b2 = b1;
        b1 = b0;
    }
    c.first().copied().unwrap_or(0.0) + x * b1 - b2
}

/// Look `name` up for a call with `argc` arguments: the function, and whether
/// ROOT accepts this spelling. `None` if no function has this name and arity.
pub(crate) fn lookup(name: &str, argc: usize) -> Option<(Function, bool)> {
    let fits = |f: Function| {
        let (lo, hi) = f.arity();
        (lo..=hi).contains(&argc)
    };
    NAMES
        .iter()
        .find(|(n, f, _)| *n == name && fits(*f))
        .map(|(_, f, root)| (*f, *root))
        .or_else(|| {
            // Any case, any of ROOT's namespaces: `SQRT`, `Math::sqrt`.
            let bare = strip_namespace(name).to_ascii_lowercase();
            NAMES
                .iter()
                .find(|(n, f, _)| strip_namespace(n).to_ascii_lowercase() == bare && fits(*f))
                .map(|(_, f, _)| (*f, false))
        })
}

/// Whether some function, of whatever arity, has the name `name` (ignoring
/// case and namespaces).
pub(crate) fn known_name(name: &str) -> Option<Function> {
    let bare = strip_namespace(name).to_ascii_lowercase();
    NAMES
        .iter()
        .find(|(n, _, _)| strip_namespace(n).to_ascii_lowercase() == bare)
        .map(|(_, f, _)| *f)
}

/// The ROOT names closest to `name`, for a "did you mean" hint.
pub(crate) fn closest(name: &str) -> Option<&'static str> {
    let bare = strip_namespace(name).to_ascii_lowercase();
    NAMES
        .iter()
        .filter(|(_, _, root)| *root)
        .map(|(n, _, _)| {
            (
                edit_distance(&bare, &strip_namespace(n).to_ascii_lowercase()),
                *n,
            )
        })
        .filter(|(d, _)| *d <= 2)
        .min_by_key(|(d, _)| *d)
        .map(|(_, n)| n)
}

fn strip_namespace(name: &str) -> &str {
    name.rsplit("::").next().unwrap_or(name)
}

/// Levenshtein distance, for short names.
pub(crate) fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut row: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut prev = row[0];
        row[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let cur = row[j + 1];
            row[j + 1] = (prev + usize::from(ca != *cb)).min(row[j] + 1).min(cur + 1);
            prev = cur;
        }
    }
    row[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_function_has_a_root_name_and_evaluates() {
        let mut seen = Vec::new();
        for (_, f, _) in NAMES {
            if seen.contains(f) {
                continue;
            }
            seen.push(*f);
            assert_ne!(f.root_name(), "?", "{f:?}");
            let (lo, hi) = f.arity();
            assert!(lo <= hi && hi <= MAX_ARGS, "{f:?}");
            let args = vec![0.4; hi];
            let _ = f.eval(&args[..lo]);
            let _ = f.eval(&args);
        }
    }

    #[test]
    fn chebyshev_is_roots() {
        // ROOT: Chebyshev3(-2.5, 1, 2, 3, 4) = -189.5
        assert!((Chebyshev(3).eval(&[-2.5, 1.0, 2.0, 3.0, 4.0]) + 189.5).abs() < 1e-12);
        assert_eq!(Chebyshev(0).eval(&[7.0, 3.0]), 3.0);
    }
}
