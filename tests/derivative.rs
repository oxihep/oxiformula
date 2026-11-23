//! Symbolic derivatives against central differences, for every function the
//! parser knows (and the operators, conditions and steps).

use oxifit_formula::Formula;

#[test]
fn derivatives_match_finite_differences() {
    let formulas = [
        "[0] + [1]*x + [2]*x^2",
        "[0]*exp(-0.5*((x-[1])/[2])^2)",
        "[0]/([1]+x) - [2]*x/[0]",
        "pow(x, [0]) + [1]^x + x^[0]^0.5",
        "sin([0]*x) + cos([1]*x) + tan([2]*0.3*x)",
        "asin([0]*0.1*x) + acos([1]*0.1*x) + atan([2]*x)",
        "sinh([0]*x) + cosh([1]*x) + tanh([2]*x)",
        "asinh([0]*x) + acosh(1.5+[1]*x*x) + atanh([2]*0.1*x)",
        "exp([0]*x) + expm1([1]*x) + log([2]+x) + log10([0]+x) + log2([1]+x)",
        "sqrt([0]+x) + cbrt([1]*x+2) + abs([2]*x-1)",
        "erf([0]*x) + erfc([1]*x)",
        "atan2([0]*x, [1]+x) + hypot([0]*x, [2]) + min([0]*x, [1]) + max([2]*x, [1])",
        "x > [1] ? [0]*(x-[1]) : [2]*x",
        "floor(x)*[0] + ceil(x)*[1] + sign(x-1)*[2] + pi*[0]",
        "sq([0]*x) + fmod([0]*x+3, [1]+2) + TMath::Sign([0]*x, x-0.2)",
        "log1p([0]*x) + expm1([1]*x) + [2]*TMath::Sign(1., x-1)",
        "ROOT::Math::Chebyshev2(x, [0], [1], [2])",
        "!(x > [0])*[1] + [2]",
        "TMath::Max([0]*x, [1]) + TMath::Min(x, [2])",
        "gaus",
        "expo + pol0(2)",
    ];
    let params: [f64; 3] = [1.3, 0.7, 2.1];
    let xs = [0.35, 1.1, 2.4];
    for src in formulas {
        let f = Formula::parse(src).unwrap_or_else(|e| panic!("{src}: {e:?}"));
        for k in 0..f.npar() {
            let d = f
                .derivative(k)
                .unwrap_or_else(|| panic!("{src}: no derivative"));
            for &x in &xs {
                let h = 1e-6 * f64::max(params[k].abs(), 1.0);
                let (mut up, mut down) = (params, params);
                up[k] += h;
                down[k] -= h;
                let numeric = (f.eval1(x, &up) - f.eval1(x, &down)) / (2.0 * h);
                let exact = d.eval1(x, &params);
                assert!(
                    (exact - numeric).abs() <= 1e-6 * numeric.abs().max(1.0),
                    "{src}, d/d[{k}] at x = {x}: {exact} vs {numeric}"
                );
            }
        }
    }
    // A function without a closed form here.
    assert!(
        Formula::parse("lgamma([0]*x)")
            .unwrap()
            .derivative(0)
            .is_none()
    );
    let f = Formula::parse("[0]*x").unwrap();
    assert_eq!(f.derivative(0).unwrap().eval1(2.0, &[5.0]), 2.0);
    assert!(f.derivative(1).is_none());
    // The densities ROOT's shortcuts call have none either; `gaus` and
    // `expo` expand to closed forms.
    assert!(Formula::parse("landau").unwrap().derivative(1).is_none());
}

#[test]
fn variable_derivatives_match_finite_differences() {
    let formulas = [
        "[0] + [1]*x + [2]*x^2",
        "[0]*exp(-0.5*((x-[1])/[2])^2)",
        "x^x + 2^x + x^0.5",
        "sin(x)*cos(x)/tan(x+0.1) + atan2(x, 2) + hypot(x, 3)",
        "sinh(x) + cosh(x) + tanh(x) + asinh(x) + atanh(x/3)",
        "log(x) + log10(x) + log2(x) + log1p(x) + sqrt(x) + cbrt(x) + erf(x)",
        "x > 1 ? x^2 : -x + abs(x - 2) + sq(x)",
        "gaus + pol1(3)",
    ];
    let params = [1.3, 0.7, 2.1, 0.4, 0.2];
    for src in formulas {
        let f = Formula::parse(src).unwrap_or_else(|e| panic!("{src}: {e:?}"));
        let d = f
            .variable_derivative(0)
            .unwrap_or_else(|| panic!("{src}: no derivative"));
        for x in [0.35, 1.6, 2.4] {
            let h = 1e-6;
            let numeric = (f.eval1(x + h, &params) - f.eval1(x - h, &params)) / (2.0 * h);
            let exact = d.eval1(x, &params);
            assert!(
                (exact - numeric).abs() <= 1e-6 * numeric.abs().max(1.0),
                "{src}, d/dx at {x}: {exact} vs {numeric}"
            );
        }
    }
}
