//! The formula API: parameters as ROOT numbers them, errors, strict parsing,
//! rendering, introspection and limits.

use oxifit_formula::{Formula, Function};

fn parse(s: &str) -> Formula {
    Formula::parse(s).unwrap_or_else(|e| panic!("{s}: {e}"))
}

#[test]
fn a_root_formula_read_back_keeps_its_parameters() {
    // ROOT stores `[1]*x+[0]` as `[p1]*x+[p0]`; `[p1]` is parameter 1, not
    // the first name seen.
    let f = parse("[p1]*x+[p0]");
    assert_eq!(f.eval1(1.3, &[1.0, 1.5]), 1.5 * 1.3 + 1.0);
    assert_eq!(f, parse("[1]*x+[0]"));
    // A formula written back and read again is the same formula.
    let back = parse(&f.to_root());
    assert_eq!(back, f);
    assert_eq!(back.eval1(1.3, &[1.0, 1.5]), f.eval1(1.3, &[1.0, 1.5]));
    // And two formulas with parameters the other way round are not equal.
    assert_ne!(parse("[0]*x+[1]"), parse("[1]*x+[0]"));
}

#[test]
fn named_parameters_are_numbered_as_root_numbers_them() {
    let f = parse("[b]*x+[a]");
    assert_eq!(f.param_names(), ["a", "b"]);
    assert_eq!(f.eval1(1.3, &[1.0, 1.5]), 1.5 * 1.3 + 1.0); // ROOT: 2.95
    assert_eq!(f.param_index("b"), Some(1));
    let f = parse("[mean]*x+[p1]");
    assert_eq!(f.param_names(), ["p1", "mean"]);
    assert_eq!(f.param_index("1"), Some(0)); // `[1]` and `[p1]` are one parameter
    assert_eq!(f.param_index("p1"), Some(0));
    let f = parse("[1]*x+[mean]");
    assert_eq!(f.param_names(), ["p0", "p1", "mean"]);
    assert_eq!(f.param_index("1"), Some(1));
    assert_eq!(f.param_index("01"), Some(1));
    assert_eq!(parse("gaus").param_names(), ["Constant", "Mean", "Sigma"]);
    assert_eq!(
        parse("gaus+pol0(3)").param_names(),
        ["p0", "p1", "p2", "p3"]
    );
    assert_eq!(
        parse("crystalball").param_names(),
        ["Constant", "Mean", "Sigma", "Alpha", "N"]
    );
}

#[test]
fn parameter_names_can_be_set_from_a_file() {
    let f = parse("[0]*x+[1]")
        .with_param_names(&["slope", "offset"])
        .unwrap();
    assert_eq!(f.param_names(), ["slope", "offset"]);
    assert_eq!(f.to_string(), "[slope]*x+[offset]");
    assert_eq!(f.to_root(), "[p0]*x+[p1]");
    let f = parse("[0]*x");
    assert!(f.clone().with_param_names(&["a", "b"]).is_err());
    assert!(f.clone().with_param_names(&["a b"]).is_err());
    assert!(parse("[0]*x+[1]").with_param_names(&["a", "a"]).is_err());
}

#[test]
fn errors_point_at_the_problem_and_suggest_a_fix() {
    let e = Formula::parse("[0]*bogus(x)").unwrap_err();
    assert_eq!(e.span(), 4..9);
    assert_eq!(e.position(), 4);
    assert_eq!(
        e.to_string(),
        "unknown name `bogus`\n  [0]*bogus(x)\n      ^^^^^\n  help: the variables are x, y, z and t; \
         parameters are written `[0]` or `[name]`"
    );
    let hint = |s: &str| Formula::parse(s).unwrap_err().hint().map(str::to_owned);
    assert_eq!(
        hint("x & y").as_deref(),
        Some("ROOT writes logical AND as `&&`")
    );
    assert_eq!(
        hint("x | y").as_deref(),
        Some("ROOT writes logical OR as `||`")
    );
    assert_eq!(
        hint("~x").as_deref(),
        Some("ROOT writes logical NOT as `!`")
    );
    assert!(hint("x % 2").unwrap().contains("fmod"));
    assert_eq!(hint("sqr(x)").as_deref(), Some("did you mean `sqrt`?"));
    assert!(hint("2x").unwrap().contains("implicit multiplication"));
    let e = Formula::parse("sin(x, 2)").unwrap_err();
    assert_eq!(e.message(), "`sin` takes 1 argument, not 2");
    let e = Formula::parse("gausn").unwrap();
    assert_eq!(e.npar(), 3);
    let e = Formula::parse("1 +").unwrap_err();
    assert_eq!(e.span(), 3..3);
    let e = Formula::parse("xygausn").unwrap_err();
    assert!(e.message().contains("normalized form"), "{e}");
}

#[test]
fn strict_parsing_takes_only_what_root_compiles() {
    for (lenient, root) in [
        ("ln(x)", "log(x)"),
        ("sign(x)", "TMath::Sign(1.,x)"),
        ("min(x,1)", "TMath::Min(x,1.)"),
        ("power(x,2)", "pow(x,2.)"),
        ("SQRT(x)", "sqrt(x)"),
        ("pi()*x", "pi*x"),
        ("infinity>x", "TMath::Infinity()>x"),
    ] {
        let f = parse(lenient);
        assert_eq!(f.to_root(), root, "{lenient}");
        assert!(Formula::parse_root(lenient).is_err(), "{lenient}");
        assert!(Formula::parse_root(root).is_ok(), "{root}");
    }
    assert_eq!(
        Formula::parse_root("sign(x)").unwrap_err().hint(),
        Some("ROOT's `sign(a, b)` takes two arguments: `TMath::Sign(1., x)` is the sign of x")
    );
}

#[test]
fn display_is_for_people_and_to_root_for_root() {
    let f = parse("[N]*exp(-0.5*((x-[mu])/[sigma])^2) + [3]");
    assert_eq!(
        f.param_names(),
        ["p0", "p1", "p2", "p3", "N", "mu", "sigma"]
    );
    assert_eq!(f.to_string(), "[N]*exp(-0.5*((x-[mu])/[sigma])^2)+[3]");
    // `p0` … `p2` are unused, so `[4]`: ROOT's `[p4]` would not define them.
    assert_eq!(f.to_root(), "[4]*exp(-0.5*pow((x-[5])/[6],2.))+[3]");
    assert_eq!(parse("[b]*x+[a]").to_root(), "[p1]*x+[p0]");
    // Precedence survives a round trip.
    for s in [
        "-x^2",
        "(-x)^2",
        "2^3^2",
        "(2^3)^2",
        "1-(x-2)",
        "x-(-2)",
        "-(-x)",
        "--x",
        "x---2",
        "!!-2",
        "x*-y",
        "2^-x",
        "2/(x/3)",
        "x<1==y<1",
        "(x<1)<2",
        "x?1:y?2:3",
        "(x?1:2)?3:4",
        "!x>1",
        "!(x>1)",
        "x>1||y>1&&x<0",
        "(x>1||y>1)&&x<0",
    ] {
        let f = parse(s);
        let (v, p) = ([1.3, 0.7], [] as [f64; 0]);
        for text in [f.to_string(), f.to_root()] {
            let g = parse(&text);
            assert_eq!(g, f, "{s} → {text}");
            assert_eq!(
                g.eval(&v, &p).to_bits(),
                f.eval(&v, &p).to_bits(),
                "{s} → {text}"
            );
        }
    }
}

#[test]
fn latex_for_plot_labels() {
    assert_eq!(
        parse("gaus").to_latex(),
        r"\mathrm{Constant} \cdot e^{-0.5 \cdot \frac{x - \mathrm{Mean}}{\mathrm{Sigma}} \cdot \frac{x - \mathrm{Mean}}{\mathrm{Sigma}}}"
    );
    assert_eq!(
        parse("[0]+[1]*x^2").to_latex(),
        r"p_{0} + p_{1} \cdot {x}^{2}"
    );
    assert_eq!(
        parse("sqrt(x)/abs(y)").to_latex(),
        r"\frac{\sqrt{x}}{\left|y\right|}"
    );
    assert_eq!(parse("[a]*pi").to_latex(), r"a \cdot \pi");
}

#[test]
fn introspection() {
    let f = parse("[0]*TMath::Gaus(x,[1],[2]) + sqrt(y) + sin(x)*sin(y)");
    assert_eq!(f.variables(), ["x", "y"]);
    assert_eq!(
        f.functions(),
        [Function::Gaus, Function::Sqrt, Function::Sin]
    );
    assert_eq!(f.ndim(), 2);
    assert_eq!(parse("x[2]+t").variables(), ["z", "t"]);
    assert_eq!(parse("x[2]+t").ndim(), 4);
}

#[test]
fn too_few_values_give_nan_not_zero() {
    let f = parse("[0]*x + [1]*y");
    assert!(f.eval1(1.0, &[1.0, 1.0]).is_nan());
    assert!(f.eval(&[1.0, 2.0], &[1.0]).is_nan());
    assert_eq!(f.eval(&[1.0, 2.0], &[1.0, 1.0]), 3.0);
    let c = f.compile();
    assert!(c.eval(&[1.0], &[1.0, 1.0]).is_nan());
    assert_eq!(c.eval(&[1.0, 2.0], &[1.0, 1.0]), 3.0);
}

#[test]
fn compiled_matches_the_tree() {
    let xs: Vec<f64> = (0..50).map(|i| -2.0 + 0.1 * f64::from(i)).collect();
    for s in [
        "gaus(0)+pol2(3)",
        "[0]/2*(1+erf((x-[1])/([2]*sqrt(2))))",
        "x>0?x^3:-x^0.5",
        "landau+expo(3)",
        "crystalball",
        "cheb4",
        "TMath::Sign([0],x)*fmod(x,0.7)+hypot(x,[1])",
        "!(x>0)&&x<1||x>1.5",
    ] {
        let f = parse(s);
        let p: Vec<f64> = (0..f.npar()).map(|i| 0.3 + 0.4 * i as f64).collect();
        let c = f.compile();
        let mut out = vec![0.0; xs.len()];
        c.eval_many(&xs, &p, &mut out);
        for (x, o) in xs.iter().zip(&out) {
            let t = f.eval1(*x, &p);
            assert!(
                t.to_bits() == o.to_bits() || (t - o).abs() <= 1e-15 * t.abs().max(1.0),
                "{s} at {x}: {t} vs {o}"
            );
        }
    }
}

#[test]
fn variable_derivatives() {
    let f = parse("[0]*x^3 + sin(x*y) + hypot(x, y)");
    let dx = f.variable_derivative(0).unwrap();
    let dy = f.variable_derivative(1).unwrap();
    let (x, y, p) = (0.7, 1.9, [1.5]);
    let fx = |x: f64, y: f64| f.eval(&[x, y], &p);
    let h = 1e-6;
    let num_dx = (fx(x + h, y) - fx(x - h, y)) / (2.0 * h);
    let num_dy = (fx(x, y + h) - fx(x, y - h)) / (2.0 * h);
    assert!((dx.eval(&[x, y], &p) - num_dx).abs() < 1e-7);
    // Chebyshev polynomials have no x-derivative here.
    assert!(parse("cheb3").variable_derivative(0).is_none());
    assert!((dy.eval(&[x, y], &p) - num_dy).abs() < 1e-7);
    assert!(parse("TMath::Landau(x)").variable_derivative(0).is_none());
    assert!(f.variable_derivative(4).is_none());
    assert!(f.derivative(1).is_none());
}

#[test]
fn deep_nesting_is_an_error_not_a_stack_overflow() {
    let deep = format!("{}x{}", "(".repeat(100_000), ")".repeat(100_000));
    assert!(Formula::parse(&deep).is_err());
    let unary = format!("{}x", "-".repeat(100_000));
    assert!(Formula::parse(&unary).is_err());
    let calls = format!("{}x{}", "sin(".repeat(100_000), ")".repeat(100_000));
    assert!(Formula::parse(&calls).is_err());
    let ok = format!("{}x{}", "(".repeat(150), ")".repeat(150));
    assert!(Formula::parse(&ok).is_ok());
}

#[test]
fn huge_parameter_indices_are_errors_not_allocations_or_overflows() {
    for f in [
        "[1000000000000]*x",
        "[18446744073709551615]",
        "[p99999]",
        "gaus(1e300)",
        "pol3(18446744073709551615)",
        "pol100000",
        "expo(9999)",
    ] {
        assert!(Formula::parse(f).is_err(), "{f} should be rejected");
    }
    assert!(Formula::parse("pol20").is_ok());
    assert_eq!(Formula::parse("pol200").unwrap().npar(), 201);
    assert!(Formula::parse("[9999]*x").is_ok());
}

#[test]
fn constants_have_roots_values() {
    assert_eq!(parse("e").eval1(0.0, &[]), std::f64::consts::E);
    assert_eq!(parse("TMath::C()").eval1(0.0, &[]), 299_792_458.0);
    assert_eq!(parse("sqrt2").eval1(0.0, &[]), std::f64::consts::SQRT_2);
    assert!(Formula::parse("TMath::Pi").is_err()); // needs ()
}
