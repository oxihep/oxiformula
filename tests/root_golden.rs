//! Formulas checked against ROOT 6.40's `TFormula`: the parameters (how many,
//! their names and order), the value, and that ROOT compiles what
//! [`Formula::to_root`] writes and gets the same value from it.
//!
//! ROOT's answers are `data/root_formulas.tsv`, made by
//! `data/root_formulas.C` (which says how). Each formula is evaluated at
//! x, y, z, t = 1.3, 0.7, 0.2, 0.9 with parameters 1, 1.5, 2, … by index.

use std::path::PathBuf;

use oxifit_formula::Formula;

const VARS: [f64; 4] = [1.3, 0.7, 0.2, 0.9];

fn params(n: usize) -> Vec<f64> {
    (0..n).map(|i| 1.0 + 0.5 * i as f64).collect()
}

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn same(a: f64, b: f64) -> bool {
    a == b || (a.is_nan() && b.is_nan()) || (a - b).abs() <= 1e-10 * a.abs().max(1.0)
}

/// Write `data/root_formulas.in.tsv`, the input of `data/root_formulas.C`.
#[test]
#[ignore = "writes the input of data/root_formulas.C"]
fn write_inputs() {
    let list = std::fs::read_to_string(data("root_formulas.txt")).unwrap();
    let mut out = String::new();
    for source in list.lines().filter(|l| !l.is_empty()) {
        let f = Formula::parse_root(source).unwrap_or_else(|e| panic!("{source}: {e}"));
        out.push_str(&format!("{source}\t{}\n", f.to_root()));
    }
    std::fs::write(data("root_formulas.in.tsv"), out).unwrap();
}

#[test]
fn formulas_are_roots() {
    let table = std::fs::read_to_string(data("root_formulas.tsv")).unwrap();
    let list = std::fs::read_to_string(data("root_formulas.txt")).unwrap();
    let mut checked = 0;
    for row in table.lines().filter(|l| !l.starts_with('#')) {
        let cols: Vec<&str> = row.split('\t').collect();
        let source = cols[0];
        assert!(
            cols.len() == 7,
            "ROOT does not accept {source:?}; take it out of root_formulas.txt"
        );
        let (ndim, npar, names) = (cols[1], cols[2], cols[3]);
        let value: f64 = cols[4].parse().unwrap();
        let f = Formula::parse_root(source).unwrap_or_else(|e| panic!("{source}: {e}"));
        let ctx = format!("{source:?}");
        assert_eq!(f.npar().to_string(), npar, "{ctx}: parameters");
        assert_eq!(f.param_names().join(","), names, "{ctx}: names");
        // ROOT counts no variables for a constant formula; this crate, one.
        let root_ndim: usize = ndim.parse().unwrap();
        assert_eq!(f.ndim(), root_ndim.max(1), "{ctx}: variables");
        let p = params(f.npar());
        let v = f.eval(&VARS[..f.ndim()], &p);
        assert!(same(v, value), "{ctx}: {v} where ROOT gives {value}");
        let compiled = f.compile().eval(&VARS[..f.ndim()], &p);
        assert!(
            same(compiled, value),
            "{ctx}: compiled {compiled}, ROOT {value}"
        );

        // What it writes is what ROOT compiled.
        assert_eq!(
            f.to_root(),
            cols[5],
            "{ctx}: to_root changed; regenerate data/root_formulas.tsv (see root_formulas.C)"
        );
        assert_ne!(
            cols[6], "INVALID",
            "{ctx}: ROOT does not compile {:?}",
            cols[5]
        );
        let root_value: f64 = cols[6].parse().unwrap();
        assert!(
            same(root_value, value),
            "{ctx}: ROOT gives {root_value} for {:?}",
            cols[5]
        );
        let again = Formula::parse_root(&f.to_root()).unwrap();
        assert_eq!(again.npar(), f.npar(), "{ctx}: to_root's parameters");
        assert!(
            same(again.eval(&VARS[..f.ndim()], &p), value),
            "{ctx}: to_root's value"
        );
        checked += 1;
    }
    assert_eq!(checked, list.lines().filter(|l| !l.is_empty()).count());
}
