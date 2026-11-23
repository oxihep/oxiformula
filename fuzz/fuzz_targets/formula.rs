//! Parse arbitrary text as a formula; whatever parses must evaluate without
//! panicking (tree and compiled alike), and its ROOT spelling must parse back
//! to the same tree.
#![no_main]

use libfuzzer_sys::fuzz_target;
use oxiformula::Formula;

fuzz_target!(|data: &[u8]| {
    let Ok(source) = std::str::from_utf8(data) else {
        return;
    };
    // Deeply nested input is the parser's stack-depth case; keep inputs to a
    // size a real formula has.
    if source.len() > 4096 {
        return;
    }
    let Ok(f) = Formula::parse(source) else {
        return;
    };
    let vars = [0.5, -1.25, 3.0, 1e-3];
    let params: Vec<f64> = (0..f.npar()).map(|i| 0.1 + i as f64).collect();
    let vars = &vars[..f.ndim().min(vars.len())];
    let tree = f.eval(vars, &params);
    let compiled = f.compile().eval(vars, &params);
    assert!(
        tree.to_bits() == compiled.to_bits() || (tree - compiled).abs() <= 1e-12 * tree.abs().max(1.0) || (tree.is_nan() && compiled.is_nan()),
        "{source:?}: tree {tree}, compiled {compiled}"
    );
    let _ = f.to_latex();
    let _ = f.to_string();
    let root = f.to_root();
    let again = Formula::parse_root(&root)
        .unwrap_or_else(|e| panic!("{source:?} parsed, but its ROOT spelling {root:?} does not: {e}"));
    assert_eq!(again.expr(), f.expr(), "{source:?} → {root:?}");
    assert_eq!(again.npar(), f.npar(), "{source:?} → {root:?}");
});
