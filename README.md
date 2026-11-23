# oxiformula

[![CI](https://github.com/oxihep/oxiformula/actions/workflows/ci.yml/badge.svg)](https://github.com/oxihep/oxiformula/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Rust 1.100+](https://img.shields.io/badge/rust-1.100%2B-orange.svg)](https://www.rust-lang.org)
[![API docs](https://img.shields.io/badge/docs-API-blue.svg)](https://oxihep.github.io/oxiformula/oxiformula/)

ROOT's [`TFormula`](https://root.cern/doc/master/classTFormula.html) in pure
Rust: parse a formula, evaluate or compile it, differentiate it symbolically,
and write it back as ROOT compiles it. It follows ROOT closely — the same
parameters in the same order, the same shortcuts (`gaus`, `expo`, `polN`,
`crystalball`, `landau`, …), the same functions and constants — and its test
suite checks it against ROOT on 151 formulas, with ROOT itself rerun in CI.

```toml
[dependencies]
oxiformula = { git = "https://github.com/oxihep/oxiformula", tag = "v0.11.0" }
```

```rust
use oxiformula::Formula;

let f = Formula::parse("gaus")?;
assert_eq!(f.param_names(), ["Constant", "Mean", "Sigma"]);
let g = f.compile();                       // about as fast as ROOT's JIT
assert_eq!(g.eval(&[0.5], &[2.0, 0.5, 1.0]), 2.0);
assert_eq!(Formula::parse("[0]*ln(x)")?.to_root(), "[p0]*log(x)");
# Ok::<(), oxiformula::ParseError>(())
```

| What | How |
| --- | --- |
| Parse | `Formula::parse` (also takes `ln`, a one-argument `sign`, `min`), or `Formula::parse_root` for exactly what ROOT accepts |
| Evaluate | `eval`, or `compile()` for repeated evaluation and `eval_many` over a column |
| Differentiate | `derivative(k)` in a parameter, `variable_derivative(v)` in `x`, `y`, … — symbolic |
| Write | `to_root()` (what ROOT compiles), `to_root_named()` (as a file stores it), `to_latex()`, `Display` |
| Inspect | `expr()` (the tree), `param_names`, `param_index`, `variables`, `functions` |
| Errors | a `ParseError` with a span and, for the common mistakes, a hint (`&&` for `&`, `log` for `ln`) |

oxiformula was `oxifit-formula`, inside [oxifit](https://github.com/oxihep/oxifit),
until 0.11.0; its history is in this repository. It is used by
[oxiroot](https://github.com/oxihep/oxiroot) (`Func1D`/`Func2D`/`Func3D`) and
oxifit (`Model::from_formula`).
