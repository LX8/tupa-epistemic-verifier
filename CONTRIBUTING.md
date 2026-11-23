# Contributing to `tupa-epistemic-verifier`

Thank you for helping strengthen open, offline-verifiable epistemic and formal-proof tooling.

## Core Engineering Invariants

1. **Zero external crates (`std`-only Rust):** Keep `Cargo.toml` dependency-free so any auditor can inspect, compile, and replay every check offline without supply-chain risk.
2. **Release verification:** Always test and benchmark with `--release`:
   ```bash
   cargo test --release
   cargo run --release -- riemann-audit
   ```
3. **Planted-defect proof for every new check:** Every new Lean 4 soundness rule, citation hygiene check, or proper-scoring invariant must include a unit test that plants a deliberate defect and proves the check fails closed.
4. **Fail-closed demarcation:** Never label an archive verified without independent external kernel replay (`leanprover/comparator` + `nanoda`). Unverified or missing evidence is always `UNVERIFIED` or `NOT_MEASURED`, never a pass.

## High-Value Contributions ("Good First Issues")

- Additional Lean 4 metaprogramming / kernel-bypass test fixtures (`tests/` or `#[cfg(test)]`).
- Edge-case tests for comment/string stripping (nested `/- ... -/` block comments, raw string literals).
- Additional ISO 7064 / DOI / arXiv synthetic-hallucination regression tests.
- WebAssembly (`wasm32-unknown-unknown`) smoke checks and documentation improvements.

## Submitting Changes

1. Fork the repository and create a focused branch (`feat/...` or `fix/...`).
2. Run `cargo test --release` and verify zero compiler warnings.
3. Open a pull request explaining **what defect or capability it covers** and **how you verified it fails closed**.

## Licence

All contributions are licensed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)** — see [LICENSE](./LICENSE).
