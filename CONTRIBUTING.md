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

## Submitting Changes & Legal Provenance (DCO 1.1)

1. Fork the repository and create a focused branch (`feat/...` or `fix/...`).
2. Run release tests, Wasm compilation check, and REUSE 3.3 compliance check with zero warnings:
   ```bash
   cargo test --release
   cargo check --release --lib --target wasm32-unknown-unknown
   reuse lint
   ```
3. **Developer Certificate of Origin (DCO 1.1):** Sign off every commit (`git commit -s`) certifying you have the legal right to contribute the work under the repository's open-source licenses (`https://developercertificate.org/`).
4. **AI-Assisted Contribution Disclosure:** If you used an AI coding/editing assistant while preparing a pull request, disclose it plainly in the PR description and confirm you personally reviewed, tested, and verified every invariant. The compiled verifier binary must remain 100% deterministic Rust/WebAssembly with zero generative AI runtime dependencies.
5. Open a pull request explaining **what defect or capability it covers** and **how you verified it fails closed**.

## Licensing & Safety Scope

- Executable source code (`src/**`, `Cargo.toml`, `.github/**`): **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)** — see [LICENSE](./LICENSE).
- JOSS manuscript (`paper/**`): **CC-BY-4.0**.
- Documentation (`README.md`, `docs/**`, `SECURITY.md`, `ETHICS.md`, `CONTRIBUTING.md`): **CC-BY-SA-4.0**.
- Metadata & build manifests (`CITATION.cff`, `codemeta.json`, `Cargo.lock`, `REUSE.toml`): **CC0-1.0**.
- See [SECURITY.md](./SECURITY.md) and [ETHICS.md](./ETHICS.md) for coordinated vulnerability disclosure, EU CRA non-commercial open-source boundaries, and non-clinical liability disclaimers.
