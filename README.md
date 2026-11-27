# Tupã Epistemic Verifier (`tupa-epistemic-verifier`)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24476043.svg)](https://doi.org/10.5281/zenodo.24476043) [![ORCID](https://img.shields.io/badge/ORCID-0009--0009--0139--3700-A6CE39)](https://orcid.org/0009-0009-0139-3700) [![CI](https://github.com/LX8/tupa-epistemic-verifier/actions/workflows/ci.yml/badge.svg)](https://github.com/LX8/tupa-epistemic-verifier/actions/workflows/ci.yml) [![REUSE status](https://img.shields.io/badge/REUSE-3.3%20compliant-2563eb)](./REUSE.toml)

**Local-First Open-Source Severe-Testing Verifier for Scientific Provenance, AI Formal Proofs, and Biomolecular Structures (`Rust` / `WebAssembly`, Zero Dependencies, `GPL-3.0-or-later`)**

**Lokale Open-Source-Prüfbank gegen KI-Halluzinationen, Zitationsfälschungen, fehlerhafte Proteinstrukturen und nicht-reproduzierbare Datenanalysen — ohne Cloud-Zwang, ohne GPU und ohne Datenabfluss.**

---

## Overview / Überblick

`tupa-epistemic-verifier` (`v0.4.0`) is a zero-dependency, 100% Rust standard-library (`wasm32-unknown-unknown` + native CLI) severe-testing engine designed for independent researchers, computational biologists, and open-science maintainers:

1. **AI-Generated Lean 4 Formal-Proof Soundness & Trust-Boundary Auditor (`audit-lean`, `riemann-audit`):** Strips comments/strings, detects soundness holes (`sorry`, `admit`, `axiom`, `native_decide`, `@[extern]`, `unsafe`, `skipKernelTC`) and elaborator bypasses (`MVarId.assign`, `kexactB`, `kcopy`), checks `lean-toolchain` advisories, and enforces multi-kernel (`leanprover/comparator` + `nanoda`) demarcation.
2. **Biomolecular Structure & AlphaFold Calibration Auditor (`audit-protein`, `structural-selftest`):** Zero-dependency mmCIF `_atom_site` backbone parser, exact IUPAC $(\phi,\psi)$ Ramachandran geometry, $\mathrm{C}_\alpha$ bond-break and steric-clash detector, and superposition-free Mariani (2013) $\mathrm{lDDT\text{-}C}_\alpha$ vs. `pLDDT` calibration auditing — catching static-model overconfidence on metamorphic fold-switchers and flagging intrinsically disordered regions.
3. **3D Molecular Physical-Field & Chirality Auditor (`audit-ligand`, `structural-selftest`):** Closed-form Grant–Pickup Gaussian 3D volumetric shape Tanimoto, Coulombic electrostatic potential (ESP) Tanimoto and signed Hodgkin index, directional 3D pharmacophore overlap (`ShEPhERD`-style physics), $SO(3)$-invariant signed tetrahedral chirality checks (PoseBusters failure mode), and `HETATM` protein–ligand steric-clash auditing.
4. **Proper Scoring, Universal Coding & Non-Monotonic Truth Maintenance (`compress`, default self-check):** Calibrated `LogSumExp` cross-entropy, bounded Savage–Bregman Tsallis-$q$ scores, 5-order Vovk/CTW compression, Doyle non-monotonic JTMS (`0` zombie truths), FIPS 180-4 SHA-256 provenance, ISO 7064 ORCID validation, and synthetic DOI/arXiv blocking.

## Quick Verification

```bash
cargo test --release
cargo run --release
cargo run --release -- riemann-audit
cargo run --release -- structural-selftest
cargo run --release -- audit-protein path/to/AF-model.cif --ref path/to/pdb-ref.cif
cargo run --release -- audit-ligand path/to/complex.cif --ligand 641 --chain A
```

See [docs/USAGE.md](./docs/USAGE.md) for full CLI flags, exit-code semantics, and Rust library API examples, and [paper/paper.md](./paper/paper.md) for the JOSS software manuscript draft.

## License, Legal Scope & Governance

- **FSFE REUSE 3.3 Compliant ([REUSE.toml](./REUSE.toml)):** Executable Rust/WebAssembly code (`src/**`, `Cargo.toml`, `.github/**`) is licensed under **GNU General Public License v3.0 or later (`GPL-3.0-or-later`, [LICENSE](./LICENSE))**; the JOSS manuscript (`paper/**`) under `CC-BY-4.0`; documentation under `CC-BY-SA-4.0`; and metadata manifests ([CITATION.cff](./CITATION.cff), [codemeta.json](./codemeta.json)) under `CC0-1.0`.
- **Security, EU CRA & Non-Clinical Boundary:** See [SECURITY.md](./SECURITY.md) (90-day vulnerability disclosure, EU Cyber Resilience Act Recital 18 non-commercial open-source scope, non-clinical / non-MDR-2017/745 disclaimer, GPLv3 §§ 15–16 warranty notice).
- **Five-Layer Ethics Scope:** See [ETHICS.md](./ETHICS.md).
- **Contributing & DCO 1.1:** See [CONTRIBUTING.md](./CONTRIBUTING.md) for zero-dependency std Rust rules, `git commit -s` DCO sign-off, and AI-assisted contribution disclosure.

## Project & Community

- **Open Empirical Study Bundles & Preprints:** [lx8-tupa-science.web.app/papers/](https://lx8-tupa-science.web.app/papers/)
- **LX8 Labs:** [lx8labs.com/research](https://lx8labs.com/research/) · [LinkedIn](https://www.linkedin.com/company/lx8labs) · [YouTube](https://www.youtube.com/@lx8labs)
