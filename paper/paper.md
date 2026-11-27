---
title: 'Tupã Epistemic Verifier: Zero-Dependency Rust and WebAssembly Severe Testing for Scientific Provenance, Formal Proofs, and Biomolecular Audits'
tags:
  - Rust
  - WebAssembly
  - reproducibility
  - provenance
  - structural biology
  - formal verification
  - severe testing
authors:
  - name: Alexei Lemos Cavalcante Ferreira
    orcid: 0009-0009-0139-3700
    corresponding: true
    affiliation: 1
affiliations:
  - name: LX8 Labs, Berlin, Germany
    index: 1
date: 27 November 2026
bibliography: paper.bib
---

# Summary

`tupa-epistemic-verifier` (`v0.4.0`) is a zero-dependency, standard-library Rust library and command-line interface (compiling natively and to `wasm32-unknown-unknown`) that performs deterministic, offline-reproducible verification of scientific artifacts. Rather than generating hypotheses or predictions via stochastic language or diffusion models, the package acts as a fail-closed epistemic and biophysical auditing layer. Its core modules implement calibrated `LogSumExp` cross-entropy and Savage–Bregman Tsallis-$q$ proper scoring, five-order Krichevsky–Trofimov / Vovk context-tree universal coding, Doyle non-monotonic justification-based truth maintenance (JTMS) [@Doyle:1979], FIPS 180-4 SHA-256 W3C RO-Crate provenance hashing [@Soiland-Reyes:2022], ISO 7064 ORCID check-digit and synthetic DOI/arXiv citation blocking, and trust-boundary auditing of AI-generated Lean 4 formal proofs. Its integrated structural module (`src/structural.rs`, CLI commands `audit-protein`, `audit-ligand`, and `structural-selftest`) provides pure-Rust mmCIF backbone/ligand parsing, superposition-free Mariani $\mathrm{lDDT\text{-}C}_\alpha$ structural calibration [@Mariani:2013], IUPAC $(\phi,\psi)$ Ramachandran geometry, and Grant–Pickup Gaussian 3D shape, Coulombic electrostatic potential (ESP), directional pharmacophore, and signed tetrahedral chirality auditing [@Grant:1995; @Buttenschoen:2024; @Adams:2025]. Every subsystem is validated against planted-defect negative-control fault mutants enforcing Mayo's severe-testing criterion ($S_{\mathrm{Mayo}} = 1 - P(\mathrm{pass} \mid \neg H) = 1.0000$) [@Mayo:2018].

# Statement of need

Computational science faces a twin verification crisis. First, data leakage, drifted metadata, fabricated references, and permissive continuous-integration fallbacks routinely allow flawed machine-learning claims to pass automated checks [@Kapoor:2023]. Second, generative deep-learning models in structural biology and chemistry—such as AlphaFold2 [@Jumper:2021], AlphaFold3 [@Abramson:2024], and 3D molecular diffusion frameworks [@Adams:2025]—produce static atomic coordinates at unprecedented scale while exhibiting documented blind spots: overconfident misprediction on metamorphic fold-switching proteins [@Chakravarty:2022], artificial secondary-structure order in intrinsically disordered regions (IDRs), unphysical ligand chirality inversions or pocket clashes [@Buttenschoen:2024], and near-zero correlation between wild-type confidence metrics and single-mutation stability ($\Delta\Delta G$) [@Pak:2023].

Independent researchers, computational biologists, and scientific software maintainers need a lightweight, self-contained verification engine that runs identically in local CI pipelines and browser WebAssembly sandboxes with zero cloud telemetry, zero GPU requirement, and zero transitive dependency drift. `tupa-epistemic-verifier` addresses this need by uniting cryptographic provenance, non-monotonic claim withdrawal, formal-proof soundness checks, and deterministic biophysical falsifiers inside a single auditable binary that refuses to promote unverified claims.

# State of the field

Existing scientific validation tools cover valuable but disjoint slices of the workflow and typically require heavy Python or C++ runtime stacks. In structural biology, `OpenStructure` computes exact $\mathrm{lDDT}$ scores [@Mariani:2013], `PoseBusters` evaluates physical and chemical validity of docked protein–ligand complexes [@Buttenschoen:2024], and `shepherd-score` provides Numba/Triton-accelerated 3D shape, electrostatic, and pharmacophore alignment [@Adams:2025]. In formal mathematics, `leanprover/comparator` performs independent kernel replay of Lean 4 proofs. In research data packaging, `ro-crate-py` and `reuse-tool` validate JSON-LD research containers [@Soiland-Reyes:2022] and SPDX licensing.

Contributing these heterogeneous checks into any single existing Python ecosystem was evaluated and rejected for three architectural reasons: (1) Python/C++ scientific stacks cannot compile cleanly to a zero-dependency client-side WebAssembly module for offline browser auditing; (2) domain-specific validators rarely share a unified non-monotonic truth-maintenance graph (JTMS) that automatically retracts downstream conclusions when an upstream dataset, citation, or structure fails audit; and (3) standard test suites verify positive cases without systematically executing planted-defect negative-control fault mutants to prove that pass-fabricating fallback paths cannot hide broken checks. `tupa-epistemic-verifier` is designed to complement—not replace—upstream tools like `PoseBusters`, `shepherd-score`, and `leanprover/comparator` by providing a portable, zero-dependency Rust/WASM verification core with uniform fail-closed demarcation semantics.

# Software design

The software is engineered around four strict invariants:

1. **Zero external crate dependencies (`[dependencies]` empty):** Both native binaries and `wasm32-unknown-unknown` builds rely exclusively on the Rust standard library (`edition = 2021`), eliminating supply-chain drift and guaranteeing deterministic floating-point and cryptographic behavior across Linux, macOS, and browsers.
2. **Explicit epistemic demarcation locks:** Every verification kernel returns structured receipts that separate software/geometric validity from physical or clinical truth. Static AlphaFold-vs-PDB geometry audits explicitly lock clinical pathogenicity adjudication out of scope, and mathematical boundary checks lock all six open Clay Millennium Problems as `OPEN PROBLEM` (`1/7` external control only) [@Ferreira:2026millennium].
3. **Non-monotonic dependency propagation (`Doyle JTMS`):** Claims record explicit justification edges to content-addressed SHA-256 dataset manifests and citation checks; invalidating any upstream node cascades deterministic withdrawal (`0` zombie truths).
4. **Severe-testing mutant batteries:** Every core module ships paired positive controls and planted negative-control mutants—such as injected Lean 4 `sorry`/`@[extern]` holes, synthetic DOIs, covalent $\mathrm{C}_\alpha\text{--}\mathrm{C}_\alpha$ chain breaks, non-bonded steric clashes ($< 3.0\text{ \AA}$), mirror-image tetrahedral chirality flips ($V_{\mathrm{model}} V_{\mathrm{ref}} < 0$), charge-inverted electrostatic decoys (Hodgkin index $-1.0$), and reversed pharmacophore vectors—requiring 100% mutant detection before release.

# Research impact statement

`tupa-epistemic-verifier` is publicly archived on Zenodo [@Ferreira:2026verifier] under `GPL-3.0-or-later` and underpins an open verification portal hosting 13 SHA-256 content-addressed empirical study bundles across five scientific domains (`https://lx8-tupa-science.web.app/papers/`). In published computational audits, the framework enabled: (i) a severe-testing demarcation benchmark across the seven Clay Millennium Prize Problems ($S_{\mathrm{Mayo}} = 1.0000$) [@Ferreira:2026millennium]; (ii) an automated self-correction of an earlier gravitomagnetic polarization estimate via independent Kerr geodesic integration; and (iii) live pure-Rust verification of AlphaFold DB v6 predictions against experimental RCSB PDB structures and 12 NCBI ClinVar GRCh38 human variants—reproducing high crystal accuracy on *Ideonella sakaiensis* PETase (`6EQE`, $\mathrm{lDDT\text{-}C}_\alpha = 0.9981$) while flagging $+32.1\%$ to $+41.9\%$ static-model overconfidence on metamorphic fold-switchers (`KaiB 2QKE`, `XCL1 2JP1`), `80.41%` whole-protein disorder in `BRCA1` (`P38398`), and low local confidence (`pLDDT = 51.22`) at actionable `EGFR p.Leu858Arg`.

# AI usage disclosure

Generative AI coding and editing assistants (including Claude and Gemini models) were used during software development for boilerplate test scaffolding, cross-checking documentation formatting, and copy-editing portions of this manuscript. All architectural decisions, mathematical formulations, biophysical invariants, negative-control fault mutants, citation DOIs, and live benchmark results were designed, reviewed, executed, and verified by the human author via deterministic `cargo test --release` suites and primary-source records. No generative AI model is listed as an author, and `tupa-epistemic-verifier` itself contains zero generative AI runtime components.

# References
