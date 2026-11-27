<!--
SPDX-FileCopyrightText: 2026 Alexei Lemos Cavalcante Ferreira <ferreira.alexei@gmail.com>
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Security & Legal Scope Policy

## Reporting a Vulnerability or Audit Defect

If you discover a security vulnerability, numerical inaccuracy, silent pass-fabricating fallback (`S_Mayo < 1.0`), Lean 4 soundness-scanner bypass, or biomolecular coordinate/clash parsing defect in `tupa-epistemic-verifier`, please report it **privately** rather than opening a public GitHub issue:

- **Email:** [ferreira.alexei@gmail.com](mailto:ferreira.alexei@gmail.com)
- **Subject:** `[security][tupa-epistemic-verifier] <one-line summary>`
- **Please include:** affected module/subcommand, minimal offline reproducer (input, expected exit code/output, observed output), and severity assessment.
- **Coordinated disclosure timeline:** Default **90 days** from acknowledgement (shorter for low-impact defects; negotiable for complex parser/verification issues).

## Scope of Covered Defects

- **Fail-closed violations:** Any input or planted defect that silently returns exit code `0` or an overstated verification status instead of failing closed.
- **Parser & memory-safety defects:** Panics, integer overflows, or malformed-input crashes in the zero-dependency Rust/WebAssembly parsers (`src/lib.rs`, `src/structural.rs`, `src/main.rs`).
- **Numerical & cryptographic defects:** Deviations in FIPS 180-4 SHA-256 digests, IUPAC `(phi, psi)` dihedrals, Mariani (2013) `lDDT-Ca`, Grant–Pickup Gaussian 3D shape Tanimoto, Coulombic ESP/Hodgkin scores, or `SO(3)`-invariant tetrahedral chirality checks.
- **Licensing, citation, & supply-chain defects:** SPDX/REUSE 3.3 header gaps, drifted `CITATION.cff` / `codemeta.json` / `Cargo.toml` metadata, or broken DOI references.

## Legal Terms, Liability Boundary & Regulatory Non-Scope

1. **Free & Open-Source Software ("AS IS"):** Source code is licensed under **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)** and provided without warranty of any kind per GPLv3 §§ 15–16 (Disclaimer of Warranty and Limitation of Liability). Documentation is licensed under `CC-BY-SA-4.0`, the JOSS manuscript under `CC-BY-4.0`, and metadata/config files under `CC0-1.0`.
2. **Not a Clinical or Medical Device:** Structural biology and ClinVar residue-auditing outputs (`audit-protein`, `audit-ligand`) are strictly for **deterministic computational research, reproducibility auditing, and educational verification**. This software is **not** a medical device or in-vitro diagnostic under EU MDR 2017/745, EU IVDR 2017/746, or US FDA regulations and must **never** be used for clinical diagnosis, patient triage, or therapeutic decision-making.
3. **Non-Commercial Open-Source Development Boundary (EU CRA):** This repository is maintained as an independent, non-commercial public-interest open-source project outside commercial product placement under the EU Cyber Resilience Act (Regulation (EU) 2024/2847, Recital 18). Commercial integrators embedding this library assume full regulatory, conformity-assessment, and civil-liability obligations in their own jurisdiction (including German BGB / ProdHaftG and Brazilian Civil Code Arts. 186/927).
4. **Offline Privacy & Zero Telemetry:** The binary performs zero network calls, collects zero user telemetry, and ships zero generative AI or cloud LLM runtime components (aligned with GDPR / LGPD data-minimization principles).
