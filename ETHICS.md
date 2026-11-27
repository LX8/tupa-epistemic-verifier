<!--
SPDX-FileCopyrightText: 2026 Alexei Lemos Cavalcante Ferreira <ferreira.alexei@gmail.com>
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Ethics & Five-Layer Safety Scope (`tupa-epistemic-verifier`)

Aligned with the canonical LX8 / Tupã Family layered-ethics framework, safety and lawful stewardship are enforced from the deepest architectural layer outward:

1. **Layer 1 — Intrinsic Harmlessness (Architectural):** `tupa-epistemic-verifier` is a passive, offline, read-only mathematical auditor. It does **not** generate novel biological agents, fold or design de novo proteins, synthesize wet-lab protocols, exfiltrate user data, or perform offensive security actions. It can only score, parse, and falsify claims against deterministic invariants and public experimental baselines.
2. **Layer 2 — Legal Copyleft & Open Licensing:** All executable verification code is licensed under `GPL-3.0-or-later`, documentation under `CC-BY-SA-4.0`, the JOSS manuscript under `CC-BY-4.0`, and machine-readable metadata under `CC0-1.0`, verified via FSFE **REUSE 3.3** (`REUSE.toml` + `LICENSES/`). Derivative modifications remain auditable in the public commons.
3. **Layer 3 — Non-Aggression & Non-Clinical Covenant:** Output is strictly restricted to defensive verification, open-science provenance, and educational auditability. Clinical, diagnostic, and harmful dual-use applications are explicitly out of scope (see [SECURITY.md](./SECURITY.md)).
4. **Layer 4 — Deterministic Audit Trail & Severe Testing:** Every verification domain carries planted-defect negative controls (`structural-selftest`, `riemann-audit`, `cargo test --release`) enforcing Mayo severe-testing power `S_Mayo = 1.0000` so false passes fail closed in CI.
5. **Layer 5 — Public Oversight & Transparency:** All source code, release tags, Zenodo archives (`doi:10.5281/zenodo.24476043`), GitHub Actions CI runs, and AI-assisted development disclosures remain publicly inspectable.
