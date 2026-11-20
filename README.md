# Tupã Epistemic Verifier (`tupa-epistemic-verifier`)

**Local-First Open-Source Verification Engine Against AI Hallucinations, Synthetic Citations, and Uncalibrated Uncertainty (`Rust` / `WebAssembly`, `GPL-3.0-or-later`)**

**Lokale Open-Source-Prüfbank gegen KI-Halluzinationen, Zitationsfälschungen und nicht-reproduzierbare Datenanalysen — ohne Cloud-Zwang, ohne GPU und ohne Datenabfluss.**

---

## Overview / Überblick

`tupa-epistemic-verifier` is a zero-dependency, 100% Rust standard-library (`wasm32-unknown-unknown` + native C-ABI) verification kernel designed for independent researchers, investigative data journalists, and civil-society organizations:

1. **AI-Generated Lean 4 Formal-Proof Soundness & Trust-Boundary Auditor (`audit-lean`, `riemann-audit`):** Strips comments/strings and scans untrusted Lean 4 proof modules for soundness holes (`sorry`, `sorryAx`, `admit`, custom `axiom`, `native_decide`, `@[extern]`, `@[implemented_by]`, `unsafe`, `skipKernelTC`), detects elaborator-bypassing metaprogramming (`MVarId.assign`, `kexactB`, `kcopy`), checks `lean-toolchain` against known runtime soundness advisories (`#14576` pre-`4.32.2`, `#15241`/`#15289` pre-`4.34.1`), and enforces multi-kernel (`leanprover/comparator` `v4.36.0` + `nanoda`) demarcation — demonstrated live on **Alpöge & Furman (arXiv:2608.13637v2, verified `67.2505%` critical-line bound)** vs the **2026-11-20 `175,187`-file RH/GRH Lean `v4.33.1` archive** (`UNVERIFIED_PENDING_EXTERNAL_KERNEL_REPLAY`).
2. **Proper Scoring & Information Geometry (`LogSumExp`, `Tsallis-q`, `sqrt-JSD`, `Fisher-Rao`) & Universal Vovk/CTW Compression (`compress`):** Evaluates calibrated cross-entropy directly on raw logits without overflow or `1e-12` clipping bias, bounded outlier-robust Savage–Bregman Tsallis $q$-scores, and 5-order Krichevsky–Trofimov / Vovk context-tree lossless byte coding.
3. **Non-Monotonic Truth Maintenance (`Doyle JTMS`):** Tracks dependencies from raw sources and datasets to downstream claims; when an upstream source is retracted or refuted, derived conclusions are deterministically withdrawn (`0` zombie claims).
4. **Cryptographic Provenance & Identifier Hygiene (`FIPS 180-4 SHA-256`, `ISO 7064 ORCID`, & Synthetic DOI/arXiv Blocker):** Verifies research container digests (W3C RO-Crate compatible), researcher identifiers, and blocks synthetic LLM citation hallucinations locally offline.

## Quick Verification

```bash
cargo test --release
cargo run --release -- riemann-audit
cargo run --release -- audit-lean path/to/Solution.lean --toolchain leanprover/lean4:v4.33.1
```

## License

Free and Open Source Software licensed under **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**. See [LICENSE](./LICENSE).
