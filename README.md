# Tupã Epistemic Verifier (`tupa-epistemic-verifier`)

**Local-First Open-Source Verification Engine Against AI Hallucinations, Synthetic Citations, and Uncalibrated Uncertainty (`Rust` / `WebAssembly`, `GPL-3.0-or-later`)**

**Lokale Open-Source-Prüfbank gegen KI-Halluzinationen, Zitationsfälschungen und nicht-reproduzierbare Datenanalysen — ohne Cloud-Zwang, ohne GPU und ohne Datenabfluss.**

---

## Overview / Überblick

`tupa-epistemic-verifier` is a zero-dependency, 100% Rust standard-library (`wasm32-unknown-unknown` + native C-ABI) verification kernel designed for independent researchers, investigative data journalists, and civil-society organizations:

1. **Proper Scoring & Information Geometry (`LogSumExp`, `Tsallis-q`, `sqrt-JSD`, `Fisher-Rao`):** Evaluates calibrated cross-entropy directly on raw logits without overflow or `1e-12` clipping bias, alongside bounded outlier-robust Savage–Bregman Tsallis $q$-scores and Riemannian Fisher–Rao distances.
2. **Non-Monotonic Truth Maintenance (`Doyle JTMS`):** Tracks dependencies from raw sources and datasets to downstream claims; when an upstream source is retracted or refuted, derived conclusions are deterministically withdrawn (`0` zombie claims).
3. **Cryptographic Provenance & Identifier Hygiene (`FIPS 180-4 SHA-256` & `ISO 7064 ORCID`):** Verifies research container digests (W3C RO-Crate compatible) and researcher identifiers locally offline.

## Quick Verification

```bash
cargo test --release
```

## License

Free and Open Source Software licensed under **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**. See [LICENSE](./LICENSE).
