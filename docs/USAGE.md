# Usage & API Guide (`tupa-epistemic-verifier v0.4.0`)

`tupa-epistemic-verifier` builds cleanly with zero external crates (`[dependencies]` is empty) on Linux, macOS, and `wasm32-unknown-unknown`.

## 1. CLI Quickstart

```bash
# Run full deterministic release test suite (including planted-defect mutants)
cargo test --release

# Self-check core epistemic scoring, JTMS retraction, citation & ORCID guards
cargo run --release

# Offline Lean 4 / Riemann benchmark trust-boundary audit
cargo run --release -- riemann-audit

# Audit an untrusted Lean 4 file against soundness holes & toolchain advisories
cargo run --release -- audit-lean path/to/Solution.lean --toolchain leanprover/lean4:v4.33.1

# Offline biomolecular & 3D interaction severe-testing battery (10/10 mutants killed, S_Mayo = 1.0000)
cargo run --release -- structural-selftest

# Audit an AlphaFold / Boltz / RFdiffusion mmCIF file against physical geometry and an experimental PDB reference
cargo run --release -- audit-protein AF-P47992-F1-model_v6.cif --ref pdb_00002jp1.cif --ref-chain A

# Audit a protein-ligand mmCIF complex for backbone steric clashes and 3D shape/ESP sanity
cargo run --release -- audit-ligand pdb_00004tzk.cif --ligand 641 --chain A

# Universal 5-order Vovk/CTW compression + FIPS 180-4 SHA-256 provenance digest
cargo run --release -- compress path/to/dataset.csv
```

### Exit Codes (Fail-Closed CI Semantics)
- `0`: All deterministic physical, cryptographic, or trust-boundary checks passed.
- `1`: I/O or CLI argument error.
- `2`: Verification gate held or failed (e.g., Lean 4 soundness hole / pending external kernel replay, severe protein/ligand steric clash, covalent chain break, or static-model overconfidence `|pLDDT/100 - lDDT-Ca| > 0.25`).

## 2. Rust Library API Example

```rust
use tupa_epistemic_verifier::{
    audit_protein_backbone, compare_model_vs_pdb_lddt_ca,
    parse_mmcif_backbone_chain, run_ligand_severe_testing_battery,
    run_protein_severe_testing_battery, sha256_hex,
};

fn verify_model_pair(af_cif: &str, pdb_cif: &str) {
    let model = parse_mmcif_backbone_chain(af_cif, Some("A"));
    let exp = parse_mmcif_backbone_chain(pdb_cif, Some("A"));
    let geom = audit_protein_backbone(&model);
    let cmp = compare_model_vs_pdb_lddt_ca(&model, &exp, 15.0).expect("aligned overlap");
    println!(
        "sha256={} lddt_ca={:.4} abs_cal_err={:.4} clashes={}",
        sha256_hex(af_cif.as_bytes()),
        cmp.lddt_ca,
        cmp.abs_calibration_error,
        geom.steric_clashes_ca
    );
    assert!(run_protein_severe_testing_battery(false, false).passed);
    assert!(run_ligand_severe_testing_battery(false, false).passed);
}
```

## 3. Epistemic Boundaries
Passing `audit-protein`, `audit-ligand`, or `audit-lean` checks deterministic geometry, calibration, or static syntax/replay conditions only. It **never** validates clinical pathogenicity, dynamic conformational completeness, binding affinity ($K_d$/$\mathrm{IC}_{50}$), or physical truth without wet-lab or independent kernel validation.
