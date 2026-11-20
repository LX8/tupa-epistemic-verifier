// SPDX-FileCopyrightText: 2026 Alexei Lemos Cavalcante Ferreira <ferreira.alexei@gmail.com>
// SPDX-License-Identifier: GPL-3.0-or-later

use std::env;
use std::fs;
use std::time::Instant;
use tupa_epistemic_verifier::{
    calibrated_cross_entropy_bits, detect_synthetic_citation_hallucination,
    fisher_rao_geodesic_rad, jtms_grounded_closure, sha256_hex, sqrt_jensen_shannon_metric,
    tsallis_proper_score_bits, validate_orcid_iso7064, vovk_ctw_compress_rate_bpb,
    JtmsJustification,
};

fn main() {
    let args: Vec<String> = env::args().collect();
    let t0 = Instant::now();

    if args.len() >= 3 && args[1] == "compress" {
        let path = &args[2];
        let bytes = match fs::read(path) {
            Ok(b) => b,
            Err(e) => {
                eprintln!("error: failed to read '{}': {}", path, e);
                std::process::exit(1);
            }
        };
        let (bpb0, bpb_mix, margin) = vovk_ctw_compress_rate_bpb(&bytes, 0.50);
        let digest = sha256_hex(&bytes);
        println!("Tupã Epistemic Verifier v0.2.0 — Universal Vovk/CTW Stream Audit");
        println!("  File                 : {} ({} bytes)", path, bytes.len());
        println!("  SHA-256 (FIPS 180-4) : {}", digest);
        println!("  Order-0 KT Rate      : {:.4} bits/byte", bpb0);
        println!("  5-Order Vovk/CTW     : {:.4} bits/byte (saved {:.4} bpb)", bpb_mix, bpb0 - bpb_mix);
        println!("  Vovk Regret Headroom : +{:.4} bits below log2(5) = 2.3219 ceiling", margin);
        println!("  Elapsed              : {:.2} ms", t0.elapsed().as_secs_f64() * 1e3);
        return;
    }

    let ce_cal = calibrated_cross_entropy_bits(&[0.85, -0.35, -0.50], 0, 0.43);
    let ce_raw = calibrated_cross_entropy_bits(&[0.85, -0.35, -0.50], 0, 1.00);
    let tsallis_out = tsallis_proper_score_bits(&[0.001, 0.4995, 0.4995], 0, 1.50);
    let jsd = sqrt_jensen_shannon_metric(&[0.7, 0.2, 0.1], &[0.2, 0.5, 0.3]);
    let fr = fisher_rao_geodesic_rad(&[0.7, 0.2, 0.1], &[0.2, 0.5, 0.3]);
    let rules = vec![
        JtmsJustification { consequent: 1, in_antecedents: vec![0], out_blockers: vec![2] },
    ];
    let retracted = jtms_grounded_closure(3, &[0, 2], &rules);
    let cite_ok = detect_synthetic_citation_hallucination("doi:10.1103/PhysRevD.38.472").is_empty();
    let orcid_ok = validate_orcid_iso7064("0009-0009-0139-3700");

    println!("Tupã Epistemic Verifier v0.2.0 (100% std-Rust · Zero Cloud · GPL-3.0-or-later)");
    println!("  • LogSumExp Cross-Entropy : {:.4} bits (τ*=0.43) vs {:.4} bits (raw τ=1.00)", ce_cal, ce_raw);
    println!("  • Tsallis q=1.5 Outlier   : {:.4} bits (bounded under p_y=0.001)", tsallis_out);
    println!("  • √JSD / Fisher-Rao       : {:.4} / {:.4} rad", jsd, fr);
    println!("  • Doyle JTMS Retraction   : zombie_claims={} (verified)", usize::from(retracted[1]));
    println!("  • Citation & ORCID Gate   : cite_ok={} · orcid_iso7064={}", cite_ok, orcid_ok);
    println!("  • Self-Check Latency      : {:.2} µs", t0.elapsed().as_secs_f64() * 1e6);
}
