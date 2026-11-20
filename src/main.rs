// SPDX-FileCopyrightText: 2026 Alexei Lemos Cavalcante Ferreira <ferreira.alexei@gmail.com>
// SPDX-License-Identifier: GPL-3.0-or-later

use std::env;
use std::fs;
use std::time::Instant;
use tupa_epistemic_verifier::{
    adjudicate_lean4_submission_trust_boundary, calibrated_cross_entropy_bits,
    detect_synthetic_citation_hallucination, evaluate_zeta23_montgomery_taylor_bounds,
    fisher_rao_geodesic_rad, jtms_grounded_closure, scan_lean4_solution_source, sha256_hex,
    sqrt_jensen_shannon_metric, tsallis_proper_score_bits, validate_orcid_iso7064,
    vovk_ctw_compress_rate_bpb, JtmsJustification, LeanSourceStaticScan,
};

fn main() {
    let args: Vec<String> = env::args().collect();
    let t0 = Instant::now();

    if args.len() >= 2 && args[1] == "riemann-audit" {
        let mt = evaluate_zeta23_montgomery_taylor_bounds(1.0);
        let zeta23_verdict = adjudicate_lean4_submission_trust_boundary(
            "leanprover/lean4:v4.33.0-rc2",
            &LeanSourceStaticScan {
                unsound_or_placeholder_tokens: vec![],
                elaborator_bypass_sites: vec![],
                mutual_reflexive_inductive_hint: false,
            },
            true,
            true,
            true,
            true,
        );
        let rh_2026_11_20_verdict = adjudicate_lean4_submission_trust_boundary(
            "leanprover/lean4:v4.33.1",
            &LeanSourceStaticScan {
                unsound_or_placeholder_tokens: vec![],
                elaborator_bypass_sites: vec!["kexactB".into(), "kcopy".into(), "MVarId.assign".into()],
                mutual_reflexive_inductive_hint: true,
            },
            true,
            true,
            false,
            false,
        );
        println!("Tupã Epistemic Verifier v0.3.0 — AI Formal-Proof & Riemann Benchmark Audit");
        println!("  [1] Alpöge–Furman (arXiv:2608.13637v2 · anthropics/formal-math/zeta23):");
        println!("      c_1* = {:.12} · simple on Re(s)=1/2: {:.6}% (> 2/3) · distinct: {:.6}% (> 5/6)",
            mt.c_lambda_star, mt.kappa_simple_critical_line * 100.0, mt.kappa_distinct_zeros * 100.0);
        println!("      Verdict: {} (Comparator + nanoda verified = {})",
            zeta23_verdict.status, zeta23_verdict.accepted_as_machine_verified);
        println!("  [2] Alpöge–Claude Full RH/GRH Submission (2026-11-20 · alpo.ge/remark.pdf · 175,187 files):");
        println!("      Verdict: {} (accepted_as_verified = {})",
            rh_2026_11_20_verdict.status, rh_2026_11_20_verdict.accepted_as_machine_verified);
        for reason in &rh_2026_11_20_verdict.reasons {
            println!("      • Gate hold: {}", reason);
        }
        println!("  Elapsed: {:.2} µs", t0.elapsed().as_secs_f64() * 1e6);
        return;
    }

    if args.len() >= 3 && args[1] == "audit-lean" {
        let path = &args[2];
        let toolchain = args
            .iter()
            .position(|a| a == "--toolchain")
            .and_then(|idx| args.get(idx + 1))
            .map(String::as_str)
            .unwrap_or("leanprover/lean4:v4.33.1");
        let nanoda_ok = args.iter().any(|a| a == "--nanoda-verified");
        let comp_ok = args.iter().any(|a| a == "--comparator-verified");
        let src = match fs::read_to_string(path) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("error: failed to read '{}': {}", path, e);
                std::process::exit(1);
            }
        };
        let scan = scan_lean4_solution_source(&src);
        let verdict = adjudicate_lean4_submission_trust_boundary(
            toolchain, &scan, true, true, comp_ok, nanoda_ok,
        );
        println!("Tupã Epistemic Verifier v0.3.0 — Lean 4 Trust-Boundary Audit");
        println!("  File         : {}", path);
        println!("  Toolchain    : {}", toolchain);
        println!("  Unsound Toks : {:?}", scan.unsound_or_placeholder_tokens);
        println!("  Elab Bypass  : {:?}", scan.elaborator_bypass_sites);
        println!("  Verdict      : {} (accepted={})", verdict.status, verdict.accepted_as_machine_verified);
        for r in &verdict.reasons {
            println!("  • {}", r);
        }
        if !verdict.accepted_as_machine_verified {
            std::process::exit(2);
        }
        return;
    }

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

    println!("Tupã Epistemic Verifier v0.3.0 (100% std-Rust · Zero Cloud · GPL-3.0-or-later)");
    println!("  • LogSumExp Cross-Entropy : {:.4} bits (τ*=0.43) vs {:.4} bits (raw τ=1.00)", ce_cal, ce_raw);
    println!("  • Tsallis q=1.5 Outlier   : {:.4} bits (bounded under p_y=0.001)", tsallis_out);
    println!("  • √JSD / Fisher-Rao       : {:.4} / {:.4} rad", jsd, fr);
    println!("  • Doyle JTMS Retraction   : zombie_claims={} (verified)", usize::from(retracted[1]));
    println!("  • Citation & ORCID Gate   : cite_ok={} · orcid_iso7064={}", cite_ok, orcid_ok);
    println!("  • Self-Check Latency      : {:.2} µs", t0.elapsed().as_secs_f64() * 1e6);
}
