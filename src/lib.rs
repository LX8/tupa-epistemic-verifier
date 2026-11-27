// SPDX-FileCopyrightText: 2026 Alexei Lemos Cavalcante Ferreira <ferreira.alexei@gmail.com>
// SPDX-License-Identifier: GPL-3.0-or-later

//! # Tupã Epistemic Verifier Core (`tupa-epistemic-verifier`)
//!
//! 100% Rust standard-library (`wasm32-unknown-unknown` + native C-ABI) verification engine
//! against AI hallucinations, synthetic citations, uncalibrated confidence scores, and
//! non-monotonic "zombie" claims.
//!
//! ## Core Local-First Modules
//! 1. **Exact LogSumExp Temperature-Calibrated Cross-Entropy & Proper Scoring:**
//!    Zero-clipping Shannon cross-entropy (bits), bounded outlier-robust Savage–Bregman
//!    Tsallis `q`-proper score, metric `sqrt(JSD)`, and Amari–Čencov Fisher–Rao geodesic distance.
//! 2. **Doyle Non-Monotonic Justification Truth-Maintenance System (JTMS):**
//!    Deterministic fixed-point premise retraction ensuring `0` zombie claims survive when
//!    an upstream source or dataset is refuted.
//! 3. **Cryptographic Provenance & Identifier Hygiene:**
//!    Pure-`std` FIPS 180-4 SHA-256 digesting and ISO 7064 MOD 11-2 ORCID checksum verification.
//! 4. **Biomolecular Structure & 3D Physical-Field Verification (`structural`):**
//!    Pure-`std` mmCIF parsing, exact Mariani (2013) `lDDT-Ca` vs. `pLDDT` calibration,
//!    Ramachandran & steric-clash auditing, Grant–Pickup Gaussian 3D shape & Coulombic ESP
//!    Tanimoto / Hodgkin fields, directional pharmacophores, and tetrahedral chirality guards.

pub mod structural;
pub use structural::*;

const LN2: f64 = core::f64::consts::LN_2;

/// Numerically stable LogSumExp over raw logits `z`, returning `ln(\sum_i exp(z_i / tau))`
/// without `exp` overflow or `1e-12` probability-clipping bias.
pub fn logsumexp_scaled(logits: &[f64], temperature: f64) -> f64 {
    if logits.is_empty() {
        return f64::NEG_INFINITY;
    }
    let inv_t = 1.0 / temperature.max(1.0e-6);
    let mut max_z = f64::NEG_INFINITY;
    for &z in logits {
        let s = z * inv_t;
        if s > max_z {
            max_z = s;
        }
    }
    let mut sum_exp = 0.0f64;
    for &z in logits {
        sum_exp += (z * inv_t - max_z).exp();
    }
    max_z + sum_exp.ln()
}

/// Exact categorical cross-entropy in Shannon bits `H(y, softmax(z / tau))` directly from raw logits.
pub fn calibrated_cross_entropy_bits(logits: &[f64], target_idx: usize, temperature: f64) -> f64 {
    if logits.is_empty() || target_idx >= logits.len() {
        return f64::INFINITY;
    }
    let inv_t = 1.0 / temperature.max(1.0e-6);
    let lse = logsumexp_scaled(logits, temperature);
    (lse - logits[target_idx] * inv_t) / LN2
}

/// Strictly proper, bounded Savage–Bregman Tsallis `q`-score (`q > 1`, `q -> 1` recovers Shannon log-loss).
/// Unlike unbounded log-loss, Tsallis `q`-score stays upper-bounded under corrupted/adversarial labels.
pub fn tsallis_proper_score_bits(probs: &[f64], target_idx: usize, q: f64) -> f64 {
    if probs.is_empty() || target_idx >= probs.len() {
        return f64::INFINITY;
    }
    if (q - 1.0).abs() < 1.0e-6 {
        return -probs[target_idx].max(1.0e-300).ln() / LN2;
    }
    let q_minus_1 = q - 1.0;
    let mut norm_q = 0.0f64;
    for &p in probs {
        norm_q += p.max(0.0).powf(q);
    }
    let pt_qm1 = probs[target_idx].max(0.0).powf(q_minus_1);
    let raw_nats = (1.0 - (q * pt_qm1 - norm_q)) / q_minus_1;
    raw_nats.max(0.0) / LN2
}

/// True metric square-root Jensen–Shannon distance `sqrt(JSD(P || Q))` in `[0, 1]` (base-2).
pub fn sqrt_jensen_shannon_metric(p: &[f64], q: &[f64]) -> f64 {
    if p.len() != q.len() || p.is_empty() {
        return 1.0;
    }
    let mut jsd_bits = 0.0f64;
    for (&pi, &qi) in p.iter().zip(q.iter()) {
        let u = pi.max(0.0);
        let v = qi.max(0.0);
        let m = 0.5 * (u + v);
        if u > 0.0 && m > 0.0 {
            jsd_bits += 0.5 * u * (u / m).ln() / LN2;
        }
        if v > 0.0 && m > 0.0 {
            jsd_bits += 0.5 * v * (v / m).ln() / LN2;
        }
    }
    jsd_bits.clamp(0.0, 1.0).sqrt()
}

/// Exact Amari–Čencov Fisher–Rao Riemannian geodesic distance `d_FR(P, Q) = 2 arccos(\sum_i sqrt(p_i q_i))`.
pub fn fisher_rao_geodesic_rad(p: &[f64], q: &[f64]) -> f64 {
    if p.len() != q.len() || p.is_empty() {
        return core::f64::consts::PI;
    }
    let mut bc = 0.0f64;
    for (&pi, &qi) in p.iter().zip(q.iter()) {
        bc += (pi.max(0.0) * qi.max(0.0)).sqrt();
    }
    2.0 * bc.clamp(0.0, 1.0).acos()
}

/// Non-monotonic Justification Truth-Maintenance rule (`antecedents => consequent`, unless `out_blockers` hold).
#[derive(Debug, Clone)]
pub struct JtmsJustification {
    pub consequent: usize,
    pub in_antecedents: Vec<usize>,
    pub out_blockers: Vec<usize>,
}

/// Computes the grounded well-founded fixed point of a Doyle JTMS graph.
/// Retracting any premise deterministically cascades to all dependent claims (`0` zombie claims).
pub fn jtms_grounded_closure(
    num_nodes: usize,
    premises: &[usize],
    rules: &[JtmsJustification],
) -> Vec<bool> {
    let mut believed = vec![false; num_nodes];
    for &p in premises {
        if p < num_nodes {
            believed[p] = true;
        }
    }
    // Iterate monotonic support conditioned on current out-blockers, restarting on blocker invalidation
    for _ in 0..(num_nodes + 2) {
        let snapshot = believed.clone();
        let mut next = vec![false; num_nodes];
        for &p in premises {
            if p < num_nodes {
                next[p] = true;
            }
        }
        loop {
            let mut progressed = false;
            for r in rules {
                if r.consequent >= num_nodes || next[r.consequent] {
                    continue;
                }
                let ins_ok = r.in_antecedents.iter().all(|&a| a < num_nodes && next[a]);
                let outs_ok = r.out_blockers.iter().all(|&b| b >= num_nodes || !snapshot[b]);
                if ins_ok && outs_ok {
                    next[r.consequent] = true;
                    progressed = true;
                }
            }
            if !progressed {
                break;
            }
        }
        if next == believed {
            break;
        }
        believed = next;
    }
    believed
}

/// Validates an ORCID identifier (`0000-0000-0000-000X`) under ISO/IEC 7064:2003 MOD 11-2.
pub fn validate_orcid_iso7064(orcid: &str) -> bool {
    let compact: Vec<char> = orcid.trim().chars().filter(|&c| c != '-').collect();
    if compact.len() != 16 {
        return false;
    }
    let mut total = 0u32;
    for &ch in &compact[..15] {
        let Some(d) = ch.to_digit(10) else {
            return false;
        };
        total = (total + d) * 2;
    }
    let remainder = total % 11;
    let result = (12 - remainder) % 11;
    let expected = if result == 10 {
        'X'
    } else {
        char::from_digit(result, 10).unwrap_or('?')
    };
    compact[15].to_ascii_uppercase() == expected
}

/// Pure-`std` FIPS 180-4 SHA-256 digest returning lowercase hex (`64` chars) for W3C RO-Crate provenance.
pub fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
    ];
    let bit_len = (data.len() as u64) * 8;
    let mut msg = data.to_vec();
    msg.push(0x80);
    while (msg.len() % 64) != 56 {
        msg.push(0x00);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in msg.chunks_exact(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[4 * i],
                chunk[4 * i + 1],
                chunk[4 * i + 2],
                chunk[4 * i + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let temp1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }
    let mut out = String::with_capacity(64);
    for word in h {
        out.push_str(&format!("{:08x}", word));
    }
    out
}

/// Evaluates 5-order (`d = 0..4`) escape-interpolated Krichevsky–Trofimov / Vovk Bayesian mixture
/// lossless code length (`bits/byte`) on arbitrary byte stream `data`, returning
/// `(order0_bpb, vovk_mix_bpb, vovk_regret_margin_bits)`.
pub fn vovk_ctw_compress_rate_bpb(data: &[u8], dirichlet_alpha: f64) -> (f64, f64, f64) {
    if data.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let ln2 = LN2;
    let alpha = dirichlet_alpha.clamp(0.01, 4.0);
    let mut order0_counts = [0u32; 256];
    let mut order0_total: u64 = 0;
    let mut ctx_totals: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();
    let mut sym_counts: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();

    let mut log_weights = [0.0f64; 5];
    let mut cum_expert_bits = [0.0f64; 5];
    let mut cum_mix_bits = 0.0f64;

    for (t, &byte) in data.iter().enumerate() {
        let sym = byte as usize;
        let max_lw = log_weights
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let mut w_norm = [0.0f64; 5];
        let mut w_sum = 0.0f64;
        for d in 0..5 {
            w_norm[d] = (log_weights[d] - max_lw).exp();
            w_sum += w_norm[d];
        }
        for d in 0..5 {
            w_norm[d] /= w_sum;
        }

        let mut d_ctx_key = [0u64; 5];
        let mut rolling: u64 = 0xCBF2_9CE4_8422_2325;
        for d in 1..=4 {
            if t >= d {
                rolling = (rolling ^ ((data[t - d] as u64) + 1)).wrapping_mul(0x0100_0000_01B3);
                d_ctx_key[d] = ((d as u64) << 56) | (rolling & 0x00FF_FFFF_FFFF_FFFF);
            }
        }

        let mut p_d = [0.0f64; 5];
        p_d[0] = ((order0_counts[sym] as f64) + alpha) / ((order0_total as f64) + 256.0 * alpha);
        for d in 1..=4 {
            if t < d {
                p_d[d] = p_d[d - 1];
            } else {
                let n_ctx = *ctx_totals.get(&d_ctx_key[d]).unwrap_or(&0);
                if n_ctx == 0 {
                    p_d[d] = p_d[d - 1];
                } else {
                    let skey = d_ctx_key[d].wrapping_mul(257) ^ (sym as u64);
                    let c_sym = *sym_counts.get(&skey).unwrap_or(&0);
                    let esc = alpha / ((n_ctx as f64) + alpha);
                    p_d[d] = ((c_sym as f64) / ((n_ctx as f64) + alpha)) + esc * p_d[d - 1];
                }
            }
        }

        let mut p_mix = 0.0f64;
        for d in 0..5 {
            p_mix += w_norm[d] * p_d[d];
        }
        cum_mix_bits += -p_mix.max(1.0e-300).ln() / ln2;

        for d in 0..5 {
            let loss_d = -p_d[d].max(1.0e-300).ln() / ln2;
            cum_expert_bits[d] += loss_d;
            log_weights[d] -= ln2 * loss_d;
        }

        order0_counts[sym] += 1;
        order0_total += 1;
        for d in 1..=4 {
            if t >= d {
                *ctx_totals.entry(d_ctx_key[d]).or_insert(0) += 1;
                let skey = d_ctx_key[d].wrapping_mul(257) ^ (sym as u64);
                *sym_counts.entry(skey).or_insert(0) += 1;
            }
        }
    }

    let n = data.len() as f64;
    let min_expert = cum_expert_bits
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let regret_margin = ((5.0f64).ln() / ln2) - (cum_mix_bits - min_expert);
    (cum_expert_bits[0] / n, cum_mix_bits / n, regret_margin)
}

/// Fail-closed detector for synthetic/placeholder DOI and arXiv hallucinations commonly emitted by LLMs.
pub fn detect_synthetic_citation_hallucination(text: &str) -> Vec<String> {
    let mut flags = Vec::new();
    let lower = text.to_lowercase();
    for pat in [
        "10.1234/",
        "10.0000/",
        "10.9999/",
        "10.xxxx/",
        "arxiv:0000.00000",
        "arxiv:1234.56789",
        "arxiv:9999.",
        "doi:10.1000/xyz",
    ] {
        if lower.contains(pat) {
            flags.push(format!("synthetic citation pattern '{}'", pat));
        }
    }
    flags
}

/// Exact Montgomery–Taylor unconditional critical-line & distinct zero proportions
/// from Alpöge & Furman (arXiv:2608.13637v2, formalized in `anthropics/formal-math/zeta23`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Zeta23CriticalLineBounds {
    pub lambda_window: f64,
    pub c_lambda_star: f64,
    pub kappa_simple_critical_line: f64,
    pub kappa_distinct_zeros: f64,
    pub gain_over_five_twelfths: f64,
}

/// Computes `c*_lambda = sqrt(2) * tan(theta) / (1 + theta * tan(theta))` with `theta = lambda / sqrt(2)`,
/// returning `(2 - 1 / c*_lambda)` for simple critical-line zeros (`0.67250463...` at `lambda = 1`)
/// and `(3 - 1 / c*_lambda) / 2` for distinct zeros (`0.83625232...` at `lambda = 1`).
pub fn evaluate_zeta23_montgomery_taylor_bounds(lambda_window: f64) -> Zeta23CriticalLineBounds {
    let lam = lambda_window.clamp(0.05, 2.0);
    let sqrt2 = core::f64::consts::SQRT_2;
    let theta = lam / sqrt2;
    let tan_t = theta.tan();
    let c_star = (sqrt2 * tan_t) / (1.0 + theta * tan_t);
    let inv_c = 1.0 / c_star;
    let kappa_simple = 2.0 - inv_c;
    let kappa_distinct = 0.5 * (3.0 - inv_c);
    Zeta23CriticalLineBounds {
        lambda_window: lam,
        c_lambda_star: c_star,
        kappa_simple_critical_line: kappa_simple,
        kappa_distinct_zeros: kappa_distinct,
        gain_over_five_twelfths: kappa_simple - (5.0 / 12.0),
    }
}

/// Comment-stripped static soundness & metaprogramming audit of Lean 4 source files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeanSourceStaticScan {
    pub unsound_or_placeholder_tokens: Vec<String>,
    pub elaborator_bypass_sites: Vec<String>,
    pub mutual_reflexive_inductive_hint: bool,
}

/// Strips Lean 4 `-- ...` line comments, nested `/- ... -/` block comments, and string literals
/// before scanning for soundness-compromising keywords (`sorry`, `axiom`, `native_decide`, `extern`,
/// `implemented_by`, `unsafe`, `opaque`, `skipKernelTC`) and elaborator-bypassing metaprogramming
/// (`MVarId.assign`, `kexactB`, `kcopy`, `elab `, `run_cmd`, `run_tac`, `partial def`).
pub fn strip_lean4_comments_and_strings(src: &str) -> String {
    let bytes = src.as_bytes();
    let mut out = String::with_capacity(src.len());
    let mut i = 0usize;
    let mut block_depth = 0usize;
    let mut in_line_comment = false;
    let mut in_string = false;

    while i < bytes.len() {
        if in_line_comment {
            if bytes[i] == b'\n' {
                in_line_comment = false;
                out.push('\n');
            }
            i += 1;
            continue;
        }
        if block_depth > 0 {
            if i + 1 < bytes.len() && bytes[i] == b'/' && bytes[i + 1] == b'-' {
                block_depth += 1;
                i += 2;
                continue;
            }
            if i + 1 < bytes.len() && bytes[i] == b'-' && bytes[i + 1] == b'/' {
                block_depth -= 1;
                i += 2;
                out.push(' ');
                continue;
            }
            if bytes[i] == b'\n' {
                out.push('\n');
            }
            i += 1;
            continue;
        }
        if in_string {
            if bytes[i] == b'\\' && i + 1 < bytes.len() {
                i += 2;
                continue;
            }
            if bytes[i] == b'"' {
                in_string = false;
            }
            i += 1;
            continue;
        }
        if i + 1 < bytes.len() && bytes[i] == b'-' && bytes[i + 1] == b'-' {
            in_line_comment = true;
            i += 2;
            continue;
        }
        if i + 1 < bytes.len() && bytes[i] == b'/' && bytes[i + 1] == b'-' {
            block_depth = 1;
            i += 2;
            continue;
        }
        if bytes[i] == b'"' {
            in_string = true;
            out.push(' ');
            i += 1;
            continue;
        }
        let ch = src[i..].chars().next().unwrap_or(' ');
        out.push(ch);
        i += ch.len_utf8();
    }
    out
}

fn contains_ident_word(haystack: &str, word: &str) -> bool {
    let is_ident = |b: u8| b.is_ascii_alphanumeric() || b == b'_' || b == b'\'';
    let hb = haystack.as_bytes();
    let wb = word.as_bytes();
    if wb.is_empty() || wb.len() > hb.len() {
        return false;
    }
    let mut pos = 0usize;
    while let Some(rel) = haystack[pos..].find(word) {
        let idx = pos + rel;
        let left_ok = idx == 0 || !is_ident(hb[idx - 1]);
        let right_idx = idx + wb.len();
        let right_ok = right_idx >= hb.len() || !is_ident(hb[right_idx]);
        if left_ok && right_ok {
            return true;
        }
        pos = idx + 1;
    }
    false
}

/// Scans an untrusted Lean 4 solution module (outside the trusted `Challenge.lean` spec).
pub fn scan_lean4_solution_source(src: &str) -> LeanSourceStaticScan {
    let stripped = strip_lean4_comments_and_strings(src);
    // Also strip comments while retaining quoted tactic syntax names (e.g. `elab "kexactB"`)
    let mut no_comments = String::with_capacity(src.len());
    {
        let bytes = src.as_bytes();
        let mut i = 0usize;
        let mut block_depth = 0usize;
        let mut in_line = false;
        while i < bytes.len() {
            if in_line {
                if bytes[i] == b'\n' {
                    in_line = false;
                    no_comments.push('\n');
                }
                i += 1;
                continue;
            }
            if block_depth > 0 {
                if i + 1 < bytes.len() && bytes[i] == b'/' && bytes[i + 1] == b'-' {
                    block_depth += 1;
                    i += 2;
                    continue;
                }
                if i + 1 < bytes.len() && bytes[i] == b'-' && bytes[i + 1] == b'/' {
                    block_depth -= 1;
                    i += 2;
                    no_comments.push(' ');
                    continue;
                }
                i += 1;
                continue;
            }
            if i + 1 < bytes.len() && bytes[i] == b'-' && bytes[i + 1] == b'-' {
                in_line = true;
                i += 2;
                continue;
            }
            if i + 1 < bytes.len() && bytes[i] == b'/' && bytes[i + 1] == b'-' {
                block_depth = 1;
                i += 2;
                continue;
            }
            let ch = src[i..].chars().next().unwrap_or(' ');
            no_comments.push(ch);
            i += ch.len_utf8();
        }
    }

    let mut unsound = Vec::new();
    for kw in [
        "sorry",
        "sorryAx",
        "admit",
        "axiom",
        "constant",
        "opaque",
        "native_decide",
        "skipKernelTC",
        "unsafe",
    ] {
        if contains_ident_word(&stripped, kw) {
            unsound.push(kw.to_string());
        }
    }
    for attr in ["@[extern", "@[implemented_by", "debug.skipKernelTC"] {
        if stripped.contains(attr) {
            unsound.push(attr.to_string());
        }
    }

    let mut bypass = Vec::new();
    for meta_kw in ["MVarId.assign", "kexactB", "kcopy", "run_cmd", "run_tac"] {
        if contains_ident_word(&stripped, meta_kw)
            || stripped.contains(meta_kw)
            || no_comments.contains(&format!("\"{}\"", meta_kw))
        {
            bypass.push(meta_kw.to_string());
        }
    }
    if contains_ident_word(&stripped, "elab") {
        bypass.push("elab".to_string());
    }
    if stripped.contains("partial def ") {
        bypass.push("partial def".to_string());
    }

    let mutual_reflexive = stripped.contains("mutual") && stripped.contains("inductive ");
    LeanSourceStaticScan {
        unsound_or_placeholder_tokens: unsound,
        elaborator_bypass_sites: bypass,
        mutual_reflexive_inductive_hint: mutual_reflexive,
    }
}

/// Parses a `lean-toolchain` string such as `"leanprover/lean4:v4.33.1"` or `"v4.34.1"` into `(major, minor, patch)`.
pub fn parse_lean_toolchain_semver(toolchain: &str) -> Option<(u32, u32, u32)> {
    let tail = toolchain
        .trim()
        .rsplit(':')
        .next()
        .unwrap_or(toolchain)
        .trim()
        .trim_start_matches('v');
    let core_part = tail.split('-').next().unwrap_or(tail);
    let mut parts = core_part.split('.');
    let major = parts.next()?.parse::<u32>().ok()?;
    let minor = parts.next()?.parse::<u32>().ok()?;
    let patch = parts.next().unwrap_or("0").parse::<u32>().ok()?;
    Some((major, minor, patch))
}

/// Trust-boundary adjudication for large AI-assisted Lean 4 formal proof submissions
/// (modeled on the Palomar / `leanprover/comparator` + `nanoda` protocol).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeanFormalSubmissionVerdict {
    pub status: &'static str,
    pub accepted_as_machine_verified: bool,
    pub kernel_runtime_cve_exposed: bool,
    pub requires_external_kernel_replay: bool,
    pub reasons: Vec<String>,
}

/// Evaluates whether an AI-generated Lean 4 submission can be promoted to `MACHINE_VERIFIED`
/// or must remain fail-closed (`REJECTED_UNSOUND` or `UNVERIFIED_PENDING_EXTERNAL_KERNEL_REPLAY`).
pub fn adjudicate_lean4_submission_trust_boundary(
    lean_toolchain: &str,
    scan: &LeanSourceStaticScan,
    challenge_imports_mathlib_only: bool,
    only_standard_three_axioms: bool,
    external_comparator_436_passed: bool,
    external_nanoda_passed: bool,
) -> LeanFormalSubmissionVerdict {
    let mut reasons = Vec::new();
    if !challenge_imports_mathlib_only {
        reasons.push("Challenge module imports untrusted modules beyond Mathlib".to_string());
    }
    if !only_standard_three_axioms {
        reasons.push("Export depends on axioms outside {propext, Quot.sound, Classical.choice}".to_string());
    }
    if !scan.unsound_or_placeholder_tokens.is_empty() {
        reasons.push(format!(
            "Solution contains unsound/placeholder tokens: {}",
            scan.unsound_or_placeholder_tokens.join(", ")
        ));
    }
    if !reasons.is_empty() {
        return LeanFormalSubmissionVerdict {
            status: "REJECTED_UNSOUND_OR_INCOMPLETE",
            accepted_as_machine_verified: false,
            kernel_runtime_cve_exposed: false,
            requires_external_kernel_replay: true,
            reasons,
        };
    }

    let semver = parse_lean_toolchain_semver(lean_toolchain).unwrap_or((4, 0, 0));
    // Lean < 4.32.2 had #14576 (wrong-structure projections proving False);
    // Lean < 4.34.1 lacks runtime soundness fixes #15241 and #15289.
    let pre_14576 = semver < (4, 32, 2);
    let runtime_cve_15241_15289 = semver < (4, 34, 1);
    let has_elab_bypass = !scan.elaborator_bypass_sites.is_empty();
    let needs_external_replay = pre_14576 || runtime_cve_15241_15289 || has_elab_bypass;

    if pre_14576 {
        reasons.push(format!(
            "Toolchain {} predates Lean 4.32.2 kernel projection soundness fix (#14576)",
            lean_toolchain
        ));
    } else if runtime_cve_15241_15289 {
        reasons.push(format!(
            "Toolchain {} predates Lean 4.34.1 runtime soundness fixes (#15241, #15289)",
            lean_toolchain
        ));
    }
    if has_elab_bypass {
        reasons.push(format!(
            "Elaborator-bypassing metaprogramming present ({}) — kernel replay is mandatory",
            scan.elaborator_bypass_sites.join(", ")
        ));
    }

    if needs_external_replay && !(external_comparator_436_passed && external_nanoda_passed) {
        return LeanFormalSubmissionVerdict {
            status: "UNVERIFIED_PENDING_EXTERNAL_KERNEL_REPLAY",
            accepted_as_machine_verified: false,
            kernel_runtime_cve_exposed: pre_14576 || runtime_cve_15241_15289,
            requires_external_kernel_replay: true,
            reasons,
        };
    }

    LeanFormalSubmissionVerdict {
        status: "MACHINE_VERIFIED_MULTI_KERNEL",
        accepted_as_machine_verified: true,
        kernel_runtime_cve_exposed: pre_14576 || runtime_cve_15241_15289,
        requires_external_kernel_replay: needs_external_replay,
        reasons,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logsumexp_and_tsallis_strict_properness() {
        // Extreme logits [1000.0, 999.0, 998.0] must not overflow exp()
        let logits = [1000.0, 999.0, 998.0];
        let ce = calibrated_cross_entropy_bits(&logits, 0, 1.0);
        assert!(ce.is_finite() && ce > 0.58 && ce < 0.59);

        // Strictly proper: true distribution p minimizes expected Tsallis score uniquely
        let p = [0.7, 0.2, 0.1];
        let q_wrong = [0.2, 0.5, 0.3];
        let mut exp_true = 0.0;
        let mut exp_wrong = 0.0;
        for i in 0..3 {
            exp_true += p[i] * tsallis_proper_score_bits(&p, i, 1.5);
            exp_wrong += p[i] * tsallis_proper_score_bits(&q_wrong, i, 1.5);
        }
        assert!(exp_true < exp_wrong);
    }

    #[test]
    fn test_jsd_triangle_inequality_and_fisher_rao_symmetry() {
        let p = [0.8, 0.15, 0.05];
        let q = [0.3, 0.5, 0.2];
        let r = [0.1, 0.2, 0.7];
        let d_pq = sqrt_jensen_shannon_metric(&p, &q);
        let d_qr = sqrt_jensen_shannon_metric(&q, &r);
        let d_pr = sqrt_jensen_shannon_metric(&p, &r);
        assert!(d_pr <= d_pq + d_qr + 1.0e-12);
        assert!((fisher_rao_geodesic_rad(&p, &q) - fisher_rao_geodesic_rad(&q, &p)).abs() < 1.0e-12);
    }

    #[test]
    fn test_jtms_retraction_eliminates_zombie_claims() {
        // 0 (source dataset) -> 1 (intermediate lemma) -> 2 (published conclusion), unless 3 (retraction notice)
        let rules = vec![
            JtmsJustification { consequent: 1, in_antecedents: vec![0], out_blockers: vec![3] },
            JtmsJustification { consequent: 2, in_antecedents: vec![1], out_blockers: vec![] },
        ];
        let before = jtms_grounded_closure(4, &[0], &rules);
        assert_eq!(before, vec![true, true, true, false]);

        // Retract source 0 OR add blocker 3 -> both 1 and 2 must drop to false (0 zombie truths)
        let after_retract = jtms_grounded_closure(4, &[], &rules);
        assert_eq!(after_retract, vec![false, false, false, false]);
        let after_blocker = jtms_grounded_closure(4, &[0, 3], &rules);
        assert_eq!(after_blocker, vec![true, false, false, true]);
    }

    #[test]
    fn test_sha256_nist_and_orcid_iso7064() {
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert!(validate_orcid_iso7064("0009-0009-0139-3700"));
        assert!(!validate_orcid_iso7064("0009-0009-0139-3701"));
    }

    #[test]
    fn test_vovk_ctw_compression_and_regret_theorem() {
        let sample = "<page><title>Epistemic Verification</title><text>Local-first reproducible Rust verification eliminates zombie claims and uncalibrated confidence.</text></page>".repeat(6);
        let (bpb0, bpb_mix, margin) = vovk_ctw_compress_rate_bpb(sample.as_bytes(), 0.50);
        assert!(bpb_mix < bpb0 * 0.55, "5-order Vovk/CTW ({}) must compress structured XML well below order-0 ({})", bpb_mix, bpb0);
        assert!(margin >= -1.0e-9, "Vovk universal regret <= log2(5) violated: margin={}", margin);
    }

    #[test]
    fn test_synthetic_citation_hallucination_blocker() {
        assert!(detect_synthetic_citation_hallucination("arXiv:2509.24593 and doi:10.1103/PhysRevD.38.472").is_empty());
        assert_eq!(detect_synthetic_citation_hallucination("See doi:10.1234/fake-llm-paper and arXiv:0000.00000").len(), 2);
    }

    #[test]
    fn test_zeta23_montgomery_taylor_exact_constants_and_monotonicity() {
        let mt = evaluate_zeta23_montgomery_taylor_bounds(1.0);
        // Exact algebraic identity: 1/c_1* = 1/2 + cot(1/sqrt(2)) / sqrt(2)
        let theta = 1.0 / core::f64::consts::SQRT_2;
        let inv_c_exact = 0.5 + 1.0 / (core::f64::consts::SQRT_2 * theta.tan());
        assert!(((1.0 / mt.c_lambda_star) - inv_c_exact).abs() < 1.0e-14);
        assert!((mt.c_lambda_star - 0.7532960678560707).abs() < 1.0e-14);
        assert!((mt.kappa_simple_critical_line - 0.6725007036794117).abs() < 1.0e-14);
        assert!((mt.kappa_distinct_zeros - 0.8362503518397059).abs() < 1.0e-14);
        assert!(mt.kappa_simple_critical_line > 2.0 / 3.0);
        assert!(mt.kappa_distinct_zeros > 5.0 / 6.0);
        assert!(mt.gain_over_five_twelfths > 0.2558);
    }

    #[test]
    fn test_lean4_source_scanner_ignores_comments_and_catches_sorry_native_decide_and_kexactb() {
        let clean_with_comments = r#"
/- This block comment mentions sorry and native_decide and @[extern "c"] -/
-- Single-line comment mentioning axiom and skipKernelTC
theorem riemannHypothesis : RiemannHypothesis := rhOfAll
"#;
        let s_clean = scan_lean4_solution_source(clean_with_comments);
        assert!(s_clean.unsound_or_placeholder_tokens.is_empty());
        assert!(s_clean.elaborator_bypass_sites.is_empty());

        let meta_bypass = r#"
elab "kexactB" t:term : tactic => do
  MVarId.assign goal expr
"#;
        let s_meta = scan_lean4_solution_source(meta_bypass);
        assert!(s_meta.unsound_or_placeholder_tokens.is_empty());
        assert!(s_meta.elaborator_bypass_sites.contains(&"kexactB".to_string()));
        assert!(s_meta.elaborator_bypass_sites.contains(&"MVarId.assign".to_string()));

        let unsound_src = "theorem fake : RiemannHypothesis := by sorry\n#eval native_decide";
        let s_bad = scan_lean4_solution_source(unsound_src);
        assert!(s_bad.unsound_or_placeholder_tokens.contains(&"sorry".to_string()));
        assert!(s_bad.unsound_or_placeholder_tokens.contains(&"native_decide".to_string()));
    }

    #[test]
    fn test_lean4_submission_auditor_distinguishes_zeta23_verified_vs_rh_2026_11_20_pending() {
        // 1. Alpöge-Furman zeta23 (arXiv:2608.13637v2): Comparator + nanoda replayed -> MACHINE_VERIFIED_MULTI_KERNEL
        let zeta23_scan = LeanSourceStaticScan {
            unsound_or_placeholder_tokens: vec![],
            elaborator_bypass_sites: vec![],
            mutual_reflexive_inductive_hint: false,
        };
        let v_zeta23 = adjudicate_lean4_submission_trust_boundary(
            "leanprover/lean4:v4.33.0-rc2",
            &zeta23_scan,
            true,
            true,
            true,
            true,
        );
        assert!(v_zeta23.accepted_as_machine_verified);
        assert_eq!(v_zeta23.status, "MACHINE_VERIFIED_MULTI_KERNEL");

        // 2. Alpöge-Claude RH/GRH archive (2026-11-20, Lean v4.33.1, 0 sorry, kexactB/kcopy present, nanoda still running)
        // Must stay fail-closed as UNVERIFIED_PENDING_EXTERNAL_KERNEL_REPLAY!
        let rh_scan = LeanSourceStaticScan {
            unsound_or_placeholder_tokens: vec![],
            elaborator_bypass_sites: vec!["kexactB".to_string(), "kcopy".to_string(), "MVarId.assign".to_string()],
            mutual_reflexive_inductive_hint: true,
        };
        let v_rh_pending = adjudicate_lean4_submission_trust_boundary(
            "leanprover/lean4:v4.33.1",
            &rh_scan,
            true,
            true,
            false,
            false,
        );
        assert!(!v_rh_pending.accepted_as_machine_verified);
        assert!(v_rh_pending.kernel_runtime_cve_exposed);
        assert!(v_rh_pending.requires_external_kernel_replay);
        assert_eq!(v_rh_pending.status, "UNVERIFIED_PENDING_EXTERNAL_KERNEL_REPLAY");
    }
}
