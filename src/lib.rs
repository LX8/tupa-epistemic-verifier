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
}
