// SPDX-FileCopyrightText: 2026 Alexei Lemos Cavalcante Ferreira <ferreira.alexei@gmail.com>
// SPDX-License-Identifier: GPL-3.0-or-later

//! Deterministic Biomolecular Structure & 3D Physical-Field Severe-Testing Module.
//!
//! Pure standard-library Rust (`wasm32-unknown-unknown` + native CLI, zero external crates):
//! - Zero-dependency mmCIF `_atom_site` backbone & `HETATM` ligand parsing.
//! - Exact IUPAC $(\phi, \psi)$ backbone dihedrals and Pauling–Corey NeRF helix kinematics.
//! - Superposition-free Mariani et al. (2013) $\mathrm{lDDT\text{-}C}_\alpha$ vs. AlphaFold
//!   `pLDDT` calibration auditing with automatic sequence-offset alignment.
//! - Grant & Pickup (1995) Gaussian 3D volumetric shape Tanimoto, Coulombic electrostatic
//!   potential (ESP) Gaussian Tanimoto + signed Hodgkin index, directional 3D pharmacophore
//!   overlap (`ShEPhERD` / `shepherd-score` physics), $SO(3)$-invariant signed tetrahedral
//!   chirality guard (PoseBusters failure mode), and protein–ligand steric clash auditing.
//! - Strict epistemic demarcation: software geometry and calibration checks NEVER establish
//!   dynamic conformational completeness, point-mutation $\Delta\Delta G$, or clinical efficacy.

use std::collections::{BTreeMap, BTreeSet};

/// Single protein backbone residue parsed from mmCIF or constructed synthetically.
#[derive(Clone, Debug)]
pub struct ProteinBackboneResidue {
    pub seq_id: i32,
    pub comp_id: [u8; 3],
    pub n: Option<[f64; 3]>,
    pub ca: [f64; 3],
    pub c: Option<[f64; 3]>,
    /// Per-residue `pLDDT` (0..100) from AlphaFold mmCIF B-factor column, or isotropic B-factor.
    pub plddt: f64,
}

/// Deterministic biophysical and confidence-calibration audit of a protein backbone.
#[derive(Clone, Debug)]
pub struct ProteinPhysicalAudit {
    pub n_residues: usize,
    pub mean_plddt: f64,
    pub frac_very_high: f64,
    pub frac_confident: f64,
    pub frac_low: f64,
    pub frac_disordered: f64,
    pub ca_bond_breaks: usize,
    pub cis_peptide_candidates: usize,
    pub steric_clashes_ca: usize,
    pub ramachandran_evaluated: usize,
    pub ramachandran_outliers: usize,
    pub low_plddt_evaluated: usize,
    pub low_plddt_alpha_helical: usize,
    pub helical_fraction_in_low_plddt: f64,
}

/// Comparison result between a predicted structure (e.g. AlphaFold DB) and an experimental
/// reference structure (RCSB PDB) on sequence-aligned backbone residues.
#[derive(Clone, Debug)]
pub struct ProteinModelVsExperiment {
    pub matched_residues: usize,
    pub seq_identity: f64,
    pub lddt_ca: f64,
    pub drmsd_ca: f64,
    pub matched_mean_plddt: f64,
    pub signed_calibration_bias: f64,
    pub abs_calibration_error: f64,
}

/// Heavy atom with 3D coordinates (Angstrom), Bondi van der Waals radius (A), and charge (e).
#[derive(Clone, Debug)]
pub struct MolAtom3D {
    pub pos: [f64; 3],
    pub vdw_radius: f64,
    pub partial_charge: f64,
}

/// Directional 3D pharmacophore feature (`ShEPhERD` / `shepherd-score` formulation).
/// `kind`: 0=Hydrophobic, 1=HBondDonor, 2=HBondAcceptor, 3=Aromatic (plane normal), 4=PosIon, 5=NegIon.
#[derive(Clone, Debug)]
pub struct Pharmacophore3D {
    pub kind: u8,
    pub pos: [f64; 3],
    pub dir: Option<[f64; 3]>,
}

/// Severe-testing self-check receipt enforcing Mayo's power criterion `S_Mayo = 1.0000`.
#[derive(Clone, Debug)]
pub struct StructuralSevereTestReceipt {
    pub positive_control_ok: bool,
    pub mutants_killed: usize,
    pub total_mutants: usize,
    pub mayo_severity_power: f64,
    pub epistemic_demarcation_locked: bool,
    pub passed: bool,
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_norm(vec3_sub(a, b))
}

/// Exact IUPAC dihedral angle (in degrees, `[-180, 180]`) defined by four 3D points.
pub fn protein_dihedral_deg(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let b1 = vec3_sub(p1, p0);
    let b2 = vec3_sub(p2, p1);
    let b3 = vec3_sub(p3, p2);
    let b2_len = vec3_norm(b2);
    if b2_len < 1.0e-12 {
        return 0.0;
    }
    let b2_hat = [b2[0] / b2_len, b2[1] / b2_len, b2[2] / b2_len];
    let n1 = vec3_cross(b1, b2);
    let n2 = vec3_cross(b2, b3);
    let m1 = vec3_cross(n1, b2_hat);
    let x = vec3_dot(n1, n2);
    let y = vec3_dot(m1, n2);
    (-y.atan2(x)).to_degrees()
}

fn split_mmcif_tokens(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let len = bytes.len();
    let mut out = Vec::with_capacity(24);
    let mut i = 0usize;
    while i < len {
        while i < len && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= len || bytes[i] == b'#' {
            break;
        }
        if bytes[i] == b'\'' || bytes[i] == b'"' {
            let quote = bytes[i];
            i += 1;
            let start = i;
            while i < len {
                if bytes[i] == quote && (i + 1 == len || bytes[i + 1].is_ascii_whitespace()) {
                    break;
                }
                i += 1;
            }
            out.push(&line[start..i]);
            if i < len {
                i += 1;
            }
        } else {
            let start = i;
            while i < len && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            out.push(&line[start..i]);
        }
    }
    out
}

fn comp_id_bytes(s: &str) -> [u8; 3] {
    let b = s.as_bytes();
    [
        if !b.is_empty() { b[0].to_ascii_uppercase() } else { b'X' },
        if b.len() > 1 { b[1].to_ascii_uppercase() } else { b'X' },
        if b.len() > 2 { b[2].to_ascii_uppercase() } else { b'X' },
    ]
}

/// Zero-dependency pure-Rust mmCIF `_atom_site.` backbone parser (extracts `N`, `CA`, `C`, `pLDDT`/B).
pub fn parse_mmcif_backbone_chain(
    mmcif_text: &str,
    preferred_chain: Option<&str>,
) -> Vec<ProteinBackboneResidue> {
    #[derive(Default)]
    struct PartialRes {
        comp_id: [u8; 3],
        n: Option<[f64; 3]>,
        ca: Option<[f64; 3]>,
        c: Option<[f64; 3]>,
        plddt: f64,
    }

    let mut in_loop = false;
    let mut reading_atom_headers = false;
    let mut in_atom_rows = false;
    let mut headers: Vec<&str> = Vec::new();

    let mut idx_group = None;
    let mut idx_atom = None;
    let mut idx_comp = None;
    let mut idx_label_asym = None;
    let mut idx_auth_asym = None;
    let mut idx_label_seq = None;
    let mut idx_auth_seq = None;
    let mut idx_alt = None;
    let mut idx_x = None;
    let mut idx_y = None;
    let mut idx_z = None;
    let mut idx_b = None;
    let mut idx_model = None;

    let mut locked_chain: Option<String> = preferred_chain.map(str::to_string);
    let mut residues: BTreeMap<i32, PartialRes> = BTreeMap::new();

    for raw_line in mmcif_text.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if line == "loop_" {
            if in_atom_rows {
                break;
            }
            in_loop = true;
            reading_atom_headers = false;
            headers.clear();
            continue;
        }
        if in_loop && line.starts_with('_') {
            if line.starts_with("_atom_site.") {
                reading_atom_headers = true;
                headers.push(line);
            } else if reading_atom_headers {
                in_loop = false;
                reading_atom_headers = false;
                headers.clear();
            } else {
                in_loop = false;
            }
            continue;
        }
        if reading_atom_headers && !line.starts_with('_') {
            reading_atom_headers = false;
            in_loop = false;
            in_atom_rows = true;
            for (i, h) in headers.iter().enumerate() {
                match *h {
                    "_atom_site.group_PDB" => idx_group = Some(i),
                    "_atom_site.label_atom_id" => idx_atom = Some(i),
                    "_atom_site.label_comp_id" => idx_comp = Some(i),
                    "_atom_site.label_asym_id" => idx_label_asym = Some(i),
                    "_atom_site.auth_asym_id" => idx_auth_asym = Some(i),
                    "_atom_site.label_seq_id" => idx_label_seq = Some(i),
                    "_atom_site.auth_seq_id" => idx_auth_seq = Some(i),
                    "_atom_site.label_alt_id" => idx_alt = Some(i),
                    "_atom_site.Cartn_x" => idx_x = Some(i),
                    "_atom_site.Cartn_y" => idx_y = Some(i),
                    "_atom_site.Cartn_z" => idx_z = Some(i),
                    "_atom_site.B_iso_or_equiv" => idx_b = Some(i),
                    "_atom_site.pdbx_PDB_model_num" => idx_model = Some(i),
                    _ => {}
                }
            }
        }
        if in_atom_rows {
            if line.starts_with('#') || line.starts_with("data_") || line.starts_with('_') {
                break;
            }
            let toks = split_mmcif_tokens(line);
            if toks.len() < headers.len() {
                continue;
            }
            if let Some(ig) = idx_group {
                if toks[ig] != "ATOM" {
                    continue;
                }
            }
            if let Some(im) = idx_model {
                if toks[im] != "1" {
                    continue;
                }
            }
            if let Some(ia) = idx_alt {
                let alt = toks[ia];
                if alt != "." && alt != "?" && alt != "A" {
                    continue;
                }
            }
            let atom_name = match idx_atom {
                Some(i) => toks[i],
                None => continue,
            };
            if atom_name != "N" && atom_name != "CA" && atom_name != "C" {
                continue;
            }
            let auth_ch = idx_auth_asym.map(|i| toks[i]).unwrap_or("");
            let label_ch = idx_label_asym.map(|i| toks[i]).unwrap_or("");
            if let Some(ref want) = locked_chain {
                if auth_ch != want.as_str() && label_ch != want.as_str() {
                    continue;
                }
            } else {
                let first = if !auth_ch.is_empty() && auth_ch != "." && auth_ch != "?" {
                    auth_ch
                } else {
                    label_ch
                };
                locked_chain = Some(first.to_string());
            }

            let seq_id = idx_auth_seq
                .and_then(|i| toks[i].parse::<i32>().ok())
                .or_else(|| idx_label_seq.and_then(|i| toks[i].parse::<i32>().ok()));
            let Some(seq_id) = seq_id else { continue };
            let (Some(ix), Some(iy), Some(iz)) = (idx_x, idx_y, idx_z) else {
                continue;
            };
            let (Ok(x), Ok(y), Ok(z)) = (
                toks[ix].parse::<f64>(),
                toks[iy].parse::<f64>(),
                toks[iz].parse::<f64>(),
            ) else {
                continue;
            };
            let b_val = idx_b
                .and_then(|i| toks[i].parse::<f64>().ok())
                .unwrap_or(0.0);
            let comp = idx_comp.map(|i| comp_id_bytes(toks[i])).unwrap_or(*b"UNK");

            let entry = residues.entry(seq_id).or_insert_with(|| PartialRes {
                comp_id: comp,
                ..Default::default()
            });
            match atom_name {
                "N" => {
                    if entry.n.is_none() {
                        entry.n = Some([x, y, z]);
                    }
                }
                "CA" => {
                    if entry.ca.is_none() {
                        entry.ca = Some([x, y, z]);
                        entry.plddt = b_val;
                    }
                }
                "C" => {
                    if entry.c.is_none() {
                        entry.c = Some([x, y, z]);
                    }
                }
                _ => {}
            }
        }
    }

    residues
        .into_iter()
        .filter_map(|(seq_id, p)| {
            p.ca.map(|ca| ProteinBackboneResidue {
                seq_id,
                comp_id: p.comp_id,
                n: p.n,
                ca,
                c: p.c,
                plddt: p.plddt,
            })
        })
        .collect()
}

/// Deterministic biophysical and confidence-calibration audit of a protein backbone.
pub fn audit_protein_backbone(residues: &[ProteinBackboneResidue]) -> ProteinPhysicalAudit {
    let n = residues.len();
    if n == 0 {
        return ProteinPhysicalAudit {
            n_residues: 0,
            mean_plddt: 0.0,
            frac_very_high: 0.0,
            frac_confident: 0.0,
            frac_low: 0.0,
            frac_disordered: 0.0,
            ca_bond_breaks: 0,
            cis_peptide_candidates: 0,
            steric_clashes_ca: 0,
            ramachandran_evaluated: 0,
            ramachandran_outliers: 0,
            low_plddt_evaluated: 0,
            low_plddt_alpha_helical: 0,
            helical_fraction_in_low_plddt: 0.0,
        };
    }

    let mut sum_p = 0.0f64;
    let (mut n_vh, mut n_conf, mut n_low, mut n_dis) = (0usize, 0usize, 0usize, 0usize);
    for r in residues {
        sum_p += r.plddt;
        if r.plddt >= 90.0 {
            n_vh += 1;
        } else if r.plddt >= 70.0 {
            n_conf += 1;
        } else if r.plddt >= 50.0 {
            n_low += 1;
        } else {
            n_dis += 1;
        }
    }

    let mut ca_bond_breaks = 0usize;
    let mut cis_peptide_candidates = 0usize;
    for i in 0..(n - 1) {
        if residues[i + 1].seq_id == residues[i].seq_id + 1 {
            let d = vec3_dist(residues[i].ca, residues[i + 1].ca);
            if d > 4.2 || d < 2.7 {
                ca_bond_breaks += 1;
            } else if d <= 3.2 {
                cis_peptide_candidates += 1;
            }
        }
    }

    let mut steric_clashes_ca = 0usize;
    for i in 0..n {
        let ca_i = residues[i].ca;
        for j in (i + 2)..n {
            let ca_j = residues[j].ca;
            let dx = ca_i[0] - ca_j[0];
            if dx.abs() >= 3.0 {
                continue;
            }
            let dy = ca_i[1] - ca_j[1];
            if dy.abs() >= 3.0 {
                continue;
            }
            let dz = ca_i[2] - ca_j[2];
            if dx * dx + dy * dy + dz * dz < 9.0 {
                steric_clashes_ca += 1;
            }
        }
    }

    let (mut rama_eval, mut rama_out, mut low_eval, mut low_helix) = (0usize, 0usize, 0usize, 0usize);
    for i in 1..(n - 1) {
        if residues[i].seq_id != residues[i - 1].seq_id + 1
            || residues[i + 1].seq_id != residues[i].seq_id + 1
        {
            continue;
        }
        let (Some(c_prev), Some(n_cur), Some(c_cur), Some(n_next)) = (
            residues[i - 1].c,
            residues[i].n,
            residues[i].c,
            residues[i + 1].n,
        ) else {
            continue;
        };
        let phi = protein_dihedral_deg(c_prev, n_cur, residues[i].ca, c_cur);
        let psi = protein_dihedral_deg(n_cur, residues[i].ca, c_cur, n_next);
        rama_eval += 1;

        if &residues[i].comp_id != b"GLY" {
            let in_left_alpha = (25.0..=105.0).contains(&phi) && (-30.0..=90.0).contains(&psi);
            if phi > 15.0 && phi < 165.0 && !in_left_alpha {
                rama_out += 1;
            }
        }
        if residues[i].plddt < 70.0 {
            low_eval += 1;
            if (-95.0..=-35.0).contains(&phi) && (-70.0..=-15.0).contains(&psi) {
                low_helix += 1;
            }
        }
    }

    let nf = n as f64;
    ProteinPhysicalAudit {
        n_residues: n,
        mean_plddt: sum_p / nf,
        frac_very_high: n_vh as f64 / nf,
        frac_confident: n_conf as f64 / nf,
        frac_low: n_low as f64 / nf,
        frac_disordered: n_dis as f64 / nf,
        ca_bond_breaks,
        cis_peptide_candidates,
        steric_clashes_ca,
        ramachandran_evaluated: rama_eval,
        ramachandran_outliers: rama_out,
        low_plddt_evaluated: low_eval,
        low_plddt_alpha_helical: low_helix,
        helical_fraction_in_low_plddt: if low_eval > 0 {
            low_helix as f64 / low_eval as f64
        } else {
            0.0
        },
    }
}

/// Sequence-offset alignment and exact superposition-free Mariani et al. (2013) `lDDT-Ca`
/// (`R_0 = 15.0 A`, thresholds `0.5 / 1.0 / 2.0 / 4.0 A`) plus `dRMSD-Ca`.
pub fn compare_model_vs_pdb_lddt_ca(
    model: &[ProteinBackboneResidue],
    reference: &[ProteinBackboneResidue],
    inclusion_radius: f64,
) -> Option<ProteinModelVsExperiment> {
    if model.is_empty() || reference.is_empty() {
        return None;
    }
    let model_by_seq: BTreeMap<i32, usize> = model
        .iter()
        .enumerate()
        .map(|(idx, r)| (r.seq_id, idx))
        .collect();

    let mut candidate_deltas: BTreeSet<i32> = BTreeSet::new();
    candidate_deltas.insert(0);
    if model.len() >= 3 && reference.len() >= 3 {
        let mut ref_3mers: BTreeMap<[u8; 9], Vec<i32>> = BTreeMap::new();
        for w in reference.windows(3) {
            if w[1].seq_id == w[0].seq_id + 1 && w[2].seq_id == w[1].seq_id + 1 {
                let mut key = [0u8; 9];
                key[0..3].copy_from_slice(&w[0].comp_id);
                key[3..6].copy_from_slice(&w[1].comp_id);
                key[6..9].copy_from_slice(&w[2].comp_id);
                ref_3mers.entry(key).or_default().push(w[0].seq_id);
            }
        }
        for w in model.windows(3) {
            if w[1].seq_id == w[0].seq_id + 1 && w[2].seq_id == w[1].seq_id + 1 {
                let mut key = [0u8; 9];
                key[0..3].copy_from_slice(&w[0].comp_id);
                key[3..6].copy_from_slice(&w[1].comp_id);
                key[6..9].copy_from_slice(&w[2].comp_id);
                if let Some(r_starts) = ref_3mers.get(&key) {
                    for &rs in r_starts {
                        candidate_deltas.insert(w[0].seq_id - rs);
                    }
                }
            }
        }
    }

    let (mut best_delta, mut best_exact, mut best_overlap) = (0i32, 0usize, 0usize);
    for &delta in &candidate_deltas {
        let (mut exact, mut overlap) = (0usize, 0usize);
        for r in reference {
            if let Some(&mi) = model_by_seq.get(&(r.seq_id + delta)) {
                overlap += 1;
                if model[mi].comp_id == r.comp_id {
                    exact += 1;
                }
            }
        }
        if exact > best_exact || (exact == best_exact && overlap > best_overlap) {
            best_exact = exact;
            best_overlap = overlap;
            best_delta = delta;
        }
    }

    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(best_overlap);
    for (ri, r) in reference.iter().enumerate() {
        if let Some(&mi) = model_by_seq.get(&(r.seq_id + best_delta)) {
            pairs.push((mi, ri));
        }
    }
    let m = pairs.len();
    if m < 10 {
        return None;
    }

    let (mut preserved, mut n_local_pairs, mut sq_diff_sum, mut n_all_pairs) =
        (0.0f64, 0u64, 0.0f64, 0u64);
    let thresholds = [0.5f64, 1.0, 2.0, 4.0];

    for i in 0..m {
        let (mi, ri) = pairs[i];
        let m_ca_i = model[mi].ca;
        let r_ca_i = reference[ri].ca;
        for j in (i + 1)..m {
            let (mj, rj) = pairs[j];
            let d_ref = vec3_dist(r_ca_i, reference[rj].ca);
            let d_mod = vec3_dist(m_ca_i, model[mj].ca);
            let err = (d_mod - d_ref).abs();
            sq_diff_sum += err * err;
            n_all_pairs += 1;
            if d_ref < inclusion_radius {
                let mut hits = 0.0f64;
                for &th in &thresholds {
                    if err < th {
                        hits += 0.25;
                    }
                }
                preserved += hits;
                n_local_pairs += 1;
            }
        }
    }
    if n_local_pairs == 0 || n_all_pairs == 0 {
        return None;
    }

    let lddt_ca = preserved / (n_local_pairs as f64);
    let drmsd_ca = (sq_diff_sum / (n_all_pairs as f64)).sqrt();
    let matched_mean_plddt =
        pairs.iter().map(|&(mi, _)| model[mi].plddt).sum::<f64>() / (m as f64);
    let signed_bias = (matched_mean_plddt / 100.0) - lddt_ca;
    let seq_identity = best_exact as f64 / (m as f64);

    Some(ProteinModelVsExperiment {
        matched_residues: m,
        seq_identity,
        lddt_ca,
        drmsd_ca,
        matched_mean_plddt,
        signed_calibration_bias: signed_bias,
        abs_calibration_error: signed_bias.abs(),
    })
}

fn nerf_place(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    bond_len: f64,
    bond_angle_deg: f64,
    dihedral_deg: f64,
) -> [f64; 3] {
    let bc_raw = vec3_sub(c, b);
    let bc_len = vec3_norm(bc_raw);
    let bc = [bc_raw[0] / bc_len, bc_raw[1] / bc_len, bc_raw[2] / bc_len];
    let ab = vec3_sub(b, a);
    let n_raw = vec3_cross(ab, bc);
    let n_len = vec3_norm(n_raw);
    let n = [n_raw[0] / n_len, n_raw[1] / n_len, n_raw[2] / n_len];
    let m = vec3_cross(n, bc);
    let (sth, cth) = bond_angle_deg.to_radians().sin_cos();
    let (schi, cchi) = dihedral_deg.to_radians().sin_cos();
    [
        c[0] + bond_len * (-cth * bc[0] + sth * (cchi * m[0] + schi * n[0])),
        c[1] + bond_len * (-cth * bc[1] + sth * (cchi * m[1] + schi * n[1])),
        c[2] + bond_len * (-cth * bc[2] + sth * (cchi * m[2] + schi * n[2])),
    ]
}

/// Generate a physically exact Pauling–Corey ideal right-handed alpha-helix via NeRF kinematics.
pub fn build_ideal_alpha_helix_backbone(n_res: usize, plddt: f64) -> Vec<ProteinBackboneResidue> {
    let mut out = Vec::with_capacity(n_res);
    if n_res == 0 {
        return out;
    }
    let (l_nca, l_cac, l_cn) = (1.458f64, 1.525f64, 1.329f64);
    let (ang_ncac, ang_cacn, ang_cnca) = (111.2f64, 116.2f64, 121.7f64);
    let mut n_pos = [0.0f64, 0.0, 0.0];
    let mut ca_pos = [l_nca, 0.0, 0.0];
    let th0 = (180.0 - ang_ncac).to_radians();
    let mut c_pos = [l_nca + l_cac * th0.cos(), l_cac * th0.sin(), 0.0];
    out.push(ProteinBackboneResidue {
        seq_id: 1,
        comp_id: *b"ALA",
        n: Some(n_pos),
        ca: ca_pos,
        c: Some(c_pos),
        plddt,
    });
    for i in 1..n_res {
        let next_n = nerf_place(n_pos, ca_pos, c_pos, l_cn, ang_cacn, -47.0);
        let next_ca = nerf_place(ca_pos, c_pos, next_n, l_nca, ang_cnca, 180.0);
        let next_c = nerf_place(c_pos, next_n, next_ca, l_cac, ang_ncac, -57.0);
        n_pos = next_n;
        ca_pos = next_ca;
        c_pos = next_c;
        out.push(ProteinBackboneResidue {
            seq_id: (i as i32) + 1,
            comp_id: *b"ALA",
            n: Some(n_pos),
            ca: ca_pos,
            c: Some(c_pos),
            plddt,
        });
    }
    out
}

/// Exact Bondi (1964) van der Waals radius (Angstrom).
pub fn bondi_vdw_radius(element: &str) -> f64 {
    match element.to_ascii_uppercase().as_str() {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "P" => 1.80,
        "S" => 1.80,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        _ => 1.70,
    }
}

/// Grant & Pickup (1995) Gaussian atomic volume overlap `V(A, B)` with `p = 2.70`.
pub fn gaussian_shape_overlap(a: &[MolAtom3D], b: &[MolAtom3D]) -> f64 {
    const P_WT: f64 = 2.70;
    let pi = core::f64::consts::PI;
    let alpha_of = |r: f64| -> f64 {
        let r_safe = r.max(0.5);
        pi * (3.0 * P_WT / (4.0 * pi * r_safe * r_safe * r_safe)).powf(2.0 / 3.0)
    };
    let mut vol = 0.0f64;
    for ai in a {
        let a_i = alpha_of(ai.vdw_radius);
        for bj in b {
            let a_j = alpha_of(bj.vdw_radius);
            let denom = a_i + a_j;
            let dx = ai.pos[0] - bj.pos[0];
            let dy = ai.pos[1] - bj.pos[1];
            let dz = ai.pos[2] - bj.pos[2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let pre = P_WT * P_WT * (pi / denom).powf(1.5);
            vol += pre * (-(a_i * a_j / denom) * r2).exp();
        }
    }
    vol
}

/// Bounded 3D Gaussian volumetric Shape Tanimoto in `[0, 1]`.
pub fn gaussian_shape_tanimoto(a: &[MolAtom3D], b: &[MolAtom3D]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let v_aa = gaussian_shape_overlap(a, a);
    let v_bb = gaussian_shape_overlap(b, b);
    let v_ab = gaussian_shape_overlap(a, b);
    let denom = v_aa + v_bb - v_ab;
    if denom <= 1.0e-15 { 0.0 } else { (v_ab / denom).clamp(0.0, 1.0) }
}

/// Regularized Coulombic electrostatic potential (`e / A`) at each atom site.
pub fn molecular_coulomb_esp_at_atoms(atoms: &[MolAtom3D]) -> Vec<f64> {
    let mut esp = vec![0.0f64; atoms.len()];
    for i in 0..atoms.len() {
        let mut phi = atoms[i].partial_charge / 0.80;
        for j in 0..atoms.len() {
            if i == j {
                continue;
            }
            let d = vec3_dist(atoms[i].pos, atoms[j].pos).max(0.80);
            phi += atoms[j].partial_charge / d;
        }
        esp[i] = phi;
    }
    esp
}

/// Returns `(esp_tanimoto in [0,1], hodgkin_index in [-1,+1])`.
pub fn electrostatic_esp_similarity(
    a: &[MolAtom3D],
    b: &[MolAtom3D],
    alpha_spatial: f64,
    lambda_esp: f64,
) -> (f64, f64) {
    if a.is_empty() || b.is_empty() {
        return (0.0, 0.0);
    }
    let esp_a = molecular_coulomb_esp_at_atoms(a);
    let esp_b = molecular_coulomb_esp_at_atoms(b);

    let gauss_esp = |m1: &[MolAtom3D], e1: &[f64], m2: &[MolAtom3D], e2: &[f64]| -> f64 {
        let mut s = 0.0f64;
        for i in 0..m1.len() {
            for j in 0..m2.len() {
                let d = vec3_sub(m1[i].pos, m2[j].pos);
                let r2 = vec3_dot(d, d);
                let de = e1[i] - e2[j];
                s += (-alpha_spatial * r2 - lambda_esp * de * de).exp();
            }
        }
        s
    };
    let w_aa = gauss_esp(a, &esp_a, a, &esp_a);
    let w_bb = gauss_esp(b, &esp_b, b, &esp_b);
    let w_ab = gauss_esp(a, &esp_a, b, &esp_b);
    let denom_t = w_aa + w_bb - w_ab;
    let esp_tanimoto = if denom_t <= 1.0e-15 { 0.0 } else { (w_ab / denom_t).clamp(0.0, 1.0) };

    let signed_esp = |m1: &[MolAtom3D], e1: &[f64], m2: &[MolAtom3D], e2: &[f64]| -> f64 {
        let mut s = 0.0f64;
        for i in 0..m1.len() {
            for j in 0..m2.len() {
                let d = vec3_sub(m1[i].pos, m2[j].pos);
                let r2 = vec3_dot(d, d);
                s += e1[i] * e2[j] * (-alpha_spatial * r2).exp();
            }
        }
        s
    };
    let h_aa = signed_esp(a, &esp_a, a, &esp_a);
    let h_bb = signed_esp(b, &esp_b, b, &esp_b);
    let h_ab = signed_esp(a, &esp_a, b, &esp_b);
    let denom_h = h_aa + h_bb;
    let hodgkin = if denom_h.abs() <= 1.0e-15 { 0.0 } else { (2.0 * h_ab / denom_h).clamp(-1.0, 1.0) };
    (esp_tanimoto, hodgkin)
}

/// Directional 3D pharmacophore Tanimoto overlap in `[0, 1]`.
pub fn directional_pharmacophore_tanimoto(
    a: &[Pharmacophore3D],
    b: &[Pharmacophore3D],
    alpha_p: f64,
) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let overlap = |p1: &[Pharmacophore3D], p2: &[Pharmacophore3D]| -> f64 {
        let mut s = 0.0f64;
        for f1 in p1 {
            for f2 in p2 {
                if f1.kind != f2.kind {
                    continue;
                }
                let dir_wt = match (f1.dir, f2.dir) {
                    (Some(u), Some(v)) => {
                        let nu = vec3_norm(u);
                        let nv = vec3_norm(v);
                        if nu < 1.0e-12 || nv < 1.0e-12 {
                            0.0
                        } else {
                            let c = vec3_dot(u, v) / (nu * nv);
                            if f1.kind == 3 { c.abs() } else { c.max(0.0) }
                        }
                    }
                    (None, None) => 1.0,
                    _ => 0.0,
                };
                if dir_wt <= 0.0 {
                    continue;
                }
                let d = vec3_sub(f1.pos, f2.pos);
                s += dir_wt * (-0.5 * alpha_p * vec3_dot(d, d)).exp();
            }
        }
        s
    };
    let p_aa = overlap(a, a);
    let p_bb = overlap(b, b);
    let p_ab = overlap(a, b);
    let denom = p_aa + p_bb - p_ab;
    if denom <= 1.0e-15 { 0.0 } else { (p_ab / denom).clamp(0.0, 1.0) }
}

/// $SO(3)$-invariant signed tetrahedral volume (flips sign under improper reflection / chirality error).
pub fn signed_tetrahedral_volume(
    a0: [f64; 3],
    a1: [f64; 3],
    a2: [f64; 3],
    a3: [f64; 3],
) -> f64 {
    let e1 = vec3_sub(a1, a0);
    let e2 = vec3_sub(a2, a0);
    let e3 = vec3_sub(a3, a0);
    vec3_dot(vec3_cross(e1, e2), e3) / 6.0
}

/// Pure-Rust mmCIF `HETATM` heavy-atom ligand parser.
pub fn parse_mmcif_hetatm_ligand(
    mmcif_text: &str,
    target_comp_id: &str,
    preferred_chain: Option<&str>,
) -> Vec<MolAtom3D> {
    let mut in_loop = false;
    let mut reading_atom_headers = false;
    let mut in_atom_rows = false;
    let mut headers: Vec<&str> = Vec::new();
    let (mut idx_group, mut idx_comp, mut idx_label_asym, mut idx_auth_asym) = (None, None, None, None);
    let (mut idx_alt, mut idx_elem, mut idx_x, mut idx_y, mut idx_z, mut idx_model) =
        (None, None, None, None, None, None);

    let want_comp = target_comp_id.to_ascii_uppercase();
    let mut locked_chain: Option<String> = preferred_chain.map(str::to_string);
    let mut out = Vec::new();

    for raw_line in mmcif_text.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        if line == "loop_" {
            if in_atom_rows {
                break;
            }
            in_loop = true;
            reading_atom_headers = false;
            headers.clear();
            continue;
        }
        if in_loop && line.starts_with('_') {
            if line.starts_with("_atom_site.") {
                reading_atom_headers = true;
                headers.push(line);
            } else if reading_atom_headers {
                in_loop = false;
                reading_atom_headers = false;
                headers.clear();
            } else {
                in_loop = false;
            }
            continue;
        }
        if reading_atom_headers && !line.starts_with('_') {
            reading_atom_headers = false;
            in_loop = false;
            in_atom_rows = true;
            for (i, h) in headers.iter().enumerate() {
                match *h {
                    "_atom_site.group_PDB" => idx_group = Some(i),
                    "_atom_site.label_comp_id" => idx_comp = Some(i),
                    "_atom_site.label_asym_id" => idx_label_asym = Some(i),
                    "_atom_site.auth_asym_id" => idx_auth_asym = Some(i),
                    "_atom_site.label_alt_id" => idx_alt = Some(i),
                    "_atom_site.type_symbol" => idx_elem = Some(i),
                    "_atom_site.Cartn_x" => idx_x = Some(i),
                    "_atom_site.Cartn_y" => idx_y = Some(i),
                    "_atom_site.Cartn_z" => idx_z = Some(i),
                    "_atom_site.pdbx_PDB_model_num" => idx_model = Some(i),
                    _ => {}
                }
            }
        }
        if in_atom_rows {
            if line.starts_with('#') || line.starts_with("data_") || line.starts_with('_') {
                break;
            }
            let toks = split_mmcif_tokens(line);
            if toks.len() < headers.len() {
                continue;
            }
            if let Some(ig) = idx_group {
                if toks[ig] != "HETATM" {
                    continue;
                }
            }
            if let Some(im) = idx_model {
                if toks[im] != "1" {
                    continue;
                }
            }
            if let Some(ia) = idx_alt {
                let alt = toks[ia];
                if alt != "." && alt != "?" && alt != "A" {
                    continue;
                }
            }
            let comp = idx_comp.map(|i| toks[i]).unwrap_or("");
            if !comp.eq_ignore_ascii_case(&want_comp) {
                continue;
            }
            let elem = idx_elem.map(|i| toks[i]).unwrap_or("C");
            if elem.eq_ignore_ascii_case("H") || elem.eq_ignore_ascii_case("D") {
                continue;
            }
            let auth_ch = idx_auth_asym.map(|i| toks[i]).unwrap_or("");
            let label_ch = idx_label_asym.map(|i| toks[i]).unwrap_or("");
            if let Some(ref want) = locked_chain {
                if auth_ch != want.as_str() && label_ch != want.as_str() {
                    continue;
                }
            } else {
                let first = if !auth_ch.is_empty() && auth_ch != "." && auth_ch != "?" {
                    auth_ch
                } else {
                    label_ch
                };
                locked_chain = Some(first.to_string());
            }
            let (Some(ix), Some(iy), Some(iz)) = (idx_x, idx_y, idx_z) else {
                continue;
            };
            let (Ok(x), Ok(y), Ok(z)) = (
                toks[ix].parse::<f64>(),
                toks[iy].parse::<f64>(),
                toks[iz].parse::<f64>(),
            ) else {
                continue;
            };
            let approx_q = match elem.to_ascii_uppercase().as_str() {
                "O" => -0.45,
                "N" => -0.30,
                "F" | "CL" => -0.15,
                "S" | "P" => 0.15,
                _ => 0.08,
            };
            out.push(MolAtom3D {
                pos: [x, y, z],
                vdw_radius: bondi_vdw_radius(elem),
                partial_charge: approx_q,
            });
        }
    }
    out
}

/// Count severe steric clashes between heavy ligand atoms and protein backbone (`N`, `CA`, `C`).
pub fn audit_protein_ligand_clashes(
    backbone: &[ProteinBackboneResidue],
    ligand: &[MolAtom3D],
    clash_threshold_angstrom: f64,
) -> (usize, f64) {
    let mut clashes = 0usize;
    let mut min_d = f64::INFINITY;
    for lat in ligand {
        for res in backbone {
            for opt_pt in [res.n, Some(res.ca), res.c] {
                if let Some(pt) = opt_pt {
                    let d = vec3_dist(lat.pos, pt);
                    if d < min_d {
                        min_d = d;
                    }
                    if d < clash_threshold_angstrom {
                        clashes += 1;
                    }
                }
            }
        }
    }
    (clashes, min_d)
}

/// Offline Wave-56 protein severe-testing battery (`5/5` Popperian mutants killed, `S_Mayo = 1.0`).
pub fn run_protein_severe_testing_battery(
    allow_static_clinical_promotion: bool,
    allow_overconfident_misfold: bool,
) -> StructuralSevereTestReceipt {
    let helix = build_ideal_alpha_helix_backbone(24, 96.0);
    let pos_audit = audit_protein_backbone(&helix);
    let pos_cmp = compare_model_vs_pdb_lddt_ca(&helix, &helix, 15.0).unwrap();
    let pos_ok = (pos_cmp.lddt_ca - 1.0).abs() < 1.0e-12
        && pos_audit.steric_clashes_ca == 0
        && pos_audit.ca_bond_breaks == 0
        && pos_audit.ramachandran_outliers == 0;

    let mut m_clash = helix.clone();
    m_clash[10].ca = [helix[2].ca[0] + 1.4, helix[2].ca[1], helix[2].ca[2]];
    let clash_caught = audit_protein_backbone(&m_clash).steric_clashes_ca >= 1;

    let mut m_break = helix.clone();
    for r in &mut m_break[12..] {
        r.ca[2] += 7.5;
    }
    let break_caught = audit_protein_backbone(&m_break).ca_bond_breaks >= 1;

    let alt_fold: Vec<ProteinBackboneResidue> = (0..24)
        .map(|i| ProteinBackboneResidue {
            seq_id: (i as i32) + 1,
            comp_id: *b"ALA",
            n: None,
            ca: [(i as f64) * 3.45, if i % 2 == 0 { 0.8 } else { -0.8 }, 0.0],
            c: None,
            plddt: 0.0,
        })
        .collect();
    let fs_cmp = compare_model_vs_pdb_lddt_ca(&helix, &alt_fold, 15.0).unwrap();
    let foldswitch_caught = fs_cmp.abs_calibration_error > 0.35 && !allow_overconfident_misfold;

    let idr_helix = build_ideal_alpha_helix_backbone(24, 45.0);
    let idr_audit = audit_protein_backbone(&idr_helix);
    let idr_caught = idr_audit.frac_disordered > 0.5 && idr_audit.helical_fraction_in_low_plddt > 0.5;

    let locked = !allow_static_clinical_promotion && !allow_overconfident_misfold;
    let mutants_killed = [clash_caught, break_caught, foldswitch_caught, idr_caught, locked]
        .iter()
        .filter(|&&b| b)
        .count();
    let s_mayo = mutants_killed as f64 / 5.0;
    StructuralSevereTestReceipt {
        positive_control_ok: pos_ok,
        mutants_killed,
        total_mutants: 5,
        mayo_severity_power: s_mayo,
        epistemic_demarcation_locked: locked,
        passed: pos_ok && mutants_killed == 5 && locked,
    }
}

/// Offline Wave-57 3D molecular field & chirality severe-testing battery (`5/5` mutants killed).
pub fn run_ligand_severe_testing_battery(
    allow_unverified_affinity_claim: bool,
    allow_chirality_inversion: bool,
) -> StructuralSevereTestReceipt {
    let ref_atoms = vec![
        MolAtom3D { pos: [0.0, 0.0, 0.0], vdw_radius: 1.70, partial_charge: 0.25 },
        MolAtom3D { pos: [1.46, 0.0, 0.0], vdw_radius: 1.55, partial_charge: -0.65 },
        MolAtom3D { pos: [-0.50, 1.44, 0.0], vdw_radius: 1.52, partial_charge: -0.70 },
        MolAtom3D { pos: [-0.50, -0.72, 1.25], vdw_radius: 1.70, partial_charge: 0.55 },
        MolAtom3D { pos: [-0.50, -0.72, -1.25], vdw_radius: 1.47, partial_charge: 0.55 },
    ];
    let ref_pharms = vec![
        Pharmacophore3D { kind: 1, pos: [1.46, 0.0, 0.0], dir: Some([1.0, 0.0, 0.0]) },
        Pharmacophore3D { kind: 2, pos: [-0.50, 1.44, 0.0], dir: Some([0.0, 1.0, 0.0]) },
    ];
    let receptor = vec![ProteinBackboneResidue {
        seq_id: 1,
        comp_id: *b"ALA",
        n: Some([6.5, 0.0, 0.0]),
        ca: [7.8, 0.0, 0.0],
        c: Some([8.5, 1.2, 0.0]),
        plddt: 95.0,
    }];

    let pos_shape = gaussian_shape_tanimoto(&ref_atoms, &ref_atoms);
    let (pos_esp, pos_hodgkin) = electrostatic_esp_similarity(&ref_atoms, &ref_atoms, 0.81, 8.0);
    let pos_pharm = directional_pharmacophore_tanimoto(&ref_pharms, &ref_pharms, 1.0);
    let (pos_clashes, _) = audit_protein_ligand_clashes(&receptor, &ref_atoms, 2.0);
    let ref_vol = signed_tetrahedral_volume(
        ref_atoms[1].pos,
        ref_atoms[2].pos,
        ref_atoms[3].pos,
        ref_atoms[4].pos,
    );
    let pos_ok = (pos_shape - 1.0).abs() < 1.0e-12
        && (pos_esp - 1.0).abs() < 1.0e-12
        && (pos_hodgkin - 1.0).abs() < 1.0e-12
        && (pos_pharm - 1.0).abs() < 1.0e-12
        && pos_clashes == 0
        && ref_vol.abs() > 0.5;

    let mut inv_q = ref_atoms.clone();
    for a in &mut inv_q {
        a.partial_charge = -a.partial_charge;
    }
    let (inv_esp, inv_hodgkin) = electrostatic_esp_similarity(&ref_atoms, &inv_q, 0.81, 8.0);
    let charge_caught = (gaussian_shape_tanimoto(&ref_atoms, &inv_q) - 1.0).abs() < 1.0e-12
        && inv_hodgkin < -0.99
        && inv_esp < 0.25;

    let mut mirror = ref_atoms.clone();
    for a in &mut mirror {
        a.pos[2] = -a.pos[2];
    }
    let mirror_vol = signed_tetrahedral_volume(
        mirror[1].pos,
        mirror[2].pos,
        mirror[3].pos,
        mirror[4].pos,
    );
    let chiral_caught = (ref_vol * mirror_vol < 0.0) && !allow_chirality_inversion;

    let rev_pharms: Vec<Pharmacophore3D> = ref_pharms
        .iter()
        .map(|p| Pharmacophore3D {
            kind: p.kind,
            pos: p.pos,
            dir: p.dir.map(|v| [-v[0], -v[1], -v[2]]),
        })
        .collect();
    let pharm_caught = directional_pharmacophore_tanimoto(&ref_pharms, &rev_pharms, 1.0).abs() < 1.0e-12;

    let mut clash_lig = ref_atoms.clone();
    clash_lig[0].pos = [7.8, 1.2, 0.0];
    let (mut_clashes, _) = audit_protein_ligand_clashes(&receptor, &clash_lig, 2.0);
    let clash_caught = mut_clashes >= 1;

    let locked = !allow_unverified_affinity_claim && !allow_chirality_inversion;
    let mutants_killed = [charge_caught, chiral_caught, pharm_caught, clash_caught, locked]
        .iter()
        .filter(|&&b| b)
        .count();
    let s_mayo = mutants_killed as f64 / 5.0;
    StructuralSevereTestReceipt {
        positive_control_ok: pos_ok,
        mutants_killed,
        total_mutants: 5,
        mayo_severity_power: s_mayo,
        epistemic_demarcation_locked: locked,
        passed: pos_ok && mutants_killed == 5 && locked,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iupac_dihedral_sign_conventions() {
        assert!((protein_dihedral_deg([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 1.0]) - 90.0).abs() < 1.0e-10);
        assert!((protein_dihedral_deg([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 1.0]) + 90.0).abs() < 1.0e-10);
        assert!(protein_dihedral_deg([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0]).abs() < 1.0e-10);
        assert!((protein_dihedral_deg([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 1.0]).abs() - 180.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_mmcif_backbone_and_hetatm_ligand_parsers() {
        let mini_cif = "\
data_MINI
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.type_symbol
_atom_site.pdbx_PDB_model_num
ATOM   N  . ALA A 1 0.0 0.0 0.0 95.0 10 A N 1
ATOM   CA . ALA A 1 1.4 0.0 0.0 95.5 10 A C 1
ATOM   C  . ALA A 1 2.0 1.4 0.0 95.0 10 A C 1
ATOM   N  . GLY A 2 3.3 1.4 0.0 92.0 11 A N 1
ATOM   CA . GLY A 2 4.9 1.4 0.0 92.2 11 A C 1
ATOM   C  . GLY A 2 5.5 2.8 0.0 92.0 11 A C 1
HETATM C1 . LIG A . 9.0 9.0 9.0 20.0 99 A C 1
HETATM O1 . LIG A . 10.2 9.0 9.0 20.0 99 A O 1
HETATM H1 . LIG A . 10.8 9.0 9.0 20.0 99 A H 1
#
";
        let prot = parse_mmcif_backbone_chain(mini_cif, Some("A"));
        assert_eq!(prot.len(), 2);
        assert_eq!(&prot[0].comp_id, b"ALA");
        assert!((prot[0].plddt - 95.5).abs() < 1.0e-12);
        let lig = parse_mmcif_hetatm_ligand(mini_cif, "LIG", Some("A"));
        assert_eq!(lig.len(), 2, "Hydrogen H1 must be skipped, leaving 2 heavy atoms");
        let (clashes, min_d) = audit_protein_ligand_clashes(&prot, &lig, 2.0);
        assert_eq!(clashes, 0);
        assert!(min_d > 5.0);
    }

    #[test]
    fn test_protein_and_ligand_severe_testing_batteries_and_fail_closed_controls() {
        let p_ok = run_protein_severe_testing_battery(false, false);
        assert!(p_ok.passed);
        assert_eq!(p_ok.mutants_killed, 5);
        assert!((p_ok.mayo_severity_power - 1.0).abs() < 1.0e-12);
        assert!(!run_protein_severe_testing_battery(true, false).passed);
        assert!(!run_protein_severe_testing_battery(false, true).passed);

        let l_ok = run_ligand_severe_testing_battery(false, false);
        assert!(l_ok.passed);
        assert_eq!(l_ok.mutants_killed, 5);
        assert!((l_ok.mayo_severity_power - 1.0).abs() < 1.0e-12);
        assert!(!run_ligand_severe_testing_battery(true, false).passed);
        assert!(!run_ligand_severe_testing_battery(false, true).passed);
    }

    #[test]
    fn test_so3_invariance_of_tetrahedral_chiral_volume() {
        let (st, ct) = 1.2345f64.sin_cos();
        let rot = |p: [f64; 3]| -> [f64; 3] {
            [ct * p[0] - st * p[1] + 4.2, st * p[0] + ct * p[1] - 1.7, p[2] + 3.1]
        };
        let (p0, p1, p2, p3) = (
            [1.46, 0.0, 0.0],
            [-0.50, 1.44, 0.0],
            [-0.50, -0.72, 1.25],
            [-0.50, -0.72, -1.25],
        );
        let v0 = signed_tetrahedral_volume(p0, p1, p2, p3);
        let v1 = signed_tetrahedral_volume(rot(p0), rot(p1), rot(p2), rot(p3));
        assert!((v0 - v1).abs() < 1.0e-12);
    }
}
