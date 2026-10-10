#!/usr/bin/env python3
"""
rank_binders.py
===============
Stage 3 of the AETh3R de novo MDM2 binder selection funnel.

Structure-based ranking of de novo MDM2 peptide binders for the BEST-BINDER
objective (non-canonical solutions allowed; pharmacophore fidelity NOT required).

WHAT THIS DID, AND WHY IT EXISTS
--------------------------------
Stages 1 and 2 (z-score_filter.py) rank on AlphaFold confidence metrics and
flag the p53 pharmacophore. Both work on sequence and predicted-geometry
scores, and neither can see whether a peptide is chemically plausible for a
hydrophobic groove. Stage 2's pick, Z2 (267.3), carries net charge -8 and
still satisfies both pharmacophore flags.

This stage reads the atomic coordinates instead, and is the first point at
which that is caught: it computes net charge from the structure and rejects
Z2 at the chemistry gate (G2 below). It is the reason the funnel has three
stages rather than two.

INPUT
-----
    --pdbdir     the 2,000 AlphaFold Multimer PDBs, chain A = MDM2 receptor,
                 chain B = peptide. NOT COMMITTED — too large for this
                 repository. Filenames must match design<NN>_n<M>.pdb for
                 the design/n identifiers to be recovered.
    --pockets    pockets_chimerax.json (default), beside this script
    --benchmark  bench_PMI_norm.pdb, beside this script

PRODUCED
--------
    structural_ranking_full.csv  all 2,000 rows, every metric and gate flag
    selected_binders.csv         the meets-or-exceeds-PMI shortlist

    Both are committed beside this script under the names the team archived
    them as: gate_filter_results_full.csv and gate_filter_top_binders.csv.
    Column order matches this script's output exactly.

    From the committed table:
        2,000 candidates
        1,275 fill all three sub-pockets
          409 pass all three gates
           18 meet or exceed PMI on both axes  -> selected_binders.csv
    Gate failures (not disjoint): 1,181 confidence, 481 chemistry, 634 site.

    The winner is Z3 (131.2), MKPEEVGEMLLSVLK: 49 contacts, 3/3 pockets,
    net charge -1, interface pLDDT 89.46, apolar BSA 42.02 Å2/residue.
    Z2 (267.3) appears in the full table with passed_gates = False.

NOT RE-RUNNABLE AS COMMITTED
----------------------------
The 2,000 input PDBs are not in this repository, so the committed tables
could not be regenerated during the 2026 documentation pass and were NOT
verified by re-running. They are taken as given. The figures downstream of
them were re-generated and do reproduce.

Ranking is an INTERSECTION of independent axes, not a weighted sum, so no
axis's blindness can dominate and the selection is a defensible logical
statement rather than an arbitrary weighting.

Pipeline
--------
GATES (must pass all):
  G1  Confidence : mean interface pLDDT >= --plddt-floor            (from B-factor col)
  G2  Chemistry  : net charge within reference envelope             (charge >= --charge-floor)
  G3  Site       : peptide bound in the p53 cleft, not a peripheral
                   patch (interface centroid within --site-tol Å of
                   the native pocket centroid)

RANK (on survivors, three structural axes):
  A1  Apolar buried surface area  (hydrophobic burial; the primary affinity proxy)
  A2  Sub-pocket occupancy        (all three of F19/W23/L26 pockets contacted)
  A3  Intermolecular contact count (interface extent / packing)

INTERSECTION SELECTION ("best binder"):
  a candidate is SELECTED iff it passes all gates AND
    - occupies all THREE sub-pockets, AND
    - is top --top-pct on apolar BSA, AND
    - is top --top-pct on contact count.

Methods / citations
-------------------
- SASA via Lee-Richards / Shrake-Rupley, computed with FreeSASA
  (Mitternacht, F1000Research 2016). BSA = (SASA_A + SASA_B - SASA_AB)/2
  (per-side interface area).
- Buried surface area as a binding proxy: Chen, Kurgan & Ruan; and the
  contacts-based affinity rationale of Vangone & Bonvin (eLife 2015),
  which also underlies PRODIGY (Xue et al., Bioinformatics 2016).
- Sub-pockets defined from the native p53-MDM2 complex 1YCR
  (Kussie et al., Science 1996). The committed pockets_chimerax.json was
  built from ChimeraX interatomic contacts (<=4.5 A), shared wall residues
  allowed in more than one pocket. An earlier geometric nearest-anchor
  variant (pockets.json, from derive_pockets.py) is referenced in places
  but was never committed to this repository and is not available.
- Apolar vs polar atom classification follows the standard
  carbon/sulfur = apolar, N/O = polar convention used in SASA-based
  interface analysis.

Usage
-----
    python3 rank_binders.py --pdbdir /path/to/pdbs --pockets pockets.json \
        --benchmark bench_PMI_norm.pdb --outdir results_struct

    (The ranker reads peptide sequence, charge and metrics from the PDBs
     themselves, so no FASTA is needed.)

Assumes each PDB has chain A = MDM2 receptor, chain B = 15-mer peptide,
matching the provided examples.
"""
import argparse
import glob
import json
import os
import re
import sys
import numpy as np
import pandas as pd

import freesasa
from Bio.PDB import PDBParser

freesasa.setVerbosity(freesasa.silent)
PARSER = PDBParser(QUIET=True)

# Carbon and sulfur count as apolar, N and O as polar. Standard convention in
# SASA-based interface analysis; hydrogens are excluded everywhere (see heavy()).
APOLAR_ELEMENTS = {"C", "S"}

# Å, residue-residue heavy-atom contact. 5.0 Å is the PRODIGY convention
# (Vangone & Bonvin, eLife 2015; Xue et al., Bioinformatics 2016), whose
# contact-count-vs-affinity relationship is the reason contacts are a ranking
# axis here at all.
CONTACT_CUTOFF = 5.0

# Å, peptide heavy atom -> sub-pocket lining residue. Matches the 4.5 Å cutoff
# the pocket definitions themselves were built at — see "cutoff" in
# pockets_chimerax.json — so occupancy is tested on the same criterion that
# decided which residues line a pocket.
POCKET_CUTOFF = 4.5


def heavy(residue):
    return [a for a in residue if a.element != "H"]


def sasa_areas(pdb_path):
    """Return (SASA_complex, SASA_receptorA_alone, SASA_peptideB_alone).

    DEAD CODE: nothing calls this. It is the whole-chain BSA calculation,
    superseded by apolar_bsa() below, which does the same subtraction
    per atom so the apolar and polar contributions can be separated.
    Kept rather than deleted so the simpler form stays on the record;
    flagged for the team to remove if they prefer.
    """
    # complex
    st = freesasa.Structure(pdb_path)
    complex_area = freesasa.calc(st).totalArea()

    # isolate chains by writing temp single-chain structures via Bio.PDB selection
    model = PARSER.get_structure("s", pdb_path)[0]

    def chain_sasa(chain_id):
        # build a freesasa Structure from atoms of one chain
        s = freesasa.Structure()
        idx = 0
        for res in model[chain_id]:
            if res.id[0] != " ":
                continue
            for atom in heavy(res):
                x, y, z = atom.coord
                s.addAtom(atom.get_name().ljust(4)[:4],
                          res.resname.ljust(3)[:3],
                          str(res.id[1]).rjust(4)[:4],
                          chain_id,
                          float(x), float(y), float(z))
                idx += 1
        if idx == 0:
            return 0.0
        return freesasa.calc(s).totalArea()

    return complex_area, chain_sasa("A"), chain_sasa("B")


def apolar_bsa(pdb_path):
    """Apolar and total per-side BSA using atom-classified SASA subtraction."""
    model = PARSER.get_structure("s", pdb_path)[0]

    def atom_sasa(structure_atoms):
        s = freesasa.Structure()
        for (nm, rn, ri, ch, coord, elem) in structure_atoms:
            s.addAtom(nm.ljust(4)[:4], rn.ljust(3)[:3], str(ri).rjust(4)[:4], ch,
                      float(coord[0]), float(coord[1]), float(coord[2]))
        r = freesasa.calc(s)
        # per-atom areas
        areas = [r.atomArea(i) for i in range(s.nAtoms())]
        return np.array(areas)

    def collect(chain_ids):
        out = []
        for ch in chain_ids:
            for res in model[ch]:
                if res.id[0] != " ":
                    continue
                for a in heavy(res):
                    out.append((a.get_name(), res.resname, res.id[1], ch, a.coord, a.element))
        return out

    atoms_AB = collect(["A", "B"])
    atoms_A = collect(["A"])
    atoms_B = collect(["B"])

    sasa_AB = atom_sasa(atoms_AB)
    sasa_A = atom_sasa(atoms_A)
    sasa_B = atom_sasa(atoms_B)

    # buried per atom = free - complex, summed; split by apolar/polar element
    nA = len(atoms_A)
    elems_A = np.array([a[5] for a in atoms_A])
    elems_B = np.array([a[5] for a in atoms_B])
    apolar_A = np.array([e in APOLAR_ELEMENTS for e in elems_A])
    apolar_B = np.array([e in APOLAR_ELEMENTS for e in elems_B])

    buried_A = sasa_A - sasa_AB[:nA]
    buried_B = sasa_B - sasa_AB[nA:]

    total_bsa = (buried_A.sum() + buried_B.sum()) / 2.0
    apolar_bsa_val = (buried_A[apolar_A].sum() + buried_B[apolar_B].sum()) / 2.0
    return total_bsa, apolar_bsa_val


def interface_metrics(pdb_path, pockets_design):
    """Contacts, pocket occupancy, site check, interface pLDDT."""
    model = PARSER.get_structure("s", pdb_path)[0]
    rec = model["A"]
    pep = model["B"]

    rec_res = [r for r in rec if r.id[0] == " "]
    pep_res = [r for r in pep if r.id[0] == " "]

    rec_coords = {r.id[1]: np.array([a.coord for a in heavy(r)]) for r in rec_res}
    pep_coords = {r.id[1]: np.array([a.coord for a in heavy(r)]) for r in pep_res}

    # residue-residue contacts (<=5 A heavy atom)
    contacts = 0
    interface_rec_ids = set()
    interface_pep_atoms = []
    for prid, pc in pep_coords.items():
        for rrid, rc in rec_coords.items():
            d = np.linalg.norm(pc[:, None, :] - rc[None, :, :], axis=2).min()
            if d <= CONTACT_CUTOFF:
                contacts += 1
                interface_rec_ids.add(rrid)

    # pocket occupancy: for each sub-pocket, is any peptide heavy atom within
    # POCKET_CUTOFF of any lining-residue heavy atom?
    all_pep = np.vstack(list(pep_coords.values()))
    occ = {}
    for pname, lining in pockets_design.items():
        hit = False
        for rrid in lining:
            if rrid not in rec_coords:
                continue
            d = np.linalg.norm(all_pep[:, None, :] - rec_coords[rrid][None, :, :], axis=2).min()
            if d <= POCKET_CUTOFF:
                hit = True
                break
        occ[pname] = hit
    pockets_occupied = sum(occ.values())

    # site check: peptide centroid vs native pocket-lining centroid
    pocket_all_ids = [r for v in pockets_design.values() for r in v]
    pocket_centroid = np.mean(
        [rec_coords[i].mean(0) for i in pocket_all_ids if i in rec_coords], axis=0
    )
    pep_centroid = all_pep.mean(0)
    site_dist = float(np.linalg.norm(pep_centroid - pocket_centroid))

    # interface pLDDT: mean B-factor over peptide residues at the interface
    pep_iface_bfac = []
    for prid, pc in pep_coords.items():
        # is this peptide residue in contact with receptor?
        in_contact = False
        for rrid, rc in rec_coords.items():
            d = np.linalg.norm(pc[:, None, :] - rc[None, :, :], axis=2).min()
            if d <= CONTACT_CUTOFF:
                in_contact = True
                break
        if in_contact:
            res = pep[(" ", prid, " ")]
            pep_iface_bfac.extend([a.get_bfactor() for a in heavy(res)])
    iface_plddt = float(np.mean(pep_iface_bfac)) if pep_iface_bfac else float("nan")

    return {
        "contacts": contacts,
        "pockets_occupied": pockets_occupied,
        "occ_F19": occ["F19"], "occ_W23": occ["W23"], "occ_L26": occ["L26"],
        "site_dist": site_dist,
        "iface_plddt": iface_plddt,
    }


def peptide_seq_and_charge(pdb_path):
    three2one = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
        'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
        'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
    model = PARSER.get_structure("s", pdb_path)[0]
    seq = "".join(three2one.get(r.resname, "X") for r in model["B"] if r.id[0] == " ")
    # Crude net charge at pH ~7.4: K,R = +1; D,E = -1; H ~ +0.1 (histidine is
    # mostly neutral at 7.4, pKa ~6). Counting only, no pKa model — this is
    # the number the G2 chemistry gate tests, and it is what rejects Z2 at
    # -8. analyze_designs.py computes a different charge for the same
    # sequences using Biopython's charge_at_pH, so the two will not agree
    # exactly; this one is the one the selection acted on.
    charge = seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E") + 0.1*seq.count("H")
    acidic = (seq.count("D") + seq.count("E")) / len(seq) if seq else 0
    return seq, charge, acidic


def id_from_filename(path):
    m = re.search(r"design[_:]?(\d+)[_ ]?n[:_]?(\d+)", os.path.basename(path))
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbdir", required=True)
    ap.add_argument("--pockets", default="pockets_chimerax.json")
    ap.add_argument("--outdir", default="results_struct")
    ap.add_argument("--plddt-floor", type=float, default=70.0,
                    help="min mean interface pLDDT (B-factor units)")
    ap.add_argument("--charge-floor", type=float, default=-4.0,
                    help="reject peptides more negative than this net charge")
    # WARNING: this default is STALE. The committed tables were produced with
    # --site-tol 14, which is what README.md's invocation passes and what
    # make_distribution_figure.py calls the "corrected" value. Running with
    # the 8.0 default gives a different site-gate count and a different
    # selection. The default is left at 8.0 because changing it would change
    # this script's behaviour; always pass --site-tol 14 to reproduce the
    # committed results.
    ap.add_argument("--site-tol", type=float, default=8.0,
                    help="max Å from native pocket centroid (site gate)")
    ap.add_argument("--top-pct", type=float, default=10.0,
                    help="top percentile threshold for BSA and contacts (used only if no --benchmark)")
    ap.add_argument("--benchmark", default=None,
                    help="PDB of a reference complex (e.g. AF3 PMI) scored identically; "
                         "switches selection to meets-or-exceeds-benchmark. Must be A=receptor, B=peptide.")
    ap.add_argument("--limit", type=int, default=0, help="debug: process only N files")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pk = json.load(open(args.pockets))
    pockets_design = {k: v for k, v in pk["pockets_design"].items()}

    files = sorted(glob.glob(os.path.join(args.pdbdir, "*.pdb")))
    if args.limit:
        files = files[:args.limit]
    if not files:
        sys.exit(f"No PDBs found in {args.pdbdir}")

    rows = []
    for i, f in enumerate(files, 1):
        try:
            d, n = id_from_filename(f)
            seq, charge, acidic = peptide_seq_and_charge(f)
            total_bsa, ap_bsa = apolar_bsa(f)
            im = interface_metrics(f, pockets_design)
            rows.append({
                "design": d, "n": n, "file": os.path.basename(f), "peptide": seq,
                "net_charge": round(charge, 2), "acidic_frac": round(acidic, 3),
                "total_bsa": round(total_bsa, 1), "apolar_bsa": round(ap_bsa, 1),
                **im,
            })
        except Exception as e:
            print(f"  ! {os.path.basename(f)}: {e}", file=sys.stderr)
        if i % 200 == 0:
            print(f"  processed {i}/{len(files)}")

    df = pd.DataFrame(rows)

    # ---- gates ----
    df["gate_conf"] = df["iface_plddt"] >= args.plddt_floor
    df["gate_chem"] = df["net_charge"] >= args.charge_floor
    df["gate_site"] = df["site_dist"] <= args.site_tol
    df["passed_gates"] = df["gate_conf"] & df["gate_chem"] & df["gate_site"]

    # per-residue apolar BSA (length-normalised so 12-mers and 15-mers compare fairly)
    df["pep_len"] = df["peptide"].str.len()
    df["apolar_bsa_per_res"] = df["apolar_bsa"] / df["pep_len"]

    # ---- optional benchmark (e.g. AF3 PMI) ----
    bench = None
    if args.benchmark:
        b_total, b_apolar = apolar_bsa(args.benchmark)
        b_im = interface_metrics(args.benchmark, pockets_design)
        b_seq, b_charge, b_acidic = peptide_seq_and_charge(args.benchmark)
        bench = {
            "peptide": b_seq, "net_charge": b_charge,
            "apolar_bsa_val": b_apolar,
            "apolar_bsa_per_res": b_apolar / len(b_seq),
            "contacts": b_im["contacts"], "pockets_occupied": b_im["pockets_occupied"],
        }
        print(f"\nBenchmark ({os.path.basename(args.benchmark)}): {b_seq}")
        print(f"  apolar_bsa={b_apolar:.1f}  per_res={bench['apolar_bsa_per_res']:.1f}"
              f"  contacts={b_im['contacts']}  pockets={b_im['pockets_occupied']}/3")

    # ---- rank axes among gate-passers ----
    g = df[df["passed_gates"]].copy()
    if len(g):
        g["pct_apolar_bsa"] = g["apolar_bsa"].rank(pct=True) * 100
        g["pct_contacts"] = g["contacts"].rank(pct=True) * 100
        if bench is not None:
            # meets-or-exceeds-PMI selection (biologically anchored intersection)
            g["selected"] = (
                (g["pockets_occupied"] == 3) &
                (g["apolar_bsa_per_res"] >= bench["apolar_bsa_per_res"]) &
                (g["contacts"] >= bench["contacts"])
            )
        else:
            # fall back to top-percentile intersection
            thr = 100 - args.top_pct
            g["selected"] = (
                (g["pockets_occupied"] == 3) &
                (g["pct_apolar_bsa"] >= thr) &
                (g["pct_contacts"] >= thr)
            )
        g = g.sort_values(["selected", "apolar_bsa_per_res", "contacts"],
                          ascending=[False, False, False])
        df = df.merge(
            g[["file", "pct_apolar_bsa", "pct_contacts", "selected"]],
            on="file", how="left"
        )
    else:
        df["selected"] = False

    df.to_csv(os.path.join(args.outdir, "structural_ranking_full.csv"), index=False)
    sel = df[df["selected"] == True].sort_values(
        ["apolar_bsa_per_res", "contacts"], ascending=False)
    sel.to_csv(os.path.join(args.outdir, "selected_binders.csv"), index=False)

    print(f"\nProcessed {len(df)} structures")
    print(f"Passed all gates: {int(df['passed_gates'].sum())}")
    print(f"  failed confidence: {int((~df['gate_conf']).sum())}")
    print(f"  failed chemistry : {int((~df['gate_chem']).sum())}")
    print(f"  failed site      : {int((~df['gate_site']).sum())}")
    print(f"All three pockets occupied: {int((df['pockets_occupied']==3).sum())}")
    if bench is not None:
        print(f"SELECTED (meets-or-exceeds benchmark): {int((df['selected']==True).sum())}")
    else:
        print(f"SELECTED (top-{args.top_pct:g}% intersection): {int((df['selected']==True).sum())}")
    if len(sel):
        cols = ["design","n","peptide","net_charge","apolar_bsa","apolar_bsa_per_res",
                "contacts","pockets_occupied","iface_plddt","site_dist"]
        print("\nTop selected binders:")
        print(sel[cols].head(15).round(1).to_string(index=False))
    print(f"\nWrote -> {args.outdir}/structural_ranking_full.csv and selected_binders.csv")


if __name__ == "__main__":
    main()
