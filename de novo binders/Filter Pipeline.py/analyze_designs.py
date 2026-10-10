#!/usr/bin/env python3
"""
analyze_designs.py
==================
Context plots for the de novo MDM2 binder pool. NOT a stage of the
selection funnel — nothing downstream consumes its output.

WHAT THIS DID
-------------
Showed where the chosen designs sit in the full 2,000-candidate
distribution, in two spaces at once: the AlphaFold-metric space the
funnel ranks on, and an independent sequence-level physicochemical space
it is blind to. The argument it was written to make is that a candidate
can rank highly on predicted geometry while being chemically wrong for a
hydrophobic groove — the same argument stage 3 (rank_binders.py) later
settled from structures rather than from sequence.

It computes:
  1. the composite Z-score for all candidates, and
  2. physicochemical descriptors (GRAVY, net charge at pH 7.4,
     aromaticity, hydrophobic-anchor fraction, acidic fraction),
then overlays the picks and the reference peptides (PMI, PDI, p53 17-29)
on the distributions.

INPUT
-----
    Collab Result/design.fasta — the same 2,000-record file stage 1 reads.
    (Referenced as design_fasta.txt in the usage line below; same file.)

PRODUCED
--------
    candidates_annotated.csv, picks_summary.csv, reference_physchem.csv,
    and eight PNGs, into --outdir. None of these are committed.

TWO THINGS TO KNOW BEFORE QUOTING IT
------------------------------------
1. ITS RANKS ARE NOT THE FUNNEL'S RANKS. The composite Z formula here
   (z_iptm + z_plddt - z_ipae - z_rmsd) is the same signed sum
   z-score_filter.py uses, but it is computed over the UNFILTERED pool of
   2,000, whereas the funnel scores only the candidates that survive the
   confidence filter. Re-centring on the survivors changes the order:

       unfiltered (this script)   rank 1 = 366.3, rank 2 = 267.3,
                                  rank 8 = 42.3
       filtered   (stage 1)       rank 1 = 42.3,  rank 2 = 267.3

   So Z1 (42.3) is the top composite pick only after the confidence
   filter. This script labels it "best / top-ranked" while its own table
   puts it 8th. Do not cite a rank from here as the funnel's rank.

2. ITS THIRD PICK IS 130.3, NOT Z3. PICKS below names
   (130, 3) = MEEFKKLWEELMKKA, under the older "full mimic / minimal
   mimic / best" framing. That is a real, distinct design — it is in the
   stage 3 table with net charge 0, 41 contacts, 3/3 pockets, passing the
   gates but not selected — and it is NOT Z3 (131.2, MKPEEVGEMLLSVLK).
   This script predates the Z1/Z2/Z3 naming and was never updated to it.

   Net charges here also differ from stage 3's: this script uses
   Biopython's charge_at_pH(7.4) (42.3 -> -2.75, 267.3 -> -8.74), while
   the chemistry gate uses a plain K+R-D-E count (-2.0, -8.0). The gate
   acted on the latter.

Usage
-----
    python3 analyze_designs.py --fasta design_fasta.txt --outdir results

Input format
------------
FASTA where each header carries the metrics and each sequence line is
"RECEPTOR/PEPTIDE", e.g.

    >design:42 n:3|mpnn:1.355|plddt:0.913|i_ptm:0.828|i_pae:3.726|rmsd:0.412
    ETLVRPKP.../SVEENMKEFEKMILG
"""

import argparse
import os
import re
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# --------------------------------------------------------------------------- #
# Reference peptides (controls). These are NOT in the design pool; they are
# overlaid on the distributions as benchmarks. Sequences from the project PDF.
# --------------------------------------------------------------------------- #
REFERENCES = {
    "PMI":        "TSFAEYWNLLSP",
    "PDI":        "LTFEHYWAQLTS",
    "p53(17-29)": "SQETFSDLWKLLP",
}

# The three designs this script marks, keyed by (design_index, seq_index).
# These labels are from the write-up as it stood when this was written.
# They are NOT the current Z1/Z2/Z3 naming, and (130, 3) is not Z3 — see
# the header.
PICKS = {
    (42, 3):  "42.3  (best / top-ranked)",
    (267, 3): "267.3 (minimal mimic)",
    (130, 3): "130.3 (full mimic)",
}

# Kyte-Doolittle is the standard hydropathy scale; higher = more hydrophobic.
HEADER_RE = re.compile(
    r">design:(?P<design>\d+)\s+n:(?P<n>\d+)\|"
    r"mpnn:(?P<mpnn>[-\d.]+)\|"
    r"plddt:(?P<plddt>[-\d.]+)\|"
    r"i_ptm:(?P<i_ptm>[-\d.]+)\|"
    r"i_pae:(?P<i_pae>[-\d.]+)\|"
    r"rmsd:(?P<rmsd>[-\d.]+)"
)


def parse_fasta(path):
    """Return a DataFrame with one row per candidate."""
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line
            else:
                if header is None:
                    continue
                m = HEADER_RE.match(header)
                if not m:
                    print(f"WARNING: could not parse header: {header}", file=sys.stderr)
                    header = None
                    continue
                d = m.groupdict()
                # sequence line is RECEPTOR/PEPTIDE
                peptide = line.split("/")[-1]
                rows.append({
                    "design": int(d["design"]),
                    "n": int(d["n"]),
                    "mpnn": float(d["mpnn"]),
                    "plddt": float(d["plddt"]),
                    "i_ptm": float(d["i_ptm"]),
                    "i_pae": float(d["i_pae"]),
                    "rmsd": float(d["rmsd"]),
                    "peptide": peptide,
                    "length": len(peptide),
                })
                header = None
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("No records parsed — check the input format.")
    return df


def add_physchem(df):
    """Add sequence-level physicochemical descriptors."""
    def props(seq):
        # ProtParam needs standard residues; all designs here are standard AAs.
        pa = ProteinAnalysis(seq)
        gravy = pa.gravy()                       # Kyte-Doolittle mean hydropathy
        charge = pa.charge_at_pH(7.4)            # net charge near physiological pH
        aromatic = pa.aromaticity()              # fraction F+Y+W
        # fraction of the canonical MDM2 anchor residue types (F/W/L/I/M/V)
        hydro_anchor = sum(seq.count(a) for a in "FWLIMV") / len(seq)
        # fraction acidic (D/E) — the property that flags the E-rich case
        acidic = sum(seq.count(a) for a in "DE") / len(seq)
        return pd.Series({
            "gravy": gravy,
            "net_charge": charge,
            "aromaticity": aromatic,
            "hydrophobic_anchor_frac": hydro_anchor,
            "acidic_frac": acidic,
        })
    return df.join(df["peptide"].apply(props))


def zscore(series):
    return (series - series.mean()) / series.std(ddof=0)


def add_composite(df):
    """Composite Z as defined in the project: Z_iptm + Z_plddt - Z_ipae - Z_rmsd.

    Computed over every row passed in — i.e. the whole 2,000-candidate pool,
    with no confidence filter first. z-score_filter.py applies the filter
    before scoring, so its composite values and ranks differ from these.
    See the header: 42.3 is 8th here and 1st there.
    """
    df["z_iptm"] = zscore(df["i_ptm"])
    df["z_plddt"] = zscore(df["plddt"])
    df["z_ipae"] = zscore(df["i_pae"])
    df["z_rmsd"] = zscore(df["rmsd"])
    df["composite_z"] = df["z_iptm"] + df["z_plddt"] - df["z_ipae"] - df["z_rmsd"]
    df["composite_rank"] = df["composite_z"].rank(ascending=False, method="min").astype(int)
    return df


def label_picks(df):
    df["pick"] = ""
    for (d, n), lab in PICKS.items():
        mask = (df["design"] == d) & (df["n"] == n)
        df.loc[mask, "pick"] = lab
    return df


def ref_frame():
    """Physicochemical properties for the reference peptides."""
    rows = []
    for name, seq in REFERENCES.items():
        pa = ProteinAnalysis(seq)
        rows.append({
            "name": name, "peptide": seq, "length": len(seq),
            "gravy": pa.gravy(),
            "net_charge": pa.charge_at_pH(7.4),
            "aromaticity": pa.aromaticity(),
            "hydrophobic_anchor_frac": sum(seq.count(a) for a in "FWLIMV") / len(seq),
            "acidic_frac": sum(seq.count(a) for a in "DE") / len(seq),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
PICK_COLORS = {
    "42.3  (best / top-ranked)": "#d62728",
    "267.3 (minimal mimic)":     "#9467bd",
    "130.3 (full mimic)":        "#2ca02c",
}


def _annotate_picks(ax, df, xcol, ycol):
    picks = df[df["pick"] != ""]
    for _, r in picks.iterrows():
        ax.scatter(r[xcol], r[ycol], s=110, marker="*",
                   color=PICK_COLORS.get(r["pick"], "black"),
                   edgecolor="black", linewidth=0.6, zorder=5)
        ax.annotate(r["pick"].split()[0], (r[xcol], r[ycol]),
                    textcoords="offset points", xytext=(6, 4),
                    fontsize=8, weight="bold")


def _annotate_refs(ax, refs, xcol, ycol):
    if xcol not in refs.columns or ycol not in refs.columns:
        return
    for _, r in refs.iterrows():
        ax.scatter(r[xcol], r[ycol], s=90, marker="D",
                   facecolor="white", edgecolor="black", linewidth=1.2, zorder=6)
        ax.annotate(r["name"], (r[xcol], r[ycol]),
                    textcoords="offset points", xytext=(6, -10),
                    fontsize=8, style="italic")


def plot_hist(df, col, outpath, refs=None, title=None, invert_good=False):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.hist(df[col], bins=50, color="#cfd8e3", edgecolor="#7f8c9b", linewidth=0.4)
    for _, r in df[df["pick"] != ""].iterrows():
        ax.axvline(r[col], color=PICK_COLORS.get(r["pick"], "black"),
                   lw=1.8, label=r["pick"])
    if refs is not None and col in refs.columns:
        for _, r in refs.iterrows():
            ax.axvline(r[col], color="black", ls="--", lw=1.0, alpha=0.7)
            ax.annotate(r["name"], (r[col], ax.get_ylim()[1]*0.92),
                        rotation=90, fontsize=7, va="top", ha="right", style="italic")
    ax.set_xlabel(col)
    ax.set_ylabel("Number of candidates")
    ax.set_title(title or f"Distribution of {col} across {len(df)} candidates")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_scatter(df, xcol, ycol, outpath, refs=None, title=None):
    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.scatter(df[xcol], df[ycol], s=8, color="#b8c4d0", alpha=0.5,
               edgecolor="none", label=f"all candidates (n={len(df)})")
    _annotate_picks(ax, df, xcol, ycol)
    if refs is not None:
        _annotate_refs(ax, refs, xcol, ycol)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(title or f"{ycol} vs {xcol}")
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = parse_fasta(args.fasta)
    df = add_physchem(df)
    df = add_composite(df)
    df = label_picks(df)
    refs = ref_frame()

    # ---- save the full annotated table ----
    out_csv = os.path.join(args.outdir, "candidates_annotated.csv")
    df.sort_values("composite_rank").to_csv(out_csv, index=False)

    # ---- percentile summary for the picks ----
    # DEAD CODE: pct() is never called. The rank_pct column below is
    # computed inline from composite_rank instead.
    def pct(col, val, higher_is_better=True):
        p = (df[col] < val).mean() * 100
        return p if higher_is_better else 100 - p

    summary_rows = []
    for _, r in df[df["pick"] != ""].iterrows():
        summary_rows.append({
            "pick": r["pick"],
            "peptide": r["peptide"],
            "composite_z": round(r["composite_z"], 3),
            "composite_rank": int(r["composite_rank"]),
            "rank_pct": round(100 * (1 - (r["composite_rank"]-1)/len(df)), 1),
            "i_ptm": r["i_ptm"], "plddt": r["plddt"],
            "i_pae": r["i_pae"], "rmsd": r["rmsd"],
            "gravy": round(r["gravy"], 3),
            "net_charge": round(r["net_charge"], 2),
            "acidic_frac": round(r["acidic_frac"], 3),
            "hydro_anchor_frac": round(r["hydrophobic_anchor_frac"], 3),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(args.outdir, "picks_summary.csv"), index=False)

    # ---- reference comparison table ----
    refs_out = refs.copy()
    refs_out.to_csv(os.path.join(args.outdir, "reference_physchem.csv"), index=False)

    # ---- plots ----
    plot_hist(df, "composite_z", os.path.join(args.outdir, "hist_composite_z.png"),
              title="Composite Z-score distribution (picks marked)")
    plot_hist(df, "gravy", os.path.join(args.outdir, "hist_gravy.png"), refs=refs,
              title="Hydrophobicity (GRAVY) — picks (solid) vs references (dashed)")
    plot_hist(df, "net_charge", os.path.join(args.outdir, "hist_net_charge.png"), refs=refs,
              title="Net charge at pH 7.4 — picks vs references")
    plot_hist(df, "acidic_frac", os.path.join(args.outdir, "hist_acidic_frac.png"), refs=refs,
              title="Acidic residue fraction (D+E) — picks vs references")

    # structural-space scatter: the two metrics the Z-score leans on most
    plot_scatter(df, "i_pae", "i_ptm", os.path.join(args.outdir, "scatter_ipae_iptm.png"),
                 title="Interface confidence: i_ptm vs i_pae (lower i_pae, higher i_ptm = better)")
    plot_scatter(df, "rmsd", "plddt", os.path.join(args.outdir, "scatter_rmsd_plddt.png"),
                 title="pLDDT vs RMSD")

    # the key figure: structural quality vs chemistry.
    # composite_z (structural) on x, hydrophobicity on y, references overlaid.
    plot_scatter(df, "composite_z", "gravy",
                 os.path.join(args.outdir, "scatter_compositez_gravy.png"),
                 refs=None,
                 title="Structural rank (composite Z) vs hydrophobicity (GRAVY)")
    plot_scatter(df, "composite_z", "acidic_frac",
                 os.path.join(args.outdir, "scatter_compositez_acidic.png"),
                 refs=None,
                 title="Structural rank (composite Z) vs acidic fraction (D+E)")

    # ---- console report ----
    print(f"Parsed {len(df)} candidates from {args.fasta}")
    print(f"Lengths present: {sorted(df['length'].unique())}")
    print("\n=== Chosen designs in context ===")
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(summary.to_string(index=False))
    print("\n=== Reference peptides (physicochemical) ===")
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(refs_out.round(3).to_string(index=False))
    print(f"\nWrote annotated table -> {out_csv}")
    print(f"Wrote plots + summaries -> {args.outdir}/")


if __name__ == "__main__":
    main()
