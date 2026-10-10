#!/usr/bin/env python3
"""
make_distribution_figure.py — the two-panel "where the designs sit" figure
for the de novo MDM2 binder strand.

WHAT THIS DID
-------------
Plotted the whole 2,000-candidate pool on the stage 3 selection axes, so
the three named candidates could be shown in context rather than quoted
as bare numbers. It computes nothing: every value is read from the stage 3
table.

  left panel   apolar BSA per residue vs intermolecular contacts, split by
               whether all three sub-pockets are filled, with the PMI
               benchmark quadrant shaded
  right panel  net-charge histogram with the chemistry gate at -4, which
               is the panel that shows why Z2 is rejected

INPUT
-----
    gate_filter_results_full.csv — the stage 3 table, all 2,000 rows.
    This is rank_binders.py's structural_ranking_full.csv, committed here
    under the name the team archived it as.

PRODUCED
--------
    design_distribution.png (600 dpi) and design_distribution.pdf, both
    with transparent backgrounds, in the working directory. Figures are
    not kept in this repository — they go through the iGEM uploads tool to
    static.igem.wiki.

RELATIONSHIP TO make_verdict_figures.py
---------------------------------------
That script's figure 2 is a later, simpler version of this same figure.
Both are kept: this one is the 600 dpi / vector-PDF version with the
reference-envelope shading, and it labels the candidates by role rather
than by Z-number. They are not interchangeable, and the two disagree on
nothing except presentation.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Publication-quality matplotlib settings
# ------------------------------------------------------------
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "axes.linewidth": 1.1,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 9.5,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.unicode_minus": True,
})

# ------------------------------------------------------------
# Input / output
# ------------------------------------------------------------
# Both relative to the working directory, so run this from this folder:
# it reads the stage 3 table beside it and writes design_distribution.png
# and design_distribution.pdf here.
INPUT = "gate_filter_results_full.csv"

OUTPUT_BASE = "design_distribution"

df = pd.read_csv(INPUT)

# ------------------------------------------------------------
# Recompute gates consistent with corrected site tolerance = 14 Å
#
# NOTE: nothing below reads gate_site or passed_gates -- neither panel uses
# them, so this block has no effect on the figure. It is kept as the record
# that 14 Å, not rank_binders.py's 8.0 Å default, is the site tolerance the
# committed stage 3 table was produced with. That is also what README.md's
# invocation passes via --site-tol.
# ------------------------------------------------------------
df["gate_site"] = df["site_dist"] <= 14
df["passed_gates"] = (
    df["gate_conf"]
    & df["gate_chem"]
    & df["gate_site"]
)

# ------------------------------------------------------------
# PMI AF3 benchmark — the meets-or-exceeds thresholds stage 3 selects on.
#
# Typed in here rather than scored from bench_PMI_norm.pdb, so these will
# not follow if the benchmark is re-scored. They match the values stage 3
# used for the committed table.
#
# The benchmark is an AF3 domain-receptor structure while the designs are
# AlphaFold Multimer, so the BSA comparison is cross-predictor — a caveat
# the write-up states. The BSA figures are NOT settled on the wiki; do not
# treat 40.0 as a fixed reference value in documentation.
# ------------------------------------------------------------
PMI = {
    "bsa": 40.0,
    "contacts": 35,
    "charge": -1.0,
}

# ------------------------------------------------------------
# Designs to highlight, keyed by (design, n).
#
# These are the three candidates the funnel produced: 42.3 = Z1,
# 267.3 = Z2, 131.2 = Z3. The labels below name them by the role each
# plays in the argument rather than by Z-number, so the legend reads as
# bare "131.2" where the canonical naming is "Z3 (131.2)".
# ------------------------------------------------------------
PICKS = {
    (131, 2): dict(
        label="131.2  (interface-optimal)",
        color="#e6194B",
        marker="*",
        size=520,
    ),
    (42, 3): dict(
        label="42.3  (confidence-optimal)",
        color="#4363d8",
        marker="s",
        size=170,
    ),
    (267, 3): dict(
        label="267.3  (chemistry-rejected)",
        color="#9A6324",
        marker="D",
        size=155,
    ),
}


def row(design, n):
    subset = df[(df["design"] == design) & (df["n"] == n)]

    if subset.empty:
        raise ValueError(
            f"Could not find design={design}, n={n} in the input CSV."
        )

    return subset.iloc[0]


# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    1,
    2,
    figsize=(16.5, 7.2),
    gridspec_kw={"width_ratios": [1.3, 1]},
)

# ============================================================
# PANEL 1
# Per-residue apolar BSA vs contacts
# ============================================================

p3 = df[df["pockets_occupied"] == 3]
pless = df[df["pockets_occupied"] < 3]

# Background population
ax1.scatter(
    pless["apolar_bsa_per_res"],
    pless["contacts"],
    s=10,
    c="#d9d9d9",
    alpha=0.50,
    edgecolors="none",
    label="< 3 pockets",
    rasterized=True,
    zorder=1,
)

ax1.scatter(
    p3["apolar_bsa_per_res"],
    p3["contacts"],
    s=12,
    c="#9ecae1",
    alpha=0.72,
    edgecolors="none",
    label="all 3 pockets",
    rasterized=True,
    zorder=2,
)

# Axis limits
xmin = df["apolar_bsa_per_res"].min() - 1
xmax = df["apolar_bsa_per_res"].max() + 1
ymin = df["contacts"].min() - 2
ymax = df["contacts"].max() + 3

ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)

# ------------------------------------------------------------
# PMI threshold quadrant
# ------------------------------------------------------------
ax1.fill_between(
    [PMI["bsa"], xmax],
    PMI["contacts"],
    ymax,
    color="#c7e9c0",
    alpha=0.38,
    zorder=0,
)

ax1.axvline(
    PMI["bsa"],
    color="#d6604d",
    linestyle="--",
    linewidth=1.8,
    zorder=3,
)

ax1.axhline(
    PMI["contacts"],
    color="#d6604d",
    linestyle="--",
    linewidth=1.8,
    zorder=3,
)

# PMI marker
ax1.scatter(
    [PMI["bsa"]],
    [PMI["contacts"]],
    marker="X",
    s=300,
    c="#d6604d",
    edgecolors="black",
    linewidths=1.0,
    zorder=7,
)

ax1.annotate(
    "PMI benchmark\n(AF3)",
    xy=(PMI["bsa"], PMI["contacts"]),
    xytext=(-85, -8),
    textcoords="offset points",
    fontsize=10.5,
    color="#8f1d1d",
    weight="bold",
    ha="right",
    va="center",
    arrowprops=dict(
        arrowstyle="->",
        color="#8f1d1d",
        lw=1.2,
    ),
)

ax1.annotate(
    "meets-or-exceeds\nPMI + 3 pockets",
    xy=(xmax - 0.5, ymax - 0.8),
    fontsize=10.5,
    color="#276221",
    weight="bold",
    ha="right",
    va="top",
)

# ------------------------------------------------------------
# Highlight designs
# ------------------------------------------------------------
for (design, n), st in PICKS.items():

    r = row(design, n)

    ax1.scatter(
        r["apolar_bsa_per_res"],
        r["contacts"],
        marker=st["marker"],
        s=st["size"],
        c=st["color"],
        edgecolors="black",
        linewidths=1.2,
        zorder=10,
        label=st["label"],
    )

# Labels
ax1.set_xlabel(
    "Apolar buried surface area per residue (Å²)",
    fontsize=13,
    labelpad=8,
)

ax1.set_ylabel(
    "Intermolecular contacts (≤ 5 Å)",
    fontsize=13,
    labelpad=8,
)

ax1.set_title(
    "Interface-quality selection axes  (n = 2000)",
    fontsize=14,
    weight="bold",
    pad=12,
)

ax1.grid(
    alpha=0.18,
    linewidth=0.8,
)

ax1.tick_params(
    direction="out",
    length=5,
    width=1,
)

ax1.legend(
    fontsize=9.5,
    loc="lower left",
    framealpha=0.95,
    edgecolor="#cccccc",
)

# ============================================================
# PANEL 2
# Net charge distribution
# ============================================================

charge_min = df["net_charge"].min()
charge_max = df["net_charge"].max()

bins = np.arange(
    charge_min - 0.5,
    charge_max + 1.5,
    1,
)

ax2.hist(
    df["net_charge"],
    bins=bins,
    color="#cfd8e3",
    edgecolor="#667788",
    linewidth=0.65,
    alpha=0.95,
)

# Chemistry rejection region
ax2.axvspan(
    charge_min - 1,
    -4,
    color="#f4cccc",
    alpha=0.40,
    zorder=0,
)

# Reference envelope
ax2.axvspan(
    -2.2,
    -0.8,
    color="#c7e9c0",
    alpha=0.42,
    zorder=0,
)

# Chemistry gate
ax2.axvline(
    -4,
    color="#333333",
    linestyle=":",
    linewidth=2.0,
    zorder=4,
)

# PMI charge
ax2.axvline(
    PMI["charge"],
    color="#d6604d",
    linestyle="--",
    linewidth=1.8,
    zorder=4,
)

# ------------------------------------------------------------
# Get histogram height for annotation placement
# ------------------------------------------------------------
hist_counts, _ = np.histogram(
    df["net_charge"],
    bins=bins,
)

hist_max = hist_counts.max()

# Gate annotation
ax2.annotate(
    "chemistry gate\n(charge ≥ −4)",
    xy=(-4, hist_max * 0.78),
    xytext=(-4.25, hist_max * 0.88),
    fontsize=10,
    ha="right",
    va="center",
    color="#333333",
    weight="bold",
    arrowprops=dict(
        arrowstyle="->",
        color="#333333",
        lw=1.0,
    ),
)

# Reference envelope annotation
ax2.annotate(
    "reference\nenvelope",
    xy=(-1.5, hist_max * 0.93),
    fontsize=9.5,
    ha="center",
    va="top",
    color="#276221",
    weight="bold",
)

# ------------------------------------------------------------
# Highlight selected designs
# ------------------------------------------------------------
for i, ((design, n), st) in enumerate(PICKS.items()):

    r = row(design, n)

    charge = r["net_charge"]

    ax2.axvline(
        charge,
        color=st["color"],
        linewidth=2.4,
        zorder=5,
    )

    # Alternate vertical positions to reduce overlap
    ypos = hist_max * (0.60 - 0.12 * i)

    ax2.annotate(
        f"{design}.{n}",
        xy=(charge, ypos),
        fontsize=9.5,
        ha="center",
        va="center",
        color=st["color"],
        weight="bold",
        bbox=dict(
            boxstyle="round,pad=0.20",
            fc="white",
            ec=st["color"],
            lw=1.0,
        ),
        zorder=8,
    )

# ------------------------------------------------------------
# Axis labels / title
# ------------------------------------------------------------
ax2.set_xlabel(
    "Net peptide charge (pH 7.4)",
    fontsize=13,
    labelpad=8,
)

ax2.set_ylabel(
    "Number of candidates",
    fontsize=13,
    labelpad=8,
)

ax2.set_title(
    "Chemistry gate — why 267.3 is rejected",
    fontsize=14,
    weight="bold",
    pad=12,
)

ax2.grid(
    alpha=0.18,
    linewidth=0.8,
    axis="y",
)

ax2.tick_params(
    direction="out",
    length=5,
    width=1,
)

# ------------------------------------------------------------
# Overall title
# ------------------------------------------------------------
fig.suptitle(
    "Where the selected designs sit in the 2000-candidate pool",
    fontsize=16,
    weight="bold",
    y=0.995,
)

# Better spacing
fig.tight_layout(
    rect=[0, 0, 1, 0.965],
    pad=1.5,
    w_pad=2.2,
)

# ------------------------------------------------------------
# Save high-resolution PNG + vector PDF with transparent background
# ------------------------------------------------------------
png_file = OUTPUT_BASE + ".png"
pdf_file = OUTPUT_BASE + ".pdf"

# Transparent figure/axes backgrounds
fig.patch.set_alpha(0)

for ax in (ax1, ax2):
    ax.set_facecolor("none")

fig.savefig(
    png_file,
    dpi=600,
    bbox_inches="tight",
    pad_inches=0.08,
    transparent=True,
)

fig.savefig(
    pdf_file,
    bbox_inches="tight",
    pad_inches=0.08,
    transparent=True,
)

plt.close(fig)

print(f"Saved transparent PNG: {png_file}")
print(f"Saved transparent PDF: {pdf_file}")