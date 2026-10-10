# Selection funnel — de novo MDM2 binders

How we got from 2,000 generated candidates to one leading warhead. This is a
record of what we ran and what each step caught; it is not a tool for anyone
else to use, and the scripts were written for our own pipeline only.

The three candidates named throughout are **Z1 (42.3)**, **Z2 (267.3)** and
**Z3 (131.2)**.

---

## The funnel, in order

The funnel has three stages but two scripts. Stages 1 and 2 are the same
script read two ways, because they compute the same thing and differ only in
which columns decide the pick. We have kept it that way rather than split it
to match the stage count.

### Stage 1 — composite Z-score · `z-score_filter.py`

Ranked the pool on AlphaFold Multimer confidence alone, as the signed sum
`z_iptm + z_plddt − z_ipae − z_rmsd`. Top pick: **Z1 (42.3)**.

**What it could not catch:** nothing here looks at chemistry, and nothing
looks at the pharmacophore. The metrics score predicted geometry.

### Stage 2 — confidence filter + pharmacophore flags · `z-score_filter.py`

Same script. Applied the confidence tier as a hard filter, then flagged each
survivor for a tryptophan and for the p53-style i, i+4, i+7 anchor register —
the spacing of Phe19 / Trp23 / Leu26 in 1YCR.

This is where the pick changes. Z1 tops the composite ranking but carries
*neither* flag; Z2 is second and carries both. Reading the flag columns moves
the pick to **Z2 (267.3)**.

**What it could not catch:** that Z2 carries net charge −8. An acidic helix
is chemically implausible for a hydrophobic groove, but it satisfies the
pharmacophore on sequence and AlphaFold predicts it confidently. No
sequence-level test in this stage can see the problem.

### Stage 3 — interface geometry and chemistry · `rank_binders.py`

Dropped sequence entirely and scored atomic coordinates: apolar buried
surface area per residue, sub-pocket occupancy, intermolecular contacts, plus
three gates (interface pLDDT, net charge, binding-site position). Selection
is an intersection of independent axes against the PMI benchmark, not a
weighted score.

```
2,000 candidates
1,275  fill all three sub-pockets
  409  pass all three gates
   18  meet or exceed PMI on both axes
```

**This is the stage that caught Z2**, which fails the chemistry gate on its
−8 charge. The winner is **Z3 (131.2)**, `MKPEEVGEMLLSVLK` — 49 contacts,
3/3 pockets, net charge −1, interface pLDDT 89.46, apolar BSA 42.02 Å²/residue.

That the funnel needed three attempts, each because the previous one missed
something, is the point of keeping all of it. A single cleaner script would
not show that Z2 was ever the leading candidate.

---

## Files

| file | role |
|---|---|
| `z-score_filter.py` | stages 1 and 2 |
| `rank_binders.py` | stage 3 |
| `z-score_results.csv` | stage 1/2 shortlist, 64 rows (`stringent` tier) |
| `gate_filter_results_full.csv` | stage 3 table, all 2,000 rows |
| `gate_filter_top_binders.csv` | stage 3 shortlist, the 18 selected |
| `pockets_chimerax.json` | the three sub-pockets, from ChimeraX contacts on 1YCR ≤4.5 Å, with the −24 renumbering offset to design numbering |
| `bench_PMI_norm.pdb` | PMI–MDM2 benchmark, chains normalised to A=receptor / B=peptide |
| `make_verdict_figures.py` | funnel diagram + selection landscape |
| `make_distribution_figure.py` | the 600 dpi / vector version of the landscape figure |
| `analyze_designs.py` | context plots; **not** a funnel stage |

The two stage 3 CSVs are `rank_binders.py`'s own `structural_ranking_full.csv`
and `selected_binders.csv`, committed under the names we archived them as.

Stage 3 also needs the 2,000 AlphaFold Multimer PDBs, published alongside
this directory and passed with `--pdbdir`. They are not in here because of
their number, not because they are unavailable.

Figures are not kept here. They go through the iGEM uploads tool to
`static.igem.wiki`.

---

## What we could and could not re-run

Checked during the October 2026 documentation pass, on Python 3.13 with
pandas 3.0.5, numpy 2.5.3, matplotlib 3.11.2, biopython 1.88.

**Reproduced exactly:**

- `z-score_results.csv`, from `../Collab Result/design.fasta` with
  `--tier stringent --anchor-motif`. Byte-identical.
- Both figure scripts run and produce their figures.
- `analyze_designs.py` runs; its outputs are not committed, so there was
  nothing to compare against.

**Not verified:** `rank_binders.py` was not re-run. The design PDB set was not
yet in place when the pass ran, and `freesasa` would not build in that
container. The two stage 3 CSVs are therefore taken as given, not
reconfirmed — every number quoted above was read back out of the committed
table rather than recomputed from structures.

With the PDBs published alongside and `freesasa` installed, this one is
re-runnable, and it is the only part of the funnel still owing that check.
Reproducing the committed tables needs `--site-tol 14`, not the stale 8.0
default.

---

## Known inconsistencies, not fixed

Left alone deliberately — correcting any of them would change what a script
does or what a figure says.

- **`rank_binders.py --site-tol` defaults to 8.0, but the committed tables
  were made with 14.** Pass `--site-tol 14` to reproduce them. The default is
  stale; `make_distribution_figure.py` calls 14 the "corrected" value and
  `README.md` at the repository root passes it.
- **`analyze_designs.py` ranks 42.3 eighth** while labelling it "best /
  top-ranked", because it scores the unfiltered pool rather than the filtered
  survivors. Its ranks are not the funnel's ranks.
- **`analyze_designs.py`'s third pick is 130.3, not Z3.** `MEEFKKLWEELMKKA`
  is a real and distinct design — it passes the stage 3 gates but is not
  selected — and it predates the Z1/Z2/Z3 naming.
- **Net charge is computed two ways.** Stage 3's gate uses a plain
  K+R−D−E count; `analyze_designs.py` uses Biopython's `charge_at_pH(7.4)`.
  For Z2 that is −8.0 against −8.74. The gate acted on the count.
- **The PMI benchmark values are typed into both figure scripts**
  (BSA 40.0, 35 contacts) rather than scored from `bench_PMI_norm.pdb`, so
  they will not follow if the benchmark is re-scored. The BSA figures are not
  settled on the wiki and should not be quoted as fixed reference values.
- **Figure 1's funnel counts are literal strings**, so they do not track the
  data. As committed they agree with the table.
- **Stage 2's survivor count is quoted inconsistently.** The funnel figure
  says 189, which is the `standard` tier; the committed CSV is the
  `stringent` tier, 64 rows. Both are real runs.
- **`pockets.json` and `derive_pockets.py`** are referenced in a few places
  but have never been in this repository.
- Dead code kept on the record rather than deleted: `sasa_areas()` in
  `rank_binders.py`, `pct()` in `analyze_designs.py`, and the gate
  recomputation in `make_distribution_figure.py`, whose result nothing reads.
- **The confidence-tier thresholds have no recorded provenance.** They are
  not attributed to a reference anywhere in the project.
