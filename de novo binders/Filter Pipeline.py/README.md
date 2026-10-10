# Selection funnel — picking the MDM2 binder

How we got from 2,000 generated candidates down to one leading warhead,
**Z3 (131.2)**.

The short version: ranking designs on how confident AlphaFold is about them
gets you a long way, and then stops working. The top-ranked candidate by
confidence has no tryptophan and no p53-like anchor spacing. The candidate
that *does* have both carries a net charge of −8, which is not something you
put in a hydrophobic groove. Neither problem is visible to the metrics that
chose them. Finding that out is what the third stage is for.

---

## The three stages

Three stages, two scripts. Stages 1 and 2 are the same script read two
different ways — they compute the same thing and differ only in which
columns decide the pick. We have kept it that way rather than split it to
match the stage count.

### Stage 1 — rank on confidence · `z-score_filter.py`

Score every candidate on AlphaFold's own confidence metrics, combined into
one number: `z_iptm + z_plddt − z_ipae − z_rmsd`.

Top pick: **Z1 (42.3)**.

*What this misses:* everything about chemistry, and the pharmacophore. These
metrics describe how sure the model is about a predicted shape.

### Stage 2 — add the pharmacophore · `z-score_filter.py`

Same script. Apply the confidence thresholds as a hard filter, then flag
each survivor for a tryptophan and for the p53 anchor register — Phe19,
Trp23 and Leu26 sit on one face of the helix at spacing i, i+4, i+7, and
that is the pattern we look for.

This is where the pick changes. Z1 tops the ranking but carries *neither*
flag; the second-place candidate carries both. Reading the flag columns
moves the pick to **Z2 (267.3)**.

*What this misses:* that Z2 has a net charge of −8. The pharmacophore is
satisfied on paper, AlphaFold is confident, and nothing at the sequence
level sees the problem.

### Stage 3 — look at the structure · `rank_binders.py`

Drop sequence entirely and measure the actual interface: how much
hydrophobic surface gets buried per residue, whether all three sub-pockets
are filled, how many contacts form, and — finally — what the peptide's net
charge is. Three gates (confidence, chemistry, binding site), then an
intersection against the PMI benchmark rather than a weighted score, so no
single axis can carry a design through on its own.

```
2,000 candidates
1,275  fill all three sub-pockets
  409  pass all three gates
   18  meet or exceed PMI on both axes
```

**This is the stage that catches Z2**, which fails the chemistry gate on its
charge. The winner is **Z3 (131.2)**, `MKPEEVGEMLLSVLK` — 49 contacts, all
three pockets, net charge −1, interface pLDDT 89.5.

That the funnel took three attempts, each written because the last one
missed something, is the point of keeping all of it. A single tidier script
would not show that Z2 was ever the leading candidate.

---

## Files

| | |
|---|---|
| `z-score_filter.py` | stages 1 and 2 |
| `rank_binders.py` | stage 3 |
| `z-score_results.csv` | the stage 1/2 shortlist |
| `gate_filter_results_full.csv` | the stage 3 table, all 2,000 rows |
| `gate_filter_top_binders.csv` | the 18 that meet or exceed PMI |
| `pockets_chimerax.json` | the three sub-pockets, from ChimeraX contacts on 1YCR |
| `bench_PMI_norm.pdb` | the PMI–MDM2 benchmark structure |
| `make_verdict_figures.py` | the funnel diagram and the selection landscape |
| `make_distribution_figure.py` | the high-resolution version of the landscape figure |
| `analyze_designs.py` | context plots — **not** a funnel stage |

Stage 3 also reads the 2,000 AlphaFold Multimer PDBs, published alongside
this directory and passed with `--pdbdir`. The two stage 3 CSVs are that
script's own output, committed under the names we archived them as.

Figures are not kept here — they go through the iGEM uploads tool to
`static.igem.wiki`.

---

## If you re-run any of this

A few things will bite otherwise.

- **Stage 3 needs `--site-tol 14`.** The script's own default is 8.0 and is
  stale; the committed results were produced at 14.
- **Stage 1/2 composite scores are relative to whichever tier you run.** The
  z-scores are computed over the candidates that survive the filter, not the
  whole pool, so the same design scores differently under a different tier.
  The ranking at the top is unaffected — Z1 comes first either way — but a
  composite number only means something next to its tier.
- **`analyze_designs.py` is a side analysis, not a stage.** It scores the
  unfiltered pool, so its ranks are not the funnel's ranks (Z1 comes 8th
  there), and the third design it marks, 130.3, is a real design but not one
  of the three the write-up names. Don't quote ranks from it.
- **Net charge is computed two ways.** Stage 3's gate counts K+R−D−E, which
  gives Z2 the −8 the write-up quotes. `analyze_designs.py` uses Biopython's
  pH-aware calculation and gets −8.74 for the same peptide.
- **The PMI benchmark values are typed into the figure scripts** rather than
  re-scored from the structure, so they will not follow if the benchmark is
  rescored. The BSA numbers are not settled — don't treat them as fixed
  reference values.
