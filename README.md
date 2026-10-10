# MDM2 de novo inhibitor — structural ranking pipeline

Scripts and data for ranking de novo MDM2-binding peptides by structure-based
interface quality, benchmarked against PMI. All structural work is on
AlphaFold-Multimer PDBs (chain A = MDM2 receptor, chain B = peptide).

## Dependencies
```
pip install biopython freesasa numpy pandas matplotlib
```

## Files

### Pocket definitions (data)
- `pockets_chimerax.json`  — **ACTIVE / recommended.** Three MDM2 sub-pockets
  (F19/W23/L26) defined from ChimeraX interatomic contacts on 1YCR (<=4.5 A),
  shared wall residues allowed in multiple pockets. This is the default the
  ranker uses.
- `pockets.json`           — earlier version, from geometric nearest-anchor
  centroid assignment. Superseded; kept for provenance. Gives identical
  occupancy calls for the four lead designs.
- `derive_pockets.py`      — regenerates `pockets.json` from 1YCR.pdb
  (the geometric method). Not needed unless re-deriving.

### Ranking (main)
- `rank_binders.py`        — scores every PDB in a folder on apolar BSA (per
  residue), sub-pocket occupancy, and intermolecular contacts; applies
  confidence / chemistry / site gates; selects designs that meet-or-exceed the
  PMI benchmark on all axes (intersection rule).
- `bench_PMI_norm.pdb`     — PMI-MDM2 benchmark structure (chains normalized
  to A=receptor/B=peptide, receptor renumbered to match designs). Pass as
  `--benchmark`.

### Analysis / figures
- `analyze_designs.py`     — sequence-level distribution analysis (composite
  Z-score + physicochemical properties: charge, hydrophobicity, acidic
  fraction) over the FASTA pool. Produces the chemistry/Z-score distribution
  plots. Input: `design_fasta.txt`.
- `make_distribution_figure.py` — the 2-panel "where the designs sit" figure
  (interface-quality axes + chemistry gate). Input: `structural_ranking_full.csv`.

## Run order

1. (optional) regenerate geometric pockets:
   `python3 derive_pockets.py`            # needs 1YCR.pdb; writes pockets.json

2. rank the full pool (uses ChimeraX pockets by default):
   ```
   python3 rank_binders.py \
     --pdbdir /path/to/2000_multimer_pdbs \
     --benchmark bench_PMI_norm.pdb \
     --charge-floor -4 --site-tol 14 --plddt-floor 70 \
     --outdir results_struct
   ```
   Outputs: `results_struct/structural_ranking_full.csv` (all designs, all
   metrics, gate flags) and `selected_binders.csv` (meets-or-exceeds-PMI shortlist).

   PDB filenames must match `design<NN>_n<M>.pdb` to carry the design/n IDs.
   To use the older geometric pockets instead: add `--pockets pockets.json`.

3. distribution figure:
   `python3 make_distribution_figure.py`  # reads structural_ranking_full.csv

4. (independent) sequence-level chemistry/Z-score analysis:
   `python3 analyze_designs.py --fasta design_fasta.txt --outdir results`

## Notes
- BSA = (SASA_A + SASA_B - SASA_AB)/2, apolar atoms (C,S) only, via FreeSASA.
- Selection is an INTERSECTION (all gates AND 3 pockets AND >= PMI on apolar
  BSA/res AND >= PMI on contacts), not a weighted score.
- PMI benchmark = domain-receptor AF3 structure (clean single-site);
  designs = Multimer. This BSA comparison is cross-predictor — stated as a
  caveat in the writeup.
- Everything is measured on predicted structures; MD is the confirmatory step.
