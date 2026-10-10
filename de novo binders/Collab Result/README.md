# Generation — the 2,000-candidate pool

What went into the selection funnel. Everything here is the output of the
generation step, run in Colab; the funnel itself is in
`../Filter Pipeline.py/`.

## What we ran

RFdiffusion generated 500 backbones against the MDM2 p53-binding cleft.
ProteinMPNN designed 4 sequences onto each, and AlphaFold Multimer folded
every one of them, giving **500 × 4 = 2,000 candidates**. That is the pool
every count in the funnel is measured against.

Confirmed from `design.fasta`: backbone indices 0–499, sequence indices
n = 0–3, 2,000 records.

## Files

| file | what it is |
|---|---|
| `input.pdb` | the RFdiffusion target — MDM2, chain A, 705 atoms, numbered from residue 25 |
| `design.fasta` | the 2,000 candidates. Header carries the AlphaFold metrics, sequence line is `receptor/peptide` |
| `mpnn_results.csv` | the same 2,000 rows as a table, with the full sequence column |
| `diffusion.ipynb` | the ColabDesign RFdiffusion notebook |

`design.fasta` is the input to stage 1. It is the file the scripts refer to
as `design_fasta.txt`.

A record looks like this:

```
>design:0 n:0|mpnn:1.381|plddt:0.867|i_ptm:0.800|i_pae:4.703|rmsd:1.423
ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMT.../GLEEEWQKLIERVEK
```

`design.fasta` and `mpnn_results.csv` hold the same data in two formats.
Stage 1 reads either.

## The notebook does not record our parameters

`diffusion.ipynb` is the stock ColabDesign notebook
(`sokrypton/ColabDesign`, `rf/examples/diffusion.ipynb`) with its widget
values untouched: `contigs = "100"`, `name = "test"`, `num_designs = 1`,
`num_seqs = 8`. None of those are what we ran — the pool is 500 designs at
4 sequences each, against the MDM2 cleft rather than a 100-residue
unconditional backbone.

So the notebook records the *method*, not the run. The actual generation
parameters are not captured anywhere in this repository; the only durable
evidence of them is the shape of `design.fasta` itself. Worth fixing by
hand if the team still has the Colab session.
