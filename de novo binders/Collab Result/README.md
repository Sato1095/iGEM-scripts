# Generation — where the 2,000 candidates came from

Everything in this folder is output from the design step, run in Colab. It
is what the selection funnel in `../Filter Pipeline.py/` starts from.

## What we ran

RFdiffusion generated **500 backbones** against the MDM2 p53-binding cleft.
ProteinMPNN designed **4 sequences** onto each, and AlphaFold Multimer
folded every one — **2,000 candidates** in total. That is the pool every
count in the funnel is measured against.

## Files

| | |
|---|---|
| `input.pdb` | the RFdiffusion target: MDM2, chain A, numbered from residue 25 |
| `design.fasta` | the 2,000 candidates, with their AlphaFold metrics in the headers |
| `mpnn_results.csv` | the same 2,000 as a table, with full sequences |
| `diffusion.ipynb` | the ColabDesign notebook the generation ran in |

`design.fasta` is what stage 1 reads. Each record pairs a metric-carrying
header with a `receptor/peptide` sequence line:

```
>design:0 n:0|mpnn:1.381|plddt:0.867|i_ptm:0.800|i_pae:4.703|rmsd:1.423
ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMT.../GLEEEWQKLIERVEK
```

`design.fasta` and `mpnn_results.csv` hold the same data in two formats;
stage 1 will read either. Some scripts refer to the FASTA as
`design_fasta.txt` — same file.

## A caveat on the notebook

`diffusion.ipynb` is the stock ColabDesign notebook with its widget values
untouched: `contigs = "100"`, `num_designs = 1`, `num_seqs = 8`. None of
those are what we ran — the pool is 500 designs at 4 sequences each, against
the MDM2 cleft rather than an unconditional 100-residue backbone.

So the notebook records the *method*, not the run. The parameters we
actually used are not captured anywhere here; the only evidence of them is
the shape of `design.fasta` itself.
