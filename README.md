# AETh3R — dry lab analysis

iGEM Thessaloniki 2026. Wiki: **https://2026.igem.wiki/thessaloniki**

These are the analysis scripts behind the Model page. They are a **record of
what we did**, kept so the numbers on the wiki can be traced back to the code
that produced them. They were written for our own pipeline and are not a
software tool for other people to use — we are not entering the Software & AI
Village and are not nominating for Best Software Tool.

---

## The project

**AETh3R** is a peptide-PROTAC (p-PROTAC) for non-small cell lung cancer,
designed to **degrade** its targets rather than inhibit them.

- **MDM2** — the negative regulator of p53. Degrading it releases
  p53-mediated tumour suppression. This is the arm the dry lab designed
  binders for.
- **KRAS** — mutated in ~30% of NSCLC.

The degrader recruits the **VHL** E3 ligase through **pHIF1**
(LA-Hyp-Y-Hle-P), a literature-derived peptide from Chen et al., joined to
the warhead by a glycine–serine linker. pHIF1's two non-canonical residues —
4-hydroxyproline and homoleucine — constrain every downstream tool and are
the source of most of the workarounds in the MD scripts.

---

## The pipeline, in order

### 1 · Warhead — de novo MDM2 binders

**Interface analysis.** Three hydrophobic sub-pockets mapped from 1YCR, the
native p53–MDM2 complex, around the p53 anchors Phe19, Trp23 and Leu26.

**Generation** → `de novo binders/Collab Result/`
RFdiffusion (500 backbones) → ProteinMPNN (4 sequences each) → AlphaFold
Multimer. **2,000 candidates.**

**Three-stage selection funnel** → `de novo binders/Filter Pipeline.py/`

| stage | what it ranks on | top candidate |
|---|---|---|
| 1 · composite Z-score | AlphaFold confidence metrics | **Z1 (42.3)** |
| 2 · confidence filter + pharmacophore flags | the above, plus Trp and the i, i+4, i+7 anchor register | **Z2 (267.3)** |
| 3 · interface geometry and chemistry | atomic coordinates: apolar BSA, pocket occupancy, contacts, net charge | **Z3 (131.2)** |

**Why the sequence matters.** Z2 carries net charge −8 — chemically
implausible for a hydrophobic groove — but confidence metrics rank it highly
because they score predicted geometry and are blind to chemistry. Stage 3
exists to catch that, and does: Z2 fails at the chemistry gate.

```
2,000 candidates → 1,275 fill 3/3 pockets → 409 pass all gates → 18 meet or exceed PMI
```

Stages 1 and 2 are one script run two ways; stage 3 is a second script. They
are kept apart deliberately. That the funnel evolved through three attempts,
each written because the previous one missed something, is the argument the
Model page makes — merging them would destroy the evidence for it.

### 2 · Assembly and validation → `de novo binders/MD simulations/`

Four constructs — **PMI-PROTAC, PZ1, PZ2, PZ3** — simulated as ternary
complexes with MDM2 and **VBC** (VHL + ElonginB + ElonginC).

**Headline result:** PZ3 carries the best MDM2 warhead and makes the most
contacts with the target, but forms the *weakest* induced MDM2–VHL interface.
PZ1 is the reverse, and is the leading degrader candidate. **Target binding
and degradation efficiency are separable.**

### 3 · Pathways — Boolean model → `de novo binders/Pathways/`

A separate strand, asking whether degrading MDM2 and KRAS together kills an
A549 cell that survives either alone.

```
117 edges → 29 nodes / 99 edges after curation → 29 rules
→ 1 provable fixed point per condition → 32 variants, 4 pass, 4 show synergy
```

---

## Layout

```
de novo binders/
  Collab Result/       generation: the 2,000-candidate pool and its inputs
  Filter Pipeline.py/  the three-stage selection funnel, and the figures
  MD simulations/      GROMACS setup for the inhibitor, ternary and VHL systems
  Pathways/            network build, curation, rules, attractors, variant screen
```

Each directory has its own README covering what the scripts did, what they
took in, what they produced, and — for the funnel — what each stage could not
catch and which later stage addressed it.

---

## Environment

**Python.** The pathways pipeline is standard library only, Python 3.8+. The
selection funnel needs more:

| script | needs |
|---|---|
| `z-score_filter.py` | standard library only |
| `rank_binders.py` | biopython, freesasa, numpy, pandas |
| `analyze_designs.py` | biopython, matplotlib, pandas |
| `make_verdict_figures.py`, `make_distribution_figure.py` | matplotlib, numpy, pandas |

**Molecular dynamics.** GROMACS 2022, AMBER99SB-ILDN, TIP3P water. pHIF1 as a
free peptide is parametrised with GAFF2 via ACPYPE, which needs AmberTools 24
and Open Babel. The SLURM scripts carry one cluster's module names and
partitions and will need adapting elsewhere.

**Structures.** `.pdb` inputs are published alongside the scripts in this
repository; the funnel's stage 3 reads the 2,000 AlphaFold Multimer PDBs via
`--pdbdir`.

---

## Things that are easy to get wrong

**Replicates.** One trajectory per system, no replicates. Stated as a
limitation on the wiki; it should stay stated.

**KRAS.** The dry lab designed MDM2 binders. The KRAS strand is an unstapled
SOS1 helix that gave a clear negative result (ipTM 0.27, no RMSD plateau) and
is presented as future work.

**Buried surface area.** The BSA figures — apolar Å² per residue, and the PMI
benchmark — are not yet settled. The numbers in the scripts are what the
selection ran on; do not treat them as fixed reference values in
documentation until the team has resolved how they are presented.

---

## What is not in here

- **Trajectories and trajectory-analysis scripts** — published separately in
  the iGEM GitLab repository.
- **Figures and `.svg` assets** — these go through the iGEM uploads tool onto
  `static.igem.wiki`. The scripts that generate them are here; their output
  is not.
- **Wiki code** — HTML/CSS/JS live in the separate wiki repository.

This repository is developed on GitHub and mirrored into the team's iGEM
GitLab repository, which is what the wiki links to.

---

## Acknowledgements

Documentation and repository organisation were done with the help of
**Claude Code** (Anthropic). The analysis itself, every threshold, every
curation decision and every result is the team's own — the assistance was
limited to writing the READMEs and code comments, reorganising the files, and
checking the scripts against each other and against the wiki. No logic,
threshold, filter, ordering or default parameter was changed in the process.
