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

### 1 · Designing the warhead → `de novo binders/`

**Interface analysis.** Three hydrophobic sub-pockets mapped from 1YCR, the
native p53–MDM2 complex, around the p53 anchors Phe19, Trp23 and Leu26.

**Generation** (`Collab Result/`). RFdiffusion produced 500 backbones,
ProteinMPNN designed 4 sequences onto each, and AlphaFold Multimer folded
every one — **2,000 candidates**.

**Selection** (`Filter Pipeline.py/`), in three stages:

| stage | ranks on | top candidate |
|---|---|---|
| 1 · composite Z-score | AlphaFold confidence metrics | **Z1 (42.3)** |
| 2 · + pharmacophore flags | the above, plus Trp and the i, i+4, i+7 anchor register | **Z2 (267.3)** |
| 3 · interface geometry and chemistry | atomic coordinates: apolar BSA, pocket occupancy, contacts, net charge | **Z3 (131.2)** |

Each stage exists because the one before it missed something. Z1 tops the
confidence ranking but has no tryptophan and no anchor register. Z2 has
both — and a net charge of −8, which is not something you put in a
hydrophobic groove. Neither problem is visible to the metrics that chose
them, which is why the third stage reads atomic coordinates instead.

```
2,000 candidates → 1,275 fill 3/3 pockets → 409 pass all gates → 18 meet or exceed PMI
```

Stages 1 and 2 are one script run two ways; stage 3 is a second script.
They are kept apart deliberately — that the funnel evolved through three
attempts is the argument the Model page makes, and merging them would
destroy the evidence for it.

### 2 · Simulating the assembled degraders → `de novo binders/MD simulations/`

Four constructs — **PMI-PROTAC, PZ1, PZ2, PZ3** — simulated as ternary
complexes with MDM2 and **VBC** (VHL + ElonginB + ElonginC). GROMACS 2022,
AMBER99SB-ILDN, 100 ns each; VHL with pHIF1 alone was run to 150 ns.

**The result that matters:** PZ3 carries the best MDM2 warhead and makes the
most contacts with the target, but forms the *weakest* induced MDM2–VHL
interface. PZ1 is the reverse, and is the leading degrader candidate.
**Target binding and degradation efficiency are separable** — a good binder
is not automatically a good degrader.

### 3 · Asking whether the combination works → `de novo binders/Pathways/`

A separate strand, and a different question: does degrading MDM2 and KRAS
*together* kill an A549 cell that survives either one alone? Answered with a
Boolean model of apoptotic commitment, built from OmniPath and curated by
hand.

```
117 edges → 29 nodes / 99 edges after curation → 29 rules
→ 1 provable fixed point per condition → 32 rule variants screened
```

**The result:** 4 of the 32 variants satisfy every pre-registered criterion.
Separately, 4 show true synergy — the combination killing where neither
single agent does. MDM2 degradation alone never kills, in any variant.

Those two sets of four overlap in only two members, and the reason is the
honest difficulty of the whole strand. Synergy requires BCL-XL to act as a
coincidence detector *and* PUMA to be p53-driven rather than also
FOXO3-driven. The criteria can distinguish the first; they cannot
distinguish the second, because FOXO3 is off in every condition they
constrain. So two of the four acceptable rule sets predict synergy and two
do not, and the data cannot say which pair is right. No dataset in hand
tests the combination either, so this is a prediction, not a result.

---

## Layout

```
de novo binders/
  Collab Result/       generation: the 2,000-candidate pool and its inputs
  Filter Pipeline.py/  the three-stage selection funnel, and the figures
  MD simulations/      GROMACS setup for the inhibitor, ternary and VHL systems
  Pathways/            network build, curation, rules, attractors, variant screen
```

Each folder has its own README explaining that part of the work in more
detail.

---

## Environment

**Python.** The pathways pipeline runs on the standard library alone
(Python 3.8+). The selection funnel needs more: `z-score_filter.py` is also
standard-library-only, the figure scripts and `analyze_designs.py` need
matplotlib, pandas, numpy and biopython, and `rank_binders.py` additionally
needs freesasa.

**Molecular dynamics.** GROMACS 2022, AMBER99SB-ILDN, TIP3P water. pHIF1 as
a free peptide is parametrised with GAFF2 via ACPYPE, which needs
AmberTools 24 and Open Babel. The SLURM scripts carry one cluster's module
names and partitions, so they will need adapting elsewhere.

**Structures.** The `.pdb` inputs are published alongside the scripts,
including the 2,000 AlphaFold Multimer designs that stage 3 scores.

---

## Limitations we state openly

**One trajectory per system, no replicates.** Stated as a limitation on the
wiki, and it should stay stated.

**The Boolean model is not validated.** Its criteria were built from the
same proteomics that prompted us to examine those rules, so passing them
constrains the model rather than confirming it. Nothing in hand tests the
KRAS arm or the combination.

**The KRAS strand of the structural work gave a negative result.** It is an
unstapled SOS1 helix with ipTM 0.27 and no RMSD plateau, presented as future
work. The binders the dry lab designed are for MDM2.

**The buried surface area figures are not settled.** The numbers in the
scripts are what the selection actually ran on, but how they are presented
is still being resolved, so they should not be quoted as fixed reference
values yet.

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
limited to writing the READMEs and code comments, reorganising the files,
and checking the scripts against each other. No logic, threshold, filter,
ordering or default parameter was changed in the process.
