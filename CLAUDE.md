# iGEM Thessaloniki 2026 — AETh3R · dry lab analysis

Context for Claude Code sessions on this repository. Read this first.

---

## What this repository is

The **dry lab analysis scripts** for AETh3R, our iGEM 2026 project. It is developed here
on GitHub and then mirrored into the team's iGEM GitLab repository, which is what the
wiki links to.

Live wiki: **https://2026.igem.wiki/thessaloniki** — the Model page is the authoritative
write-up of everything in here, and the two must stay consistent. If a script and the
wiki disagree, stop and flag it rather than changing either silently.

### What is NOT in this repository

- **Trajectories and trajectory-analysis scripts** — published separately in the iGEM
  GitLab repository. Nothing here produces the induced-contact counts or RMSD plateaus.
- **Figures and `.svg` assets** — these must go through the iGEM uploads tool onto
  `static.igem.wiki`. They do not belong here.
- **Wiki code** — HTML/CSS/JS live in the separate wiki repository.

Structures (`.pdb`) **are** published here, alongside the scripts that read them: the
RFdiffusion target, the PMI benchmark, and the 2,000 AlphaFold Multimer designs that
stage 3 scores.

---

## What the project is

**AETh3R** is a peptide-PROTAC (p-PROTAC) for non-small cell lung cancer, designed to
degrade its targets rather than inhibit them:

- **MDM2** — the negative regulator of p53. Degrading it releases p53-mediated tumour
  suppression. This is the arm the dry lab designed binders for.
- **KRAS** — mutated in ~30% of NSCLC.

The degrader recruits the **VHL** E3 ligase through a literature-derived peptide (pHIF1),
joined to the warhead by a glycine–serine linker.

---

## Hard rules

### 1. These scripts are a record, not a product

They were written for our own pipeline. They are **not** a software tool for others to
use, and nothing here or on the wiki should suggest otherwise. We are not entering the
Software & AI Village and are not nominating for Best Software Tool. READMEs should
describe *what we did*, not offer an API.

### 2. Never change script behaviour

Every number on the Model page came out of these scripts. A refactor that alters output
silently invalidates the wiki.

Permitted: comments, docstrings, variable and function names, formatting, import order,
dead-code removal.

**Not permitted:** any change to logic, thresholds, filters, ordering, default parameters,
output format or output file names.

After editing, re-run the script where inputs still exist and confirm the output is
unchanged. If it cannot be re-run, record that in its README rather than assuming.

### 3. The selection funnel stays as it is

Three stages, two scripts: stages 1 and 2 are `z-score_filter.py` read two ways, stage 3
is `rank_binders.py`. They could be merged into one cleaner script. **Do not.** That the
funnel evolved through three attempts — each written because the previous one missed
something — is the team's main modelling contribution. Merging them destroys the evidence
for the argument the wiki makes.

### 4. Do not mention MM/PBSA

It is not part of the public write-up. Leave it out of READMEs, comments and
documentation. If a script contains MM/PBSA code, leave the code alone but do not
document or surface it.

---

## The dry lab pipeline, in order

### Warhead: de novo MDM2 binders

1. **Interface analysis** — three hydrophobic sub-pockets mapped from 1YCR (the p53
   anchors Phe19, Trp23, Leu26).
2. **Generation** — RFdiffusion (500 backbones) → ProteinMPNN (4 sequences each) →
   AlphaFold Multimer. 2,000 candidates.
3. **Three-stage selection funnel:**
   - *Stage 1* — composite Z-score of AlphaFold confidence metrics → **Z1 (42.3)**
   - *Stage 2* — confidence filter + pharmacophore flags → **Z2 (267.3)**
   - *Stage 3* — interface geometry and chemistry from atomic coordinates → **Z3 (131.2)**

**Why the sequence matters:** Stage 2's top pick, Z2, carries net charge −8 — chemically
implausible for a hydrophobic groove, but confidence metrics rank it highly because they
score predicted geometry and are blind to chemistry. Stage 3 exists to catch that.

### Recruiter

**pHIF1** (LA-Hyp-Y-Hle-P), from Chen et al. Two non-canonical residues —
4-hydroxyproline and homoleucine — which constrain every downstream tool and are the
source of most of the workarounds in the MD and docking scripts.

### Assembly and validation

Four constructs: PMI-PROTAC, PZ1, PZ2, PZ3. Simulated as ternary complexes with MDM2 and
**VBC** (VHL + ElonginB + ElonginC), 100 ns each.

**Headline result:** PZ3 carries the best MDM2 warhead and makes the most contacts with
the target, but forms the *weakest* induced MDM2–VHL interface. PZ1 is the reverse and is
the leading degrader candidate. Target binding and degradation efficiency are separable.

### Molecular dynamics scripts

The MD **setup** scripts are in this repository; the trajectory analysis is in GitLab.
Protocol, which the READMEs should state consistently:

- GROMACS 2022, AMBER99SB-ILDN force field, TIP3P water, 300 K, 1 bar
- **neutralising counter-ions only — no background salt.** `genion` is called with
  `-neutral` and no `-conc`. These are not 0.15 M NaCl simulations.
- cubic box, minimum 1.0 nm solute-to-edge, energy-minimised
- staged equilibration: NVT then NPT with position restraints, then unrestrained production
- nonbonded settings are CHARMM-style and deliberate: rcoulomb = rvdw = 1.2 nm,
  Force-switch from 1.0 nm, no DispCorr. Not the AMBER reference values; used across every
  system so they stay comparable. Report them as written.
- `md.mdp` is **100 ns in both routes**. The published VHL–pHIF1 trajectory is 150 ns —
  that run was continued past what the script specifies.
- pHIF1 as a free peptide is parametrised with GAFF2 via ACPYPE, because homoleucine is
  absent from AMBER99SB-ILDN and pdb2gmx cannot build it. The assembled PROTACs were built
  without homoleucine, leaving 4-hydroxyproline as their only non-canonical residue —
  which the force field already carries — so they needed no ligand parametrisation.
- `run_pipeline.sh` (the general route) passes no `-ff`, so its force field was answered at
  the prompt rather than recorded. Those runs used AMBER99SB-ILDN.

**Check for discrepancies.** These scripts were written at different times. Confirm the
force field, water model, salt, temperature, pressure, box size and run length are
consistent between scripts and with the list above — and flag any that differ rather than
quietly harmonising them. A genuine difference in protocol needs reporting, not hiding.

---

## Pathways (Boolean model)

A separate strand: a Boolean model of apoptotic commitment in **A549** (TP53 wild-type,
KRAS G12S, CDKN2A deleted), asking whether degrading MDM2 and KRAS together kills a cell
that survives either alone.

**Scripts here:** network construction from OmniPath, the curation tooling, the rules file
and its validation tool, value percolation, and the 32-variant screen.

Points worth preserving in the READMEs, because they are what make the work defensible:

- **The OmniPath snapshot is pinned and archived** with URL, timestamp and checksum.
  OmniPath updates continuously, so the network is only reproducible against that
  snapshot. Do not re-query it.
- **Evidence thresholds** — directed, ≥2 independent resources, ≥1 reference, curation
  effort ≥5. Every rejected edge is written to an exclusion file with its reason.
- **One edge was restored by hand against the thresholds**: KRAS→PI3K, supported by a
  single resource. Without it the whole PI3K/AKT arm is disconnected from the perturbed
  node. It is flagged in the network as a manual addition — keep that flag visible.
- **Curation** adjudicated 40 edges by hand across the build and curation steps — 25
  reviewed sign calls in `decisions_final.csv`, 15 corrections in `curate_network.py`.
  Thirty-seven fall into five recognisable error classes (reverse records; signs
  contradicting their own annotations; reciprocal sequestration; text-mining artefacts;
  transcriptional bypasses); the other three are clamped-input drops. PIK3CA/PIK3R1
  merged into one PI3K node; PARP1 removed. All 28 resulting changes are in
  `network_changelog.csv`.
- **Rules**: 29 total — 3 clamped inputs, 15 left at the default convention, 11 structural
  overrides. The validation tool checks each rule back against the edge signs and reports
  **17 edges that no rule uses**, so a deliberate omission stays distinguishable from an
  accident. (Those 17 are unused *edges*, not overrides — easy to conflate.)
- **Asynchronous updating**, one node per step. Fixed points are identical under either
  scheme; only cyclic attractors differ.
- **Value percolation is exact, not sampled.** All 29 nodes fix within six rounds in every
  condition, so each condition has a single provable attractor. Random-start sampling was
  run alongside and agreed in every case.
- **Criteria were frozen before any variant was run.** That ordering is the point — never
  describe them as chosen afterwards.
- **R2 is the PUMA rule** (BBC3): form A is `TP53 OR FOXO3`, form B is `TP53` only. The
  whole synergy argument turns on it, and the data cannot settle it because FOXO3 is
  inactive in every condition the criteria constrain.
- **R1 is reported, not resolved** — the screen needs form B, the binding literature leans
  toward form A. Keep the tension stated.
- **Nothing tests the KRAS arm.** The only proteomics the screen uses is the
  MDM2-inhibition (Nutlin) dataset, so the KRAS-degraded and both-degraded outcomes are
  predictions, not results, and neither contributes a pass criterion. Bortezomib and
  seven other drugs in the panel are held back as the only genuine out-of-sample test.
- Four nodes (BCL2L11, PMAIP1, MCL1, MDM4) went undetected and are unconstrained by any
  criterion.

---

## Canonical naming — use these consistently

| Use | Not |
|---|---|
| AETh3R | AETHER, AETH3R, Aether |
| Z1 (42.3) · Z2 (267.3) · Z3 (131.2) | bare 42.3 / 267.3 / 131.2; "Z4"; 130.3 |
| PZ1 · PZ2 · PZ3 · PMI-PROTAC | — |
| VBC | "the VHL complex" where VBC is meant |
| Dr. Alessio Ciulli | Alessandro Ciulli |

---

## Numbers that appear on the wiki

These act as a checksum. If a cleanup changes any of them, something has gone wrong.

- 2,000 candidates → 1,275 fill 3/3 pockets → 409 pass all gates → **18** meet or exceed PMI
- Z3: 49 contacts, 3/3 pockets, net charge −1, pLDDT 89.5
- Z2: net charge −8 (rejected at the chemistry gate)
- Induced MDM2–VHL contacts: PZ1 ~450–600, PMI-PROTAC ~300–380, PZ2 ~200–300, PZ3 ~130–200
- Pathways: 117 edges after filtering → **29 nodes, 99 edges** after curation; 29 rules
  (3 inputs, 15 default, 11 overrides); 1 fixed point per condition out of 6.7 × 10⁷ states;
  32 variants, 4 satisfy all criteria, **4 of 32** show true synergy (28 kill, but 24 of
  those are KRAS-alone)

### Buried surface area — unsettled

The BSA figures (apolar Å² per residue, and the PMI benchmark) are **not yet settled on
the wiki**. Leave the numbers in the scripts untouched, but do not treat them as fixed
reference values in documentation until the team has resolved how they are presented.

---

## Settled — do not re-open

- **KRAS.** The wiki's framing is correct as written: the dry lab designed MDM2 binders;
  the KRAS strand is an unstapled SOS1 helix that gave a clear negative result (ipTM 0.27,
  no RMSD plateau) and is presented as future work. Do not change this.
- **MDMX.** The MDMX-directed PROTAC is complete and its final outputs exist. Any text
  still saying "under construction for simulation" is outdated, and the page will be
  updated to match.
- **Ternary scope.** We simulated MDM2·PROTAC·**VBC**. We did *not* model the Cullin-2 /
  Rbx1 scaffold. Any claim of the full CRL2^VHL assembly is wrong.
- **Single trajectories.** One run per system, no replicates. Stated as a limitation on
  the wiki; keep it stated.

---

## Layout

```
de novo binders/
  Collab Result/       generation: the 2,000-candidate pool and its inputs
  Filter Pipeline.py/  the three-stage selection funnel, and the figure scripts
  MD simulations/      GROMACS setup for the inhibitor, ternary and VHL systems
  Pathways/            network build, curation, rules, attractors, variant screen
```

Each directory has a README stating what the scripts did, what they took as input, what
they produced, and — for the funnel — what each stage could not catch and which later
stage addressed it. The top-level README records the environment: Python and package
requirements, GROMACS 2022, AMBER99SB-ILDN.
