# Molecular dynamics — building the systems

The GROMACS setup for every system we simulated: the MDM2 inhibitor
complexes, the ternary MDM2·PROTAC·VBC complexes, and VHL with pHIF1.

This directory is **setup only**. The trajectories themselves, and the
scripts that analyse them, are published in the iGEM GitLab repository.

The scripts are written for one HPC cluster, with SLURM directives and
module names baked in — they will need adapting to run anywhere else.

---

## Two routes, and why there are two

Which script a system goes through depends on one thing: whether `pdb2gmx`
can build its topology at all.

### Most systems — `run_pipeline.sh`

Takes a PDB, builds the topology, solvates, adds ions, and writes the five
`.mdp` files and four SLURM stage scripts. This is the route for the MDM2
inhibitor complexes and the ternary complexes.

The assembled PROTACs (PZ1–PZ3, PMI-PROTAC) come through here. They were
built without homoleucine, which leaves 4-hydroxyproline as their only
non-canonical residue — and AMBER99SB-ILDN already carries that one — so
`pdb2gmx` handles them directly.

### VHL with pHIF1 — the three-step route

pHIF1 is **LA-Hyp-Y-Hle-P**, and homoleucine is simply absent from
AMBER99SB-ILDN, so `pdb2gmx` cannot build the peptide as a polymer at all.
It has to be treated as a small molecule instead:

1. **`prep_vhl_peptide.sh`** splits the Boltz-2 prediction into protein and
   peptide, adds hydrogens with Open Babel, and parametrises the peptide
   with GAFF2 via ACPYPE
2. **`combine_topology.py`** merges the AMBER99SB-ILDN protein and the GAFF2
   peptide into one system
3. **`run_pipeline_lig.sh`** takes it from there

4-hydroxyproline alone would not have forced this — GROMACS ships an HYP
entry. It is homoleucine that does, which is why only the free peptide needs
the detour.

---

## The protocol

| | |
|---|---|
| GROMACS | 2022 |
| force field | AMBER99SB-ILDN, plus GAFF2 for the free pHIF1 peptide |
| water | TIP3P |
| ions | neutralising counter-ions only — no background salt |
| box | cubic, at least 1.0 nm from solute to edge |
| minimisation | steepest descent, up to 50,000 steps |
| equilibration | 100 ps NVT then 100 ps NPT, 300 K / 1 bar, with position restraints |
| production | 100 ns, unrestrained, 2 fs steps, h-bonds constrained |
| frames | one every 10 ps — 10,000 per trajectory |
| replicates | none; one run per system |

Production runs in 24 h chunks, which is the cluster's wall-time limit.
`mdrun -maxh 23.5` makes each chunk checkpoint and exit half an hour before
SLURM would kill it, and the next chunk picks up from there. A 100 ns run
takes several chunks and GROMACS stops on its own once it reaches the step
count.

**On run length:** `md.mdp` specifies 100 ns in both routes. The published
VHL–pHIF1 trajectory is 150 ns — that run was continued past what the
script says, so running the script as committed gives you the first 100 ns
of it.

**On ions:** `genion` is called with `-neutral` and no `-conc`. The box
carries only the counter-ions needed to balance the charge. These are not
0.15 M NaCl simulations.

### The nonbonded settings are CHARMM-style on purpose

```
rcoulomb = rvdw = 1.2 nm,  vdw-modifier = Force-switch,  rvdw-switch = 1.0,
no DispCorr
```

These are not the AMBER reference values (1.0 nm plain cutoff with
`DispCorr = EnerPres`), even though the force field is AMBER99SB-ILDN. The
same settings are used for every simulation in the project so the systems
stay comparable with one another. The longer cutoff is more accurate in
direct space and force-switching conserves energy better; dropping
`DispCorr` slightly underestimates density.

The effect on RMSD, RMSF, Rg, contacts and buried surface area is
negligible. Absolute energies would be affected, so these trajectories
support structural comparison and not energetic estimates.

**Report these values as written, not the AMBER defaults.**

---

## If you re-run any of this

- **`run_pipeline.sh` will stop and ask you which force field to use.** It
  passes no `-ff`, so GROMACS prompts. Our runs answered AMBER99SB-ILDN.
  (`prep_vhl_peptide.sh` pins the equivalent choice instead, by piping the
  menu selections in.) Adding `-ff amber99sb-ildn` would let it run
  unattended.
- **Nothing gets submitted automatically.** Both scripts write the four
  SLURM stage scripts and stop. Submitting them, and resubmitting `md.slurm`
  for each production chunk, was done by hand.
- **`run_pipeline.sh` describes `-ignh` but does not pass it**, so existing
  hydrogens are kept rather than rebuilt.
- **`combine_topology.py` edits `topol.top` in place** with no backup, and
  appends the ligand to `[ molecules ]` at end-of-file — which is correct
  only while that is the last section, as it is in normal `pdb2gmx` output.
- **pHIF1 uses Gasteiger charges** (`acpype -c gas`) rather than GAFF2's
  usual AM1-BCC partner. Fine for the structural observables these
  trajectories were analysed for, and the reason no binding free energy
  should be quoted from that system.
