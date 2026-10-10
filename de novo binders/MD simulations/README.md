# Molecular dynamics — setup

How the MD systems were built and run. This is a record of our protocol, not
a general-purpose pipeline; the scripts are written for one HPC cluster, with
SLURM and the module names baked in.

Everything here is **setup only**. Trajectory analysis, and the trajectories
themselves, are published separately in the iGEM GitLab repository.

---

## Two routes

Which script a system goes through depends on one thing: whether `pdb2gmx`
can build its topology.

### `run_pipeline.sh` — the general route

For anything `pdb2gmx` handles on its own: the MDM2 inhibitor complexes and
the ternary MDM2·PROTAC·**VBC** complexes. Takes a PDB, builds the topology,
solvates, ionises, writes the five `.mdp` files and four SLURM stage scripts.

The assembled PROTACs (PZ1–PZ3, PMI-PROTAC) come through here. They were
built without homoleucine, leaving 4-hydroxyproline as their only
non-canonical residue, and AMBER99SB-ILDN already carries that one as a
heteroatom residue — so they need no ligand parametrisation.

### `prep_vhl_peptide.sh` → `combine_topology.py` → `run_pipeline_lig.sh`

For the VHL–pHIF1 system. pHIF1 is **LA-Hyp-Y-Hle-P**, and homoleucine is
absent from AMBER99SB-ILDN entirely, so `pdb2gmx` cannot build the peptide as
a polymer. Instead:

1. `prep_vhl_peptide.sh` splits the Boltz-2 prediction on ATOM/HETATM, adds
   hydrogens with Open Babel, and parametrises the peptide as a single small
   molecule with GAFF2 via ACPYPE
2. `combine_topology.py` merges the AMBER99SB-ILDN protein and the GAFF2
   peptide into one `complex.gro` / `topol.top`
3. `run_pipeline_lig.sh` picks up from there

4-hydroxyproline on its own would not have forced this. It is homoleucine
that does, which is why only pHIF1 as a free peptide needs this route.

---

## The protocol

| | |
|---|---|
| GROMACS | 2022 (`gromacs/2022-47qrtrj`) |
| force field | AMBER99SB-ILDN (protein) + GAFF2 via ACPYPE (free pHIF1 only) |
| water | TIP3P, `spc216.gro` starting configuration |
| ions | neutralising counter-ions only — no background salt |
| box | cubic, minimum 1.0 nm solute-to-edge |
| minimisation | steepest descent, 50,000 steps, `emtol` 1000 kJ/mol/nm |
| NVT | 100 ps, 300 K, V-rescale (`tau_t` 0.1), position restraints |
| NPT | 100 ps, 300 K / 1 bar, Parrinello-Rahman (`tau_p` 2.0), restraints |
| production | 100 ns, unrestrained, `dt` 2 fs, h-bonds constrained |
| frames | one per 10 ps → 10,000 per trajectory |
| replicates | none — one run per system |

Production is submitted in 24 h chunks, the cluster's wall-time limit, with
`mdrun -maxh 23.5` so each chunk checkpoints and exits cleanly half an hour
before SLURM would kill it. The next chunk restarts from `md.cpt`. A 100 ns
run takes several chunks; GROMACS exits on its own once `nsteps` is reached.

### Ionic strength

Both scripts call `genion` with `-neutral` and no `-conc`: the box carries
only the counter-ions needed to neutralise the system, and no background
salt. These are not 0.15 M NaCl simulations.

### Run length

`md.mdp` is 100 ns in both routes. The published VHL–pHIF1 trajectory is
150 ns — that run was continued past what the script specifies, so the script
alone gives the first 100 ns of it. Every other system is 100 ns as written.

### Nonbonded settings are CHARMM-style, deliberately

```
rcoulomb = rvdw = 1.2 nm,  vdw-modifier = Force-switch,  rvdw-switch = 1.0,
no DispCorr
```

These are **not** the AMBER reference values (1.0 nm plain cutoff with
`DispCorr = EnerPres`), even though the force field is AMBER99SB-ILDN. The
same settings are used for every simulation in the project so the systems
stay comparable with each other. The longer cutoff is more accurate in direct
space and force-switching improves energy conservation; omitting `DispCorr`
slightly underestimates density. Structural observables — RMSD, RMSF, Rg,
contacts, buried surface area — are affected negligibly; absolute free
energies would not be.

Report these values as written, not the AMBER defaults.

### The general route does not pin its force field

`run_pipeline.sh` calls `pdb2gmx` without `-ff`, so GROMACS prompts and the
choice is made at the terminal. Those runs used AMBER99SB-ILDN.
`prep_vhl_peptide.sh` pins the equivalent choice by piping `6\n1` into
`pdb2gmx` instead. Passing `-ff amber99sb-ildn` would let the general route
run unattended.

---

## Things to know before changing anything here

- **Neither script submits the SLURM chain it describes.** The four stage
  scripts are written and `chmod +x`'d; submitting them, and resubmitting
  `md.slurm` for each production chunk, was done by hand.
- **`run_pipeline.sh` documents `-ignh` but does not pass it.** The comment
  explains why stripping and rebuilding hydrogens matters for AlphaFold
  output; the `pdb2gmx` call below it omits the flag, so existing hydrogens
  are kept.
- **`combine_topology.py` edits `topol.top` in place with no backup**, and
  appends the `[ molecules ]` entry at end-of-file rather than inside the
  section. Correct only while `[ molecules ]` is last, which is the normal
  `pdb2gmx` layout and is true here.
- **pHIF1 uses Gasteiger charges** (`acpype -c gas`), not GAFF2's usual
  AM1-BCC. Adequate for the structural observables these trajectories were
  analysed for, and the reason no binding free energy should be quoted from
  that system.
- **`PROD_CHUNKS` in `run_pipeline.sh` is set and never read.** Chunking is
  manual.
