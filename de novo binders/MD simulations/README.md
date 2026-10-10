# Molecular dynamics — setup

How the MD systems were built and run. This is a record of our protocol, not
a general-purpose pipeline; the scripts are written for one HPC cluster with
SLURM and the module names baked in.

Everything here is **setup only**. The trajectory analysis is not in this
repository — see *What is missing* at the end.

---

## Two routes

Which script a system goes through depends on one thing: whether `pdb2gmx`
can build its topology.

### `run_pipeline.sh` — the general route

For anything `pdb2gmx` handles on its own: the MDM2 inhibitor complexes and
the ternary MDM2·PROTAC·**VBC** complexes. Takes a PDB, builds the topology,
solvates, ionises, writes the five `.mdp` files and four SLURM stage scripts.

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

4-hydroxyproline on its own would not have forced this — GROMACS ships an HYP
entry for AMBER99SB-ILDN. It is homoleucine that does.

---

## The protocol as the scripts actually implement it

| | |
|---|---|
| GROMACS | 2022 (`gromacs/2022-47qrtrj`) |
| force field | AMBER99SB-ILDN (protein) + GAFF2 via ACPYPE (pHIF1 only) |
| water | TIP3P, `spc216.gro` starting configuration |
| ions | **neutralising counter-ions only — no background salt** |
| box | cubic, minimum 1.0 nm solute-to-edge |
| minimisation | steepest descent, 50,000 steps, `emtol` 1000 kJ/mol/nm |
| NVT | 100 ps, 300 K, V-rescale (`tau_t` 0.1), position restraints |
| NPT | 100 ps, 300 K / 1 bar, Parrinello-Rahman (`tau_p` 2.0), restraints |
| production | 100 ns, unrestrained, `dt` 2 fs, h-bonds constrained |
| frames | one per 10 ps → 10,000 per trajectory |
| replicates | **none — one run per system** |

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

---

## Three places the scripts and the wiki disagree

All three were put to the team in October 2026 and answered. **The wiki and
CLAUDE.md need correcting, not the scripts.**

### 1. There is no 0.15 M NaCl

Both scripts call `genion` with `-neutral` and no `-conc`, so the box carries
only the counter-ions needed to neutralise the system. The wiki's "0.15 M
NaCl" does not describe these simulations. Confirmed: the scripts are right.

### 2. VHL–pHIF1 is 150 ns on the wiki, 100 ns in the script

`run_pipeline_lig.sh`'s `md.mdp` is `nsteps = 50000000` at `dt = 0.002` —
exactly 100 ns, and labelled as such. The published 150 ns trajectory was
obtained by continuing the run past what the file specifies. **Running this
script as committed reproduces the first 100 ns, not the published
trajectory**, and the extension step is not recorded anywhere here.

### 3. The general route never pinned its force field

`run_pipeline.sh` calls `pdb2gmx` without `-ff`, so GROMACS prompted and
someone answered by hand. Those runs used AMBER99SB-ILDN (confirmed), but the
script does not record it and cannot run unattended.

An earlier copy of that header claimed CHARMM36, in the title and in four
generated `.mdp` comments. That was wrong and has been corrected throughout.
The sibling copy in `igem-scripts-no-claude` said only "FF", claiming
nothing.

Adding `-ff amber99sb-ildn` would fix the reproducibility gap. Not done here,
because it changes what the script does — recommended for a future revision.

---

## Defects left in place

Flagged rather than fixed; each would change behaviour.

- **`run_pipeline_lig.sh`'s `md.slurm` has a wall-time mismatch.**
  `--time=10:00:00` but `mdrun -maxh 23.5`, with a comment claiming it exits
  cleanly 30 minutes before SLURM's kill. SLURM kills at 10 h, ~14 h before
  `-maxh` would fire, so the clean-exit logic never runs. Periodic
  checkpoints meant runs still completed, just over more chunks than
  intended. Either `--time` should be `24:00:00` (matching the unused
  `PROD_TIME`) or `-maxh` should be ~9.5.
- **Neither script submits anything.** Both headers describe a SLURM
  dependency chain; neither contains an `sbatch` call. The four stage scripts
  are written and `chmod +x`'d, and submitting them was a manual step.
  `PROD_CHUNKS` in `run_pipeline.sh` and `PROD_TIME` in `run_pipeline_lig.sh`
  are set and never read.
- **`run_pipeline.sh` documents `-ignh` but does not pass it.** The comment
  explains why stripping and rebuilding hydrogens matters for AlphaFold
  output; the `pdb2gmx` call below it omits the flag, so existing hydrogens
  are kept.
- **`combine_topology.py` edits `topol.top` in place with no backup**, and
  appends the `[ molecules ]` entry at end-of-file rather than inside the
  section. Correct only while `[ molecules ]` is last, which is the normal
  `pdb2gmx` layout and was true here.
- **pHIF1 uses Gasteiger charges** (`acpype -c gas`), not GAFF2's usual
  AM1-BCC. Defensible for structural observables, and the reason no binding
  free energy should be quoted from that system — a stated limitation, not an
  oversight.

### An unresolved question about the PROTACs

CLAUDE.md says the assembled PROTACs (PZ1–PZ3, PMI-PROTAC) "needed no special
ligand parametrisation" because hydroxyproline is in the force field. By
`prep_vhl_peptide.sh`'s own argument that reasoning does not follow: it is
homoleucine, not hydroxyproline, that forces the ligand route, and the
assembled PROTACs contain pHIF1 — homoleucine included. Either those
constructs were built another way, or by a route not recorded here. **Needs
settling before the protocol is published.**

---

## What we could not verify

**Nothing in this directory was re-run.** The documentation container had no
GROMACS, ACPYPE, Open Babel or SLURM, and no cluster access. All four scripts
were checked for syntax only — `bash -n` on the three shell scripts,
`py_compile` on the Python one, plus a line-by-line confirmation that no
executable line changed:

- `run_pipeline_lig.sh` and `combine_topology.py` — **zero** non-comment
  changes (the latter verified by AST comparison with docstrings excluded)
- `prep_vhl_peptide.sh` — one `echo` corrected, which pointed readers at
  `run_pipeline.sh` instead of `run_pipeline_lig.sh`
- `run_pipeline.sh` — four `.mdp` comment lines and one `echo` corrected, all
  for the CHARMM36 error. **No `.mdp` parameter changed.**

## What is missing

- **The trajectory-analysis scripts.** CLAUDE.md says the analysis was done
  with native GROMACS tools and MDAnalysis, but nothing in this repository
  produces the induced MDM2–VHL contact counts (PZ1 ~450–600, PMI-PROTAC
  ~300–380, PZ2 ~200–300, PZ3 ~130–200), the RMSD plateaus, or the buried
  surface area figures. Those numbers cannot currently be traced to code.
- **The input structures** for the production runs, and the step that
  extended VHL–pHIF1 to 150 ns.
