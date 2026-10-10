#!/bin/bash
# =============================================================================
# run_pipeline.sh
#
# GROMACS 2022, 100 ns MD -- AMBER99SB-ILDN / cubic box / TIP3P /
# neutralising counter-ions only
#
# Protocol based on: doi:10.1186/s43088-026-00765-8 (100 ns instead of 300 ns)
#
# WHAT THIS DID
#   The general-purpose route: any system pdb2gmx can build on its own.
#   That is the MDM2 inhibitor complexes and the ternary
#   MDM2 - PROTAC - VBC complexes, 100 ns each, one run per system with no
#   replicates. Systems containing pHIF1 cannot go through here -- see
#   run_pipeline_lig.sh and prep_vhl_peptide.sh for why.
#
# INPUT
#   one PDB, passed as $1
#
# PRODUCED
#   ions.mdp minim.mdp nvt.mdp npt.mdp md.mdp   run parameters
#   processed.gro boxed.gro solvated.gro solv_ions.gro
#   topol.top
#   em.slurm nvt.slurm npt.slurm md.slurm       stage scripts, NOT submitted
#
# THE FORCE FIELD IS NOT PINNED BY THIS SCRIPT
#   pdb2gmx below is called without -ff, so GROMACS prompts for the force
#   field and it is answered at the prompt. The runs used AMBER99SB-ILDN,
#   matching prep_vhl_peptide.sh, which pins it rather than asking.
#   Passing "-ff amber99sb-ildn" would let this run unattended.
#
# IONIC STRENGTH: NEUTRALISING IONS ONLY
#   genion is called with -neutral and no -conc, so the box carries only
#   the counter-ions needed to neutralise the system. There is no
#   background salt in these simulations.
#
# USAGE (run from the folder containing your PDB):
#   bash run_pipeline.sh PMI-MDM2.pdb
#
# Preprocessing runs here directly (<1 min). The four SLURM stage scripts
# are then WRITTEN BUT NOT SUBMITTED -- despite the note below, there is no
# sbatch call anywhere in this file. Submit them yourself as a dependency
# chain, em -> nvt -> npt -> md.
# =============================================================================
set -euo pipefail

PDB_FILE="${1:?Usage: bash run_pipeline.sh <input.pdb>}"

# ---------------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------------

WATER_MODEL="tip3p"
BOX_TYPE="cubic"      # paper specifies cubic box
BOX_DIST="1.0"        # minimum solute-to-wall distance (nm)
ION_PNAME="NA"
ION_NNAME="CL"        # -neutral picks whichever of the two is needed; for a
                      # net-negative solute that is Na+ only, for a
                      # net-positive one Cl- only. No background salt either
                      # way -- see the header.

GMX_MODULES="gcc/9.4.0-eewq4j6 openmpi/4.1.2-akxtzzl gromacs/2022-47qrtrj"
PARTITION="gpu"
EQ_TIME="02:00:00"    # NVT/NPT are 100 ps each -- finish in minutes
PROD_TIME="24:00:00"  # per production chunk
PROD_CHUNKS=2         # 2 x 24h is a safe ceiling for 100 ns;
                       # GROMACS exits cleanly once nsteps is reached
                       # NOTE: PROD_CHUNKS is never read. Nothing in this
                       # script loops or submits, so the chunking it
                       # describes is done by hand.

if [ ! -f "$PDB_FILE" ]; then
    echo "ERROR: '$PDB_FILE' not found in $(pwd)." >&2
    exit 1
fi

module purge
module load $GMX_MODULES

# ---------------------------------------------------------------------------
# NONBONDED SETTINGS
#   1. rcoulomb = rvdw = 1.2 nm  (not 1.0)
#   2. vdwtype = Cut-off  +  vdw-modifier = Force-switch  +  rvdw-switch = 1.0
#   3. NO DispCorr
#
# These three are CHARMM-style settings, and they are NOT the AMBER
# reference values (1.0 nm plain cutoff with DispCorr = EnerPres) even
# though the force field here is AMBER99SB-ILDN. That is deliberate: the
# same settings are used for every simulation in this project, so the
# systems stay comparable with each other. run_pipeline_lig.sh carries the
# same block and the same reasoning.
#
# The trade-off, as stated there: the longer cutoff is more accurate in
# direct space and force-switching improves energy conservation, while
# omitting DispCorr slightly underestimates density. Structural
# observables -- RMSD, RMSF, Rg, contacts, buried surface area -- are
# affected negligibly; absolute free energies would not be.
#
# Report these values as written, not the AMBER defaults.
# ---------------------------------------------------------------------------
cat > ions.mdp << 'EOF'
; ions.mdp 
integrator        = steep
emtol             = 1000.0
emstep            = 0.01
nsteps            = 1000
nstlist           = 1
cutoff-scheme     = Verlet
coulombtype       = cutoff
rcoulomb          = 1.2
vdwtype           = Cut-off
vdw-modifier      = Force-switch
rvdw-switch       = 1.0
rvdw              = 1.2
pbc               = xyz
EOF

cat > minim.mdp << 'EOF'
; minim.mdp - 50,000 steps steepest descent
integrator        = steep
emtol             = 1000.0
emstep            = 0.01
nsteps            = 50000
nstlist           = 1
cutoff-scheme     = Verlet
coulombtype       = PME
rcoulomb          = 1.2
vdwtype           = Cut-off
vdw-modifier      = Force-switch
rvdw-switch       = 1.0
rvdw              = 1.2
pbc               = xyz
EOF

cat > nvt.mdp << 'EOF'
; nvt.mdp - 100 ps NVT at 300 K
define                = -DPOSRES
integrator            = md
nsteps                = 50000
dt                    = 0.002
nstxout-compressed    = 500
nstenergy             = 500
nstlog                = 500
cutoff-scheme         = Verlet
coulombtype           = PME
rcoulomb              = 1.2
vdwtype               = Cut-off
vdw-modifier          = Force-switch
rvdw-switch           = 1.0
rvdw                  = 1.2
tcoupl                = V-rescale
tc-grps               = Protein Non-Protein
tau_t                 = 0.1 0.1
ref_t                 = 300 300
constraints           = h-bonds
continuation          = no
gen_vel               = yes
gen_temp              = 300
gen_seed              = -1
EOF

cat > npt.mdp << 'EOF'
; npt.mdp - 100 ps NPT at 300 K / 1 bar
define                = -DPOSRES
integrator            = md
nsteps                = 50000
dt                    = 0.002
cutoff-scheme         = Verlet
coulombtype           = PME
rcoulomb              = 1.2
vdwtype               = Cut-off
vdw-modifier          = Force-switch
rvdw-switch           = 1.0
rvdw                  = 1.2
tcoupl                = V-rescale
tc-grps               = Protein Non-Protein
tau_t                 = 0.1 0.1
ref_t                 = 300 300
pcoupl                = Parrinello-Rahman
pcoupltype            = isotropic
tau_p                 = 2.0
ref_p                 = 1.0
compressibility       = 4.5e-5
constraints           = h-bonds
continuation          = yes
refcoord_scaling      = com
EOF

cat > md.mdp << 'EOF'
; md.mdp - 100 ns production
integrator            = md
dt                    = 0.002
nsteps                = 50000000   ; 100 ns
cutoff-scheme         = Verlet
coulombtype           = PME
rcoulomb              = 1.2
vdwtype               = Cut-off
vdw-modifier          = Force-switch
rvdw-switch           = 1.0
rvdw                  = 1.2
tcoupl                = V-rescale
tc-grps               = Protein Non-Protein
tau_t                 = 0.1 0.1
ref_t                 = 300 300
pcoupl                = Parrinello-Rahman
pcoupltype            = isotropic
tau_p                 = 2.0
ref_p                 = 1.0
compressibility       = 4.5e-5
constraints           = h-bonds
continuation          = yes
nstxout-compressed    = 5000   ; 1 frame per 10 ps -> 10,000 frames total
nstenergy             = 5000
nstlog                = 5000
EOF

echo "[1/3] .mdp files written (rcoulomb=rvdw=1.2 nm, Force-switch vdW, no DispCorr)"

# ---------------------------------------------------------------------------
# PREPROCESSING (runs directly here, no GPU needed, takes <1 min)
# ---------------------------------------------------------------------------

# TWO THINGS TO KNOW ABOUT THIS CALL:
#
#   1. NO -ff. GROMACS prompts for the force field and waits for an answer,
#      so this does not run unattended and does not record what was
#      chosen. The runs used AMBER99SB-ILDN. prep_vhl_peptide.sh pins the
#      equivalent choice by piping "6\n1" into pdb2gmx.
#
#   2. NO -ignh either, despite the note that follows -- it describes a
#      flag this call does not pass. Existing hydrogens are kept rather
#      than rebuilt, the opposite of what that note implies.
#
# -ignh strips any existing H atoms and lets GROMACS re-add them correctly
# (important for AlphaFold outputs or PDB structures with non-standard H names)
gmx_mpi pdb2gmx \
    -f "$PDB_FILE" \
    -o processed.gro \
    -p topol.top \
    -water "$WATER_MODEL" \

gmx_mpi editconf \
    -f processed.gro \
    -o boxed.gro \
    -c \
    -d "$BOX_DIST" \
    -bt "$BOX_TYPE"

gmx_mpi solvate \
    -cp boxed.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p topol.top

gmx_mpi grompp \
    -f ions.mdp \
    -c solvated.gro \
    -p topol.top \
    -o ions.tpr

# SOL: select the solvent group for genion to replace with ions.
# -neutral and no -conc: counter-ions to neutralise the net charge, and
# nothing else. There is no 0.15 M background salt in these systems.
echo "SOL" | gmx_mpi genion \
    -s ions.tpr \
    -o solv_ions.gro \
    -p topol.top \
    -pname "$ION_PNAME" \
    -nname "$ION_NNAME" \
    -neutral

echo "[2/3] Preprocessing complete -> solv_ions.gro"

# ---------------------------------------------------------------------------
# WRITE SLURM STAGE SCRIPTS
# ---------------------------------------------------------------------------
cat > em.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_em
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --time=$EQ_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

gmx_mpi grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx_mpi mdrun -deffnm em -nb gpu -pme cpu -pin on
# Note: -update gpu is intentionally omitted here; the GPU leapfrog update
# only applies to the md integrator, not steepest-descent minimization.
EOF

cat > nvt.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_nvt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --time=$EQ_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx_mpi mdrun -deffnm nvt -nb gpu -pme gpu -update gpu -pin on
EOF

cat > npt.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_npt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus=1
#SBATCH --time=$EQ_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx_mpi mdrun -deffnm npt -nb gpu -pme gpu -update gpu -pin on
EOF

cat > md.slurm << EOF
#!/bin/bash
#SBATCH --partition=$PARTITION
#SBATCH --job-name=gmx_md
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1
#SBATCH --time=$PROD_TIME
#SBATCH --gpu-bind=closest

module purge
module load $GMX_MODULES
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK
cd \$SLURM_SUBMIT_DIR

# grompp only on the first chunk; subsequent chunks reuse md.tpr
if [ ! -f md.tpr ]; then
    gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
fi

# -maxh: saves a checkpoint and exits cleanly 30 min before SLURM's hard kill
# so the next chunk can pick up exactly where this one left off
if [ -f md.cpt ]; then
    gmx_mpi mdrun -deffnm md -cpi md.cpt -nb gpu -pme gpu -update gpu -pin on -maxh 23.5
else
    gmx_mpi mdrun -deffnm md -nb gpu -pme gpu -update gpu -pin on -maxh 23.5
fi
EOF

chmod +x em.slurm nvt.slurm npt.slurm md.slurm
echo "[3/3] SLURM stage scripts written: em.slurm nvt.slurm npt.slurm md.slurm"

# NOT SUBMITTED. The header describes a dependency chain, but this script
# ends here -- there is no sbatch call. Submit by hand, e.g.
#   J=$(sbatch --parsable em.slurm)
#   J=$(sbatch --parsable --dependency=afterok:$J nvt.slurm)
#   J=$(sbatch --parsable --dependency=afterok:$J npt.slurm)
#   sbatch --dependency=afterok:$J md.slurm
