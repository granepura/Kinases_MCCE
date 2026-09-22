#!/bin/bash
# Shared configuration for the kinase MCCE trials.
# Sourced by setup_trial.sh and derive_apo.sh.  Edit here, not in the scripts.

# Derived from this file's own location, so the repo works wherever it is
# cloned.  Override with KINASES_ROOT=... if you relocate scripts_Kinases_MCCE.
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${KINASES_ROOT:-$(cd "$SCRIPTS/.." && pwd)}"

# The canonical stepC script.  submit_mcce4*.sh points STEPC at THIS path in
# every trial, so all trials share one implementation instead of drifting copies.
PRUNE="$SCRIPTS/prune_kin-inhib_head3.py"

# The canonical stepB scripts.  run_holo's step1-2 job splits the finished
# step2_out.pdb into holo_step2_out.pdb (exact copy) and apo_step2_out.pdb
# (inhibitor deleted); run_apo's step3-4 job checks that pair and points its own
# step2_out.pdb at the apo one.
VARIANTS="$SCRIPTS/make_holo_apo_step2_out.py"
INSTALL="$SCRIPTS/install_apo_step2_out.py"

# Seeds run_apo from run_holo.  setup_trial.sh copies it into each trial as
# 0-prepare_run_apo.py, so a trial records the version it was actually run with.
PREPARE="$SCRIPTS/prepare_run_apo.py"

# Runs xts_corr.py in every structure of all three trees, after step4.  Copied
# into each trial as 1-run_xts_corr.py.
XTSRUN="$SCRIPTS/run_xts_corr.py"

# Figure scripts, copied into each trial as 2- and 3-.  They read
# xts_sum_crg.out and write per-trial PNGs + CSVs; they never touch the runs.
FIG3="$SCRIPTS/plot_sumcrg_inhibitors_xts_Fig3.py"
FIG4A="$SCRIPTS/plot_sumcrg_comparison_xts_Fig4A.py"

# CPUs for step3's PBE solver.  Only the step3-4 scripts benefit -- step1/2 and
# the tiny ligand-only runs are single-threaded anyway.
CPUS_S34=5
CPUS_DEFAULT=1

# Masters at the repo root, shared by every trial (read-only sources, never run
# in place).  They used to live in cof_tpl_GR, which has been removed.
MASTER_KIN_PDB="$ROOT/kin-pdb"            # holo PDBs        (37)
MASTER_COF_PDB="$ROOT/cof-pdb"            # ligand-only PDBs (37)

# One list for all trials: symlinked into each trial rather than copied, so a
# change to the kinase/inhibitor mapping cannot leave trials disagreeing.
MASTER_INHIB_LST="$ROOT/pdb_inhibitor.lst"   # PDB -> inhibitor code -> kinase

# SBATCH/env template, kept in the repo so a fresh clone can scaffold a trial
# without depending on an existing one.  It is Trial01's step1-2 script with the
# job name, seed, CPUS and hooks neutralised; gen_submit patches all of those
# per output script.
MASTER_SUBMIT_TPL="$SCRIPTS/submit_mcce4_template.sh"

# MCCE step commands.  $PYEX/$MCBIN/$EPS/$CPUS/$TMP are expanded inside the
# submit script at run time, so they stay single-quoted here.
STEP1_CMD='$PYEX $MCBIN/step1.py -d $EPS --noter --dry'
STEP2_CMD='$PYEX $MCBIN/step2.py -d $EPS -l 1'
STEP3_CMD='$PYEX $MCBIN/step3.py -d $EPS -s delphi -salt 0.15 --fly -p $CPUS -t $TMP'
STEP4_BASE='$PYEX $MCBIN/step4.py --xts -i 7.4 -n 1'

# Monte Carlo seed.  MCCE defaults to MONTE_SEED=-1 (time-based), which makes a
# trial unreproducible and can collide when two jobs start in the same second.
# An explicit per-trial seed keeps trials independent AND repeatable.
# Set USE_EXPLICIT_SEED=0 to fall back to MCCE's -1 default.
USE_EXPLICIT_SEED=1
SEED_BASE=1000          # Trial01 -> 1001, Trial02 -> 1002, Trial03 -> 1003

# The list of holo step3/step4 products NOT carried into run_apo lives in
# prepare_run_apo.py (EXCLUDE), which is the script that does the seeding.

trial_seed() { echo $(( SEED_BASE + $1 )); }
trial_dir()  { printf "%s/Trial%02d" "$ROOT" "$1"; }
