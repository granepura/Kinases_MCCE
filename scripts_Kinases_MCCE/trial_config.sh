#!/bin/bash
# Shared configuration for the kinase MCCE trials.
# Sourced by setup_trial.sh and derive_apo.sh.  Edit here, not in the scripts.

ROOT="/data/home/granepura/5-Kinases/Kinases_MCCE"
SCRIPTS="$ROOT/scripts_Kinases_MCCE"

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

# Masters copied into a trial at setup (read-only sources, never run in place)
MASTER_KIN_PDB="$ROOT/cof_tpl_GR/run_kin/kin-pdb"        # holo PDBs        (37)
MASTER_COF_PDB="$ROOT/cof_tpl_GR/run_cof2/cof-pdb"       # ligand-only PDBs (37)
MASTER_INHIB_LST="$ROOT/cof_tpl_GR/pdb_inhibitor.lst"    # PDB -> inhibitor code
MASTER_SUBMIT_TPL="$ROOT/cof_tpl_GR/run_kin/submit_mcce4_s1s2.sh"  # SBATCH/env template

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
