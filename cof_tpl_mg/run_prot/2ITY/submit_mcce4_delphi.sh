#!/bin/bash
#SBATCH --job-name=run_prot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12G                 # Adjust memory if needed

# ==============================================================================
# Script Name   : submit_mcce4.sh
# Purpose       : Automate and control the full MCCE4 simulation pipeline including optional custom preprocessing steps.
#
# Description   :
#   This script manages the sequential execution of MCCE4 simulation steps (1 to 4), with optional hooks (stepM, stepA, stepB, stepC)
#   that allow the user to insert custom membrane generation and intermediate processing scripts.
#   It records the timing and success/failure of each step in a detailed log file (`mcce_timing.log`).
#   The script supports flexible control through flags to enable/disable specific MCCE steps or custom steps.
#
# Main Features :
#   - Step 1: Convert standard PDB to MCCE-compatible format
#   - Step 2: Generate rotamers
#   - Step 3: Perform energy calculations
#   - Step 4: Run Monte Carlo sampling
#   - Step clean: Clean temporary pbe_data
#   - Optional StepM: Add membrane-specific conformers (e.g., using IPECE)
#   - Optional StepA/B/C: Insert custom preprocessing scripts between core steps
#   - Intelligent skip logic and output checking to prevent redundant work
#   - Runtime logging with timestamps for each phase
#
# Input Requirements :
#   - input_pdb       : A protein PDB file named `prot.pdb`
#   - MCCE_HOME       : Path to MCCE4 installation directory
#   - USER_PARAM      : Directory with user-defined MCCE parameters (optional)
#   - EXTRA           : Custom extra.tpl file (optional)
#   - TMP             : Path to store temporary pbe_data (default: /tmp)
#   - CPUS            : Multiprocessing for step3
#   - Optional scripts for stepM/A/B/C must exist and be executable if enabled
#
# Output Files :
#   - step1_out.pdb, step2_out.pdb, head3.lst, pK.out
#   - Timing report: mcce_timing.log
#   - Logs for each step: step1.log, step2.log, etc.
#
# Usage :
#   Set control flags (`step1`, `step2`, etc.) to "t" or "f" to enable/disable each step.
#   Set paths to optional scripts as needed.
#   Submit this script to a SLURM cluster or run locally if sbatch is not used.
#
# Author        : Gehan A. Ranepura
# Date Created  : July 2025
# ==============================================================================

#=============================================================================
#-----------------------------------------------------------------------------
# Input and Output:
input_pdb="prot_new.pdb"    # (INPDB)

# Set MCCE4 Parameters
MCCE_HOME="/home/granepura/MCCE4"
USER_PARAM="./user_param"
EXTRA="./user_param/extra.tpl"
TMP="/tmp"
CPUS=1

# Step control flags
step1="t"               # STEP1: pre-run, pdb-> mcce pdb  (DO_PREMCCE)
step2="t"               # STEP2: make rotamers            (DO_ROTAMERS)
step3="t"               # STEP3: Energy calculations      (DO_ENERGY)
step4="t"               # STEP4: Monte Carlo Sampling     (DO_MONTE)
step_clean="t"          # Clean PBE data                  (BACKUP CLEAN) : Set to f if step3 --debug option is used

# Optional step controls
stepM="f"               # Generate Partial Membranes                    : If true, user MUST satisfy condidtions of stepM.sh, which can be be obtained on MCCE4/inhouse/stepM.sh
stepA="f"               # Run a custom script between step1 and step2   : If true, user MUST satisfy condidtions of their custom script
stepB="f"               # Run a custom script between step2 and step3   : If true, user MUST satisfy condidtions of their custom script
stepC="t"               # Run a custom script between step3 and step4   : If true, user MUST satisfy condidtions of their custom script

# MCCE Simulation
STEP1="step1.py -d 4 --noter --dry"
STEP2="step2.py -d 4 -l 1"
STEP3="step3.py -d 4 -s delphi -salt 0.15 --fly -p $CPUS -t \$TMP"
STEP4="step4.py --xts --ms -i 7.4 -n 1"

# Optional MCCE script locations
STEPM="/path/to/stepM.sh"         # Optional StepM: Bash script
STEPA="/path/to/stepA_script.py"  # Optional StepA: Python script to run between step1 and step2.
STEPB="/path/to/stepB_script.py"  # Optional StepB: Python script to run between step2 and step3.
#STEPC="/home/gunner/source/update_vdw0.py"   # Optional StepC: Python script to run between step3 and step4.
STEPC="/home/granepura/5-kinases/vdw0_update2zero.py"  # Optional StepC: Python script to run between step3 and step4.

# NOTE: User is responsible to precheck if custom scripts work properly and efficiently
# NO USER INPUT NECCESARY BELOW THIS LINE
#------------------------------------------------------------------------------
#==============================================================================

# Check MCCE_HOME exists before PATH export
if [[ ! -d "$MCCE_HOME/bin" ]]; then
    echo -e "\033[0;31m[ERROR]\033[0m MCCE_HOME is not set correctly or $MCCE_HOME/bin does not exist."
    echo "Please check your MCCE_HOME path in submit_mcce4.sh."
    exit 1
fi

# Export environment for downstream script
export input_pdb MCCE_HOME USER_PARAM EXTRA TMP
export step1 step2 step3 step4 step_clean
export stepM stepA stepB stepC
export STEP1 STEP2 STEP3 STEP4
export STEPM STEPA STEPB STEPC

# Inititiate MCCE_HOME PATH and call driver_mcce4.sh
export PATH="$MCCE_HOME/bin:$PATH"  # Add MCCE_HOME/bin to PATH for mcce step executables
driver_mcce4.sh

