#!/bin/bash
#
# setup_trial.sh -- scaffold one kinase MCCE trial.
#
#   TrialNN/
#   ├── run_holo/   kin-pdb/  + submit_mcce4_s1s2.sh, submit_mcce4_s3s4.sh
#   ├── run_apo/    kin-pdb/  + submit_mcce4_s3s4.sh
#   │                          (seeded from run_holo by 0-prepare_run_apo.py;
#   │                           stepB then links step2_out.pdb -> apo_step2_out.pdb)
#   ├── run_inhib/  cof-pdb/  + submit_mcce4.sh        (full step1-4)
#   ├── 0-prepare_run_apo.py  seeds run_apo from run_holo
#   ├── 1-run_xts_corr.py     entropy-corrects all three trees after step4
#   └── RUNBOOK.md  the order to run things + sha256 of the shared scripts
#
# Per-structure directories are NOT created here -- pro_batch builds them from
# the kin-pdb / cof-pdb folders when you launch.
#
# Usage:  ./setup_trial.sh 1 [--force]
#
set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/trial_config.sh"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; RED='\033[0;31m'; RESET='\033[0m'

N="${1:-}"; FORCE=0; RUNBOOK_ONLY=0
case "${2:-}" in
    --force)   FORCE=1 ;;
    --runbook) RUNBOOK_ONLY=1 ;;   # rewrite RUNBOOK.md only, touch nothing else
    "")        ;;
    *) echo "Unknown option: $2"; exit 1 ;;
esac
[[ "$N" =~ ^[0-9]+$ ]] || {
    echo "Usage: $0 <trial-number> [--force | --runbook]"
    echo "  --force    (re)generate scripts and PDB folders"
    echo "  --runbook  rewrite RUNBOOK.md only, leaving submit scripts untouched"
    exit 1
}

TRIAL=$(trial_dir "$N")
TAG=$(printf "T%02d" "$N")
SEED=$(trial_seed "$N")

echo -e "${CYAN}=== setup_trial.sh  ->  $TRIAL ===${RESET}"

for src in "$MASTER_KIN_PDB" "$MASTER_COF_PDB" "$MASTER_INHIB_LST" "$MASTER_SUBMIT_TPL" \
           "$PRUNE" "$VARIANTS" "$INSTALL" "$PREPARE" "$XTSRUN" "$FIG3" "$FIG4A"; do
    [[ -e "$src" ]] || { echo -e "${RED}[FATAL] missing master: $src${RESET}"; exit 1; }
done

# ---------------------------------------------------------------- runbook
# Every trial runs the same seven steps, so this is the same document each time
# apart from the trial number, the seed and the script fingerprints.  The why
# lives in CLAUDE.md; this file is the order of operations plus the provenance
# of what was actually run.
write_runbook() {
{
    echo "# Trial$(printf '%02d' "$N") runbook"
    echo
    echo "Written: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "MONTE_SEED: $([[ $USE_EXPLICIT_SEED -eq 1 ]] && echo "$SEED (explicit)" || echo '-1 (time-based)')"
    cat <<BODY

Every trial is built and run the same way -- only the seed and the job names
differ.  The numbered scripts sit here at the trial root and work the rest out
for themselves; the pro_batch launches are run from inside each run directory.
Invariants and the reasoning behind them are in ../CLAUDE.md.

## Order of operations

Job names match the #SBATCH --job-name in each submit script so a trial's jobs are distinguishable in squeue.
--skip-prerun is used throughout: pro_batch's pre-run check is not needed here,
the PDBs are already curated.

To start a tree from scratch, clear it first (this deletes all results in it):
       rm -rf 1* 2* 3* 4* 5* meta_bench pro_batch_* book.txt

1. inhib, steps 1-4   (shortest; independent of holo/apo, so a good first check)
       cd run_inhib
       pro_batch cof-pdb -custom submit_mcce4.sh -job-name ${TAG}_inhib --skip-prerun
       cat */mcce_timing.log | grep STEP4 | wc -l        # 37 when done

2. holo, steps 1-2
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s1s2.sh -job-name ${TAG}_holo_s1s2 --skip-prerun
       cat */mcce_timing.log | grep STEP2 | wc -l        # 37 when done
   stepB = make_holo_apo_step2_out.py: splits the finished step2_out.pdb into
   holo_step2_out.pdb (exact copy) and apo_step2_out.pdb (inhibitor deleted).

3. seed run_apo from run_holo
       ./0-prepare_run_apo.py                 # --dry-run | --keep | 1XKK 2ITZ
   Copies each whole run_holo/<PDBID>, skipping holo's step2_out.pdb and any
   step3/4 products, then links step2_out.pdb -> apo_step2_out.pdb.
   Needs only step 2, so it can run while holo's step3/4 is still going.

4. holo, steps 3-4
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name ${TAG}_holo_s3s4 --skip-prerun
       cat */mcce_timing.log | grep STEP4 | wc -l        # 37 when done

5. apo, steps 3-4
       cd run_apo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name ${TAG}_apo_s3s4 --skip-prerun
       cat */mcce_timing.log | grep STEP4 | wc -l        # 37 when done
   stepB = install_apo_step2_out.py: re-checks the step2 pair and the link.

   Steps 4 and 5 are independent of each other; run them concurrently.
   Status:  pro_batch --check -job-name <name>     (r pending, c done, e error)

6. entropy correction -- step4 does NOT do this
       ./1-run_xts_corr.py                    # --dry-run | -t run_inhib | --force
   Runs xts_corr.py in all three trees, producing xts_sum_crg.out (plus
   xts_fort.38, entropy_correction.log).  Every figure reads xts_sum_crg.out, so
   this must be done in all three trees or corrected numbers would be compared
   against uncorrected ones.  It prints one line per structure and lists the
   cause of any failure.

7. figures
       ./2-plot_sumcrg_inhibitors_xts_Fig3.py    # inhibitor: bound vs in solution
       ./3-plot_sumcrg_comparison_xts_Fig4A.py   # per residue: holo vs apo
   Add --title to draw titles on the PNGs.  Both write PNGs and a CSV carrying
   the trial name, into plots_Fig3_inhibitors_xts/ and
   plots_Fig4A_holo_vs_apo_xts/.  They only read the runs; nothing is modified.

## Checking one structure

       cat run_holo/1XKK/mcce_timing.log     per-step wall time, success/failure
       cat run_holo/1XKK/stepB.log           the holo/apo step2 split
       cat run_apo/1XKK/stepB.log            the pair check + the step2_out link
       cat run_holo/1XKK/stepC.log           the head3.lst edit
       ls  run_holo/1XKK/pK.out              exists => step4 finished

BODY
    echo '## Provenance -- sha256 of the scripts this trial was set up with'
    echo '```'
    sha256sum "$TRIAL"/[0-3]-*.py "$VARIANTS" "$INSTALL" "$PRUNE" \
              "$SCRIPT_DIR/trial_config.sh" 2>/dev/null \
        | sed "s|$ROOT/||"
    echo '```'
    echo
    echo "Re-check these before comparing trials: the canonical scripts in"
    echo "scripts_Kinases_MCCE/ change over time, and a trial's copies are the"
    echo "record of what it actually ran."
} > "$TRIAL/RUNBOOK.md"
}

if [[ $RUNBOOK_ONLY -eq 1 ]]; then
    [[ -d "$TRIAL" ]] || { echo -e "${RED}[FATAL] $TRIAL does not exist${RESET}"; exit 1; }
    write_runbook
    echo -e "${GREEN}[GEN]${RESET} RUNBOOK.md  (only; nothing else touched)"
    exit 0
fi

if [[ -d "$TRIAL" && $FORCE -eq 0 ]]; then
    echo -e "${RED}[FATAL] $TRIAL already exists. Re-run with --force to overwrite its scripts.${RESET}"
    exit 1
fi

# ---------------------------------------------------------------- submit files
# Rewrites the proven submit_mcce4 template: only the job name, the step flags,
# centering, the stepB/stepC hooks and the step4 seed change.  All the apptainer/
# PATH setup below the "NO USER INPUT" line is inherited untouched.
# The hook lines are padded to one column so STEPB and STEPC line up whichever
# script each points at.
# The hook lines are padded to one column so STEPB and STEPC line up whichever
# script each points at.
HOOKPAD=0
for h in "$PRUNE" "$VARIANTS" "$INSTALL"; do (( ${#h} > HOOKPAD )) && HOOKPAD=${#h}; done

hook_line() {   # hook_line STEPB|STEPC <script path>
    local name=$1 path=$2
    printf '%s="%s"%*s  # Central %s: see scripts_Kinases_MCCE/' \
        "$name" "$path" $(( HOOKPAD - ${#path} )) "" "${name/STEP/Step}"
}

# stepb/stepc are the script each hook runs, or "-" for "hook off".  A hook that
# is off keeps the template's placeholder path, so an unused STEPB never points
# at a real script that someone could switch on by accident.
gen_submit() {
    local out=$1 job=$2 s1=$3 s2=$4 s3=$5 s4=$6 center=$7 stepb=$8 stepc=$9
    local cpus=${10:-$CPUS_DEFAULT}
    # SLURM's stdout file follows the script's own name, so an s3s4 job cannot
    # overwrite the s1s2 log it inherited from the template.
    local logname="$(basename "$out" .sh).log"
    local fb="f" fc="f"
    [[ $stepb != "-" ]] && fb="t"
    [[ $stepc != "-" ]] && fc="t"

    local -a sedargs=(
        -e "s|^#SBATCH --job-name=.*|#SBATCH --job-name=${job}|"
        -e "s|^#SBATCH -o .*|#SBATCH -o ${logname}|"
        -e "s|^CPUS=[0-9]*|CPUS=${cpus}|"
        -e "s|^step1=\"[tf]\"|step1=\"${s1}\"|"
        -e "s|^step2=\"[tf]\"|step2=\"${s2}\"|"
        -e "s|^step3=\"[tf]\"|step3=\"${s3}\"|"
        -e "s|^step4=\"[tf]\"|step4=\"${s4}\"|"
        -e "s|^center=\"[tf]\"|center=\"${center}\"|"
        -e "s|^stepB=\"[tf]\"|stepB=\"${fb}\"|"
        -e "s|^stepC=\"[tf]\"|stepC=\"${fc}\"|"
    )
    # Always rewrite both hook lines: a real path when the hook runs, the
    # template's placeholder when it does not.  That way an unused STEPB never
    # points at a real script, whatever the template happened to contain.
    if [[ $fb == "t" ]]; then
        sedargs+=( -e "s|^STEPB=.*|$(hook_line STEPB "$stepb")|" )
    else
        sedargs+=( -e "s|^STEPB=.*|STEPB=\"/path/to/stepB_script.py\"  # Optional StepB: Python script to run between step2 and step3.|" )
    fi
    if [[ $fc == "t" ]]; then
        sedargs+=( -e "s|^STEPC=.*|$(hook_line STEPC "$stepc")|" )
    else
        sedargs+=( -e "s|^STEPC=.*|STEPC=\"/path/to/stepC_script.py\"  # Optional StepC: Python script to run between step3 and step4.|" )
    fi

    # STEP4 is patched only where step4 actually runs: leaving the template's
    # command alone in the step1-2 script keeps an unused seed out of it.
    if [[ $s4 == "t" ]]; then
        local step4cmd="$STEP4_BASE"
        [[ $USE_EXPLICIT_SEED -eq 1 ]] && step4cmd="$step4cmd -u MONTE_SEED=$SEED"
        sedargs+=( -e "s|^STEP4=.*|STEP4=\"${step4cmd}\"|" )
    fi

    sed "${sedargs[@]}" "$MASTER_SUBMIT_TPL" > "$out"
    chmod +x "$out"

    # Fail loudly rather than silently submitting a half-patched script.
    local bad=0
    grep -q "^#SBATCH --job-name=${job}$"  "$out" || bad=1
    grep -q "^step3=\"${s3}\""             "$out" || bad=1
    grep -q "^#SBATCH -o ${logname}$"      "$out" || bad=1
    grep -q "^CPUS=${cpus} "                "$out" || bad=1
    grep -q "^stepB=\"${fb}\""             "$out" || bad=1
    grep -q "^stepC=\"${fc}\""             "$out" || bad=1
    [[ $fb == "t" ]] && { grep -q "^STEPB=\"${stepb}\"" "$out" || bad=1; }
    [[ $fc == "t" ]] && { grep -q "^STEPC=\"${stepc}\"" "$out" || bad=1; }
    [[ $s4 == "t" && $USE_EXPLICIT_SEED -eq 1 ]] && \
        { grep -q "MONTE_SEED=$SEED" "$out" || bad=1; }
    if [[ $bad -eq 1 ]]; then
        echo -e "${RED}[FATAL] $out did not patch cleanly -- template format changed?${RESET}"; exit 1
    fi
    local bname="-" cname="-"
    [[ $fb == "t" ]] && bname=$(basename "$stepb")
    [[ $fc == "t" ]] && cname=$(basename "$stepc")
    printf "${GREEN}[GEN]${RESET} %-28s job=%-16s steps=%s%s%s%s center=%s cpus=%s\n" \
        "${out#$TRIAL/}" "$job" "$s1" "$s2" "$s3" "$s4" "$center" "$cpus"
    printf "      stepB=%-26s stepC=%s\n" "$bname" "$cname"
}

mkdir -p "$TRIAL"/{run_holo,run_apo,run_inhib}

# ---------------------------------------------------------------- PDB folders
# Only *.pdb -- kin-pdb also holds a p_info/ that pro_batch regenerates itself.
for d in run_holo run_apo; do
    mkdir -p "$TRIAL/$d/kin-pdb"
    cp "$MASTER_KIN_PDB"/*.pdb "$TRIAL/$d/kin-pdb/"
done
mkdir -p "$TRIAL/run_inhib/cof-pdb"
cp "$MASTER_COF_PDB"/*.pdb "$TRIAL/run_inhib/cof-pdb/"
echo -e "${GREEN}[PDB]${RESET} kin-pdb -> run_holo, run_apo ($(ls "$TRIAL/run_holo/kin-pdb"/*.pdb | wc -l) files); cof-pdb -> run_inhib ($(ls "$TRIAL/run_inhib/cof-pdb"/*.pdb | wc -l) files)"

# Symlinked, not copied: one list at the repo root is the single source of
# truth, so trials cannot drift apart on the kinase/inhibitor mapping.
ln -sfn "$(realpath --relative-to="$TRIAL" "$MASTER_INHIB_LST")" "$TRIAL/pdb_inhibitor.lst"

# Copied, not symlinked: the trial keeps the versions it was run with, the same
# way it keeps its own generated submit scripts.  Both numbered scripts sit at
# the trial root so they read as one sequence: 0- then 1-.
cp "$PREPARE" "$TRIAL/0-prepare_run_apo.py"
cp "$XTSRUN"  "$TRIAL/1-run_xts_corr.py"
cp "$FIG3"    "$TRIAL/2-plot_sumcrg_inhibitors_xts_Fig3.py"
cp "$FIG4A"   "$TRIAL/3-plot_sumcrg_comparison_xts_Fig4A.py"
chmod +x "$TRIAL"/[0-3]-*.py
echo -e "${GREEN}[GEN]${RESET} 0-prepare_run_apo.py                   (from ${PREPARE#$ROOT/})"
echo -e "${GREEN}[GEN]${RESET} 1-run_xts_corr.py                      (from ${XTSRUN#$ROOT/})"
echo -e "${GREEN}[GEN]${RESET} 2-plot_sumcrg_inhibitors_xts_Fig3.py   (from ${FIG3#$ROOT/})"
echo -e "${GREEN}[GEN]${RESET} 3-plot_sumcrg_comparison_xts_Fig4A.py  (from ${FIG4A#$ROOT/})"

# stepB does a different job in each tree, so each script names its own:
#   holo s1s2 -> make_holo_apo_step2_out.py  split step2_out.pdb into the
#                                            holo_ and apo_ variants
#   apo  s3s4 -> install_apo_step2_out.py    check the pair, point step2_out.pdb
#                                            at apo_step2_out.pdb
# stepC (prune_kin-inhib_head3.py) runs everywhere step4 does.  "-" = hook off.
#                out                                  job              1 2 3 4  ctr  stepB       stepC      cpus
gen_submit "$TRIAL/run_holo/submit_mcce4_s1s2.sh"  "${TAG}_holo_s1s2"  t t f f  t    "$VARIANTS" "-"        $CPUS_DEFAULT
gen_submit "$TRIAL/run_holo/submit_mcce4_s3s4.sh"  "${TAG}_holo_s3s4"  f f t t  f    "-"         "$PRUNE"   $CPUS_S34
gen_submit "$TRIAL/run_apo/submit_mcce4_s3s4.sh"   "${TAG}_apo_s3s4"   f f t t  f    "$INSTALL"  "$PRUNE"   $CPUS_S34
gen_submit "$TRIAL/run_inhib/submit_mcce4.sh"      "${TAG}_inhib"      t t t t  t    "-"         "$PRUNE"   $CPUS_DEFAULT



write_runbook
echo -e "${GREEN}[GEN]${RESET} RUNBOOK.md"
echo -e "\n${CYAN}Layout:${RESET}"
find "$TRIAL" -maxdepth 2 -not -path '*/kin-pdb/*' -not -path '*/cof-pdb/*' | sed "s|$TRIAL|Trial$(printf '%02d' "$N")|" | sort
echo -e "\n${GREEN}Done.${RESET}  Next: see $TRIAL/RUNBOOK.md"
