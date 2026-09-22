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

N="${1:-}"; FORCE=0
[[ "${2:-}" == "--force" ]] && FORCE=1
[[ "$N" =~ ^[0-9]+$ ]] || { echo "Usage: $0 <trial-number> [--force]"; exit 1; }

TRIAL=$(trial_dir "$N")
TAG=$(printf "T%02d" "$N")
SEED=$(trial_seed "$N")

echo -e "${CYAN}=== setup_trial.sh  ->  $TRIAL ===${RESET}"

for src in "$MASTER_KIN_PDB" "$MASTER_COF_PDB" "$MASTER_INHIB_LST" "$MASTER_SUBMIT_TPL" \
           "$PRUNE" "$VARIANTS" "$INSTALL" "$PREPARE" "$XTSRUN"; do
    [[ -e "$src" ]] || { echo -e "${RED}[FATAL] missing master: $src${RESET}"; exit 1; }
done

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
    local fb="f" fc="f"
    [[ $stepb != "-" ]] && fb="t"
    [[ $stepc != "-" ]] && fc="t"

    local -a sedargs=(
        -e "s|^#SBATCH --job-name=.*|#SBATCH --job-name=${job}|"
        -e "s|^step1=\"[tf]\"|step1=\"${s1}\"|"
        -e "s|^step2=\"[tf]\"|step2=\"${s2}\"|"
        -e "s|^step3=\"[tf]\"|step3=\"${s3}\"|"
        -e "s|^step4=\"[tf]\"|step4=\"${s4}\"|"
        -e "s|^center=\"[tf]\"|center=\"${center}\"|"
        -e "s|^stepB=\"[tf]\"|stepB=\"${fb}\"|"
        -e "s|^stepC=\"[tf]\"|stepC=\"${fc}\"|"
    )
    [[ $fb == "t" ]] && sedargs+=( -e "s|^STEPB=.*|$(hook_line STEPB "$stepb")|" )
    [[ $fc == "t" ]] && sedargs+=( -e "s|^STEPC=.*|$(hook_line STEPC "$stepc")|" )

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
    printf "${GREEN}[GEN]${RESET} %-28s job=%-16s steps=%s%s%s%s center=%s\n" \
        "${out#$TRIAL/}" "$job" "$s1" "$s2" "$s3" "$s4" "$center"
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

cp "$MASTER_INHIB_LST" "$TRIAL/pdb_inhibitor.lst"

# Copied, not symlinked: the trial keeps the versions it was run with, the same
# way it keeps its own generated submit scripts.  Both numbered scripts sit at
# the trial root so they read as one sequence: 0- then 1-.
cp "$PREPARE" "$TRIAL/0-prepare_run_apo.py"
cp "$XTSRUN"  "$TRIAL/1-run_xts_corr.py"
chmod +x "$TRIAL/0-prepare_run_apo.py" "$TRIAL/1-run_xts_corr.py"
echo -e "${GREEN}[GEN]${RESET} 0-prepare_run_apo.py  (from ${PREPARE#$ROOT/})"
echo -e "${GREEN}[GEN]${RESET} 1-run_xts_corr.py     (from ${XTSRUN#$ROOT/})"

# stepB does a different job in each tree, so each script names its own:
#   holo s1s2 -> make_holo_apo_step2_out.py  split step2_out.pdb into the
#                                            holo_ and apo_ variants
#   apo  s3s4 -> install_apo_step2_out.py    check the pair, point step2_out.pdb
#                                            at apo_step2_out.pdb
# stepC (prune_kin-inhib_head3.py) runs everywhere step4 does.  "-" = hook off.
#                out                                  job              1 2 3 4  ctr  stepB      stepC
gen_submit "$TRIAL/run_holo/submit_mcce4_s1s2.sh"  "${TAG}_holo_s1s2"  t t f f  t    "$VARIANTS" "-"
gen_submit "$TRIAL/run_holo/submit_mcce4_s3s4.sh"  "${TAG}_holo_s3s4"  f f t t  f    "-"        "$PRUNE"
gen_submit "$TRIAL/run_apo/submit_mcce4_s3s4.sh"   "${TAG}_apo_s3s4"   f f t t  f    "$INSTALL"  "$PRUNE"
gen_submit "$TRIAL/run_inhib/submit_mcce4.sh"      "${TAG}_inhib"      t t t t  t    "-"        "$PRUNE"

# ---------------------------------------------------------------- runbook
{
    echo "# Trial$(printf '%02d' "$N") runbook"
    echo
    echo "Created: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "MONTE_SEED: $([[ $USE_EXPLICIT_SEED -eq 1 ]] && echo "$SEED (explicit)" || echo '-1 (time-based)')"
    echo
    echo '## Shared scripts (sha256 at setup time -- re-check before comparing trials)'
    echo '```'
    sha256sum "$TRIAL/0-prepare_run_apo.py" "$TRIAL/1-run_xts_corr.py" \
              "$VARIANTS" "$INSTALL" "$PRUNE" \
              "$SCRIPT_DIR/trial_config.sh" 2>/dev/null
    echo '```'
    cat <<BODY

## Order of operations

1. holo, steps 1-2 -- builds the conformers and the coordinate frame everything else inherits
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s1s2.sh -job-name holo_s1s2 -j 15
   stepB here is make_holo_apo_step2_out.py: step3/step4 are off, so it runs last
   and splits the finished step2_out.pdb into holo_step2_out.pdb (exact copy) and
   apo_step2_out.pdb (inhibitor deleted).  Check each structure's stepB.log.

2. Seed apo from holo (needs only holo's split step2 variants, so it can run as
   soon as holo's steps 1-2 finish -- no need to wait for holo's step3/4).
   Copies each whole run_holo/<PDBID> to run_apo/<PDBID> -- every file steps 1-2
   left, plus both step2 variants -- replacing any that is already there, and
   skipping holo's step2_out.pdb and any step3/4 products.  It then links
   step2_out.pdb -> apo_step2_out.pdb, resetting the link if one exists.
       ./0-prepare_run_apo.py
       ./0-prepare_run_apo.py --dry-run   # inspect without writing
       ./0-prepare_run_apo.py 1XKK 2ITZ   # re-seed just these

3. holo, steps 3-4
       cd run_holo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name holo_s3s4 -j 15

4. apo, steps 3-4
       cd run_apo
       pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name apo_s3s4 -j 15

5. inhib, steps 1-4 (independent: starts from cof-pdb, not carved from holo)
       cd run_inhib
       pro_batch cof-pdb -custom submit_mcce4.sh -job-name inhib -j 15

6. entropy-correct all three trees, once their step4 has finished
       cd $TRIAL
       ./1-run_xts_corr.py                # xts_fort.38 + xts_sum_crg.out per structure
       ./1-run_xts_corr.py --dry-run
   Every published figure reads xts_sum_crg.out, and step4 does NOT produce it --
   xts_corr.py is a separate pass.  Run it in all three trees or you would be
   comparing corrected numbers against uncorrected ones.

Steps 3 and 5 are independent of each other and of step 4; run them concurrently.
Check progress with:  pro_batch --check -job-name <name>

## The step2 chain

The ligand is deleted once, in the holo job, and both trees then share the same
two files.  Only step2_out.pdb differs between them:

    run_holo/<ID>/step2_out.pdb          what steps 1-2 produced = holo
         |  holo stepB: make_holo_apo_step2_out.py
         +-> holo_step2_out.pdb          exact copy of it
         +-> apo_step2_out.pdb           same file, inhibitor deleted

    0-prepare_run_apo.py copies the whole directory across, then:

    run_apo/<ID>/holo_step2_out.pdb      copied, for reference and the pair check
    run_apo/<ID>/apo_step2_out.pdb       copied -- the apo structure
    run_apo/<ID>/step2_out.pdb  ->  apo_step2_out.pdb   (relative symlink)

Deleting the ligand is the ONLY difference between the two structures -- step1
and step2 are off in run_apo, so the pocket is never repacked.  The code comes
from pdb_inhibitor.lst keyed on the directory name and is matched on the
residue-name columns 18-20, so other heteroatoms stay (2ITZ keeps its _CL).

apo's stepB is install_apo_step2_out.py.  It rewrites nothing in the normal
case: it checks that apo_step2_out.pdb holds none of this structure's inhibitor,
that it is exactly holo_step2_out.pdb minus those lines, and that step2_out.pdb
points at it -- resetting the link if it is missing, a plain file, or pointing
elsewhere.  That check runs immediately before the ~8 minutes of step3 that
depend on it.

driver_mcce4.sh logs a stepB failure without aborting the run, so when the pair
does not check out the script REMOVES step2_out.pdb.  step3 with step2="f" only
runs when that file exists, so the structure ends without a pK.out and
pro_batch --check flags it, instead of step3 computing something unvouched-for.
Check stepB.log and mcce_timing.log.

## Why stepC runs in all three

prune_kin-inhib_head3.py does two independent edits: the inhibitor conformer
pruning AND forcing ARG positive (neutral ARG -> FL=t).  apo has no inhibitor but
does have ARG, so skipping stepC there would give holo and apo different ARG
treatments and invalidate the comparison.  inhib has no ARG and gets only the
inhibitor edits.  All three therefore run stepC="t".

## Frames

run_holo and run_apo share one coordinate frame: apo reuses holo's step2_out.pdb
and its submit script has step1/step2 off, so nothing re-centers it.
run_inhib is built from cof-pdb through its own step1/step2 and therefore sits in
its own centered frame.  That is fine for pKa/charge -- it is an isolated-ligand
reference state -- but do not compare its coordinates to holo's.
BODY
} > "$TRIAL/RUNBOOK.md"

echo -e "${GREEN}[GEN]${RESET} RUNBOOK.md"
echo -e "\n${CYAN}Layout:${RESET}"
find "$TRIAL" -maxdepth 2 -not -path '*/kin-pdb/*' -not -path '*/cof-pdb/*' | sed "s|$TRIAL|Trial$(printf '%02d' "$N")|" | sort
echo -e "\n${GREEN}Done.${RESET}  Next: see $TRIAL/RUNBOOK.md"
