#!/bin/bash
#
# derive_apo.sh -- build TrialNN/run_apo from TrialNN/run_holo.
#
# SUPERSEDED, and it no longer matches the current layout -- do not run it.
# prepare_run_apo.py, which setup_trial.sh drops into each trial as
# 1-prepare_run_apo.py, is what seeds run_apo now: same idea in Python (no rsync
# dependency), living inside the trial so the trial records the version it ran.
# It also names holo's step2_out.pdb as holo_step2_out.pdb in the apo directory
# and leaves no step2_out.pdb, which is what apo's stepB expects and what makes
# a failed strip unable to reach step3.  This script leaves step2_out.pdb in
# place instead.  Kept as the reference implementation; edit the .py, not this.
#
# For every structure that has a holo step2_out.pdb:
#   1. copy run_holo/<ID> -> run_apo/<ID>, EXCLUDING all step3/step4 artifacts
#      (see STEP34_ARTIFACTS in trial_config.sh).  This matters: head3.lst_BK is
#      never overwritten by prune_kin-inhib_head3.py once it exists, so an
#      inherited copy would leave apo holding holo's head3.lst as its "pristine
#      step3 output" backup.
#   2. write provenance.txt recording the source and its sha256, so a stale apo
#      built from an older holo is detectable instead of silent.
#
# It does NOT delete the inhibitor.  That happens inside the apo run, as stepB:
# run_apo/submit_mcce4_s3s4.sh has stepB="t" and STEPB=rm_inhib_step2_out.py, so
# each structure strips its own step2_out.pdb immediately before its own step3,
# and the strip is logged per structure in stepB.log next to the step it feeds.
# Pass --strip to do it here instead (stepB then reports SKIP -- it keys off the
# holo_step2_out.pdb backup -- so the two are safe to combine, not duplicated).
#
# Only holo's step2_out.pdb is needed, so this can run as soon as holo's steps
# 1-2 finish -- no need to wait for holo's step3/4.
#
# Usage:  ./derive_apo.sh 1 [--force] [--dry-run] [--strip]
#
set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/trial_config.sh"

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; RED='\033[0;31m'; RESET='\033[0m'

N="${1:-}"; FORCE=0; DRY=0; STRIP_HERE=0
shift || true
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)   FORCE=1; shift ;;
        --dry-run) DRY=1; shift ;;
        --strip)   STRIP_HERE=1; shift ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done
[[ "$N" =~ ^[0-9]+$ ]] || { echo "Usage: $0 <trial-number> [--force] [--dry-run] [--strip]"; exit 1; }

TRIAL=$(trial_dir "$N")
HOLO="$TRIAL/run_holo"
APO="$TRIAL/run_apo"
LOG="$TRIAL/derive_apo.log"

[[ -d "$HOLO" ]] || { echo -e "${RED}[FATAL] $HOLO not found -- run setup_trial.sh $N first${RESET}"; exit 1; }
mkdir -p "$APO"
exec > >(tee -a "$LOG") 2>&1

echo -e "\n${CYAN}=== derive_apo.sh  Trial$(printf '%02d' "$N")  $(date '+%Y-%m-%d %H:%M:%S') ===${RESET}"
[[ $DRY -eq 1 ]] && echo -e "${YELLOW}*** DRY RUN - nothing written ***${RESET}"

# rsync excludes are cheaper than copying 8.9M of energies/ per structure only to delete it.
EXCL=()
for a in "${STEP34_ARTIFACTS[@]}"; do EXCL+=( --exclude="$a" ); done
HAVE_RSYNC=0; command -v rsync >/dev/null && HAVE_RSYNC=1

ok=0; skip=0; fail=0
for d in "$HOLO"/*/; do
    id=$(basename "${d%/}")
    [[ "$id" =~ ^[0-9][A-Za-z0-9]{3}$ ]] || continue

    src_step2="$d/step2_out.pdb"
    dst="$APO/$id"

    if [[ ! -f "$src_step2" ]]; then
        echo -e "${YELLOW}[WAIT]     $id -- holo has no step2_out.pdb yet${RESET}"; ((skip++)); continue
    fi
    if [[ -d "$dst" && $FORCE -eq 0 ]]; then
        echo -e "${YELLOW}[SKIP]     $id -- run_apo/$id exists (use --force)${RESET}"; ((skip++)); continue
    fi

    sum=$(sha256sum "$src_step2" | awk '{print $1}')

    if [[ $DRY -eq 1 ]]; then
        echo -e "${GREEN}[WOULD]${RESET}    $id  copy holo -> apo (minus step3/4), src sha ${sum:0:12}"
        ((ok++)); continue
    fi

    [[ -d "$dst" ]] && rm -rf "$dst"
    mkdir -p "$dst"
    if [[ $HAVE_RSYNC -eq 1 ]]; then
        rsync -a "${EXCL[@]}" "$d" "$dst/" || { echo -e "${RED}[FAILED]   $id rsync${RESET}"; ((fail++)); continue; }
    else
        cp -a "$d." "$dst/"  || { echo -e "${RED}[FAILED]   $id cp${RESET}"; ((fail++)); continue; }
        # $a is unquoted on purpose: some entries are globs (submit_mcce4*.sh).
        for a in "${STEP34_ARTIFACTS[@]}"; do rm -rf ${dst:?}/$a; done
    fi

    # Guard: the purge list must actually have removed holo's step3 products.
    for a in energies head3.lst head3.lst_BK xts_sum_crg.out pK.out holo_step2_out.pdb; do
        if [[ -e "$dst/$a" ]]; then
            echo -e "${RED}[FATAL]    $id still holds $a after purge${RESET}"; ((fail++)); continue 2
        fi
    done

    cat > "$dst/provenance.txt" <<PROV
derived_by   : derive_apo.sh
derived_at   : $(date '+%Y-%m-%d %H:%M:%S')
source_dir   : $d
source_file  : step2_out.pdb
source_sha256: $sum
note         : apo reuses holo's step2 conformers; the inhibitor is deleted by
               stepB (rm_inhib_step2_out.py) at the start of the apo step3 job.
               The pocket is NOT repacked -- side chains stay in their holo
               positions, which is what isolates the ligand's electrostatic
               contribution.  Re-derive if source_sha256 no longer matches
               $d/step2_out.pdb
PROV

    echo -e "${GREEN}[SEEDED]${RESET}   $id  (holo step2 sha ${sum:0:12})"
    ((ok++))
done

echo -e "\n${CYAN}Seeded: $ok   skipped: $skip   failed: $fail${RESET}"

if [[ $STRIP_HERE -eq 1 && $ok -gt 0 && $DRY -eq 0 ]]; then
    echo -e "\n${CYAN}--- stripping inhibitors here (--strip) ---${RESET}"
    "$STRIP" -d "$APO" -l "$TRIAL/pdb_inhibitor.lst" || fail=$(( fail + 1 ))
elif [[ $STRIP_HERE -eq 1 && $DRY -eq 1 ]]; then
    echo -e "${YELLOW}(dry run: reporting the strip step only)${RESET}"
    "$STRIP" -d "$APO" -l "$TRIAL/pdb_inhibitor.lst" --dry-run || true
else
    echo -e "\n${CYAN}Inhibitors NOT stripped here -- run_apo's stepB does it in-job"
    echo -e "(rm_inhib_step2_out.py, between step2 and step3).${RESET}"
fi

echo -e "\n${GREEN}Done.${RESET}  Next: cd $APO && pro_batch kin-pdb -custom submit_mcce4_s3s4.sh -job-name apo_s3s4 -j 15"
[[ $fail -gt 0 ]] && exit 1 || exit 0
