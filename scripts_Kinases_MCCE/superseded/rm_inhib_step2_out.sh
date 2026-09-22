#!/bin/bash
#
# rm_inhib_step2_out.sh -- strip the kinase inhibitor from step2_out.pdb.
#
# SUPERSEDED, and it no longer matches the current layout -- do not run it.
# rm_inhib_step2_out.py is what run_apo actually runs, as the driver's stepB
# hook, and it goes the other way round: it READS holo_step2_out.pdb (seeded by
# 1-prepare_run_apo.py) and WRITES step2_out.pdb, so its input is never edited.
# This script instead edits step2_out.pdb in place and writes the backup itself.
# The resulting step2_out.pdb is byte-identical; kept only as the reference
# implementation of the column 18-20 matching.  Edit the .py, not this.
#
# Adapted from rm_cofs_step2_out.sh, with three deliberate changes:
#   1. Inhibitor codes are read from pdb_inhibitor.lst instead of a hardcoded
#      list, so the list has one source of truth.
#   2. Matching is anchored to the residue-name columns (18-20) instead of an
#      unanchored `egrep "$COFACTORS"` over the whole line.  The old pattern
#      also carried "END", so it removed the END record and could in principle
#      match a code appearing anywhere else on a line.
#   3. The backup is holo_step2_out.pdb (was run_kin_step2_out.pdb) to match the
#      run_holo / run_apo / run_inhib naming.
#
# Non-inhibitor heteroatoms are deliberately KEPT -- e.g. 2ITZ's _CL chloride
# stays in the apo structure, exactly as the original script left it.
#
# Run from inside run_apo (loops over */), or pass -d <dir>.
#
#   ./rm_inhib_step2_out.sh [-l pdb_inhibitor.lst] [-d run_apo] [--dry-run]
#
set -u

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; RED='\033[0;31m'; RESET='\033[0m'

LST=""; BASE="."; DRY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -l|--list) LST="$2"; shift 2 ;;
        -d|--dir)  BASE="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

BASE=$(cd "$BASE" && pwd) || exit 1
[[ -n "$LST" ]] || LST="$BASE/pdb_inhibitor.lst"
[[ -f "$LST" ]] || LST="$BASE/../pdb_inhibitor.lst"
[[ -f "$LST" ]] || { echo -e "${RED}[FATAL] pdb_inhibitor.lst not found (use -l)${RESET}"; exit 1; }

LOG="$BASE/rm_inhib_step2_out.log"
exec > >(tee -a "$LOG") 2>&1

echo -e "\n${CYAN}=== rm_inhib_step2_out.sh  $(date '+%Y-%m-%d %H:%M:%S') ===${RESET}"
echo -e "${CYAN}dir: $BASE${RESET}"
echo -e "${CYAN}lst: $LST${RESET}"
[[ $DRY -eq 1 ]] && echo -e "${YELLOW}*** DRY RUN - nothing written ***${RESET}"

# PDB -> inhibitor code
declare -A CODE
while read -r pdb code; do
    [[ -n "$pdb" && -n "$code" ]] && CODE["$pdb"]="$code"
done < <(awk 'NR>1 && NF>=3 {print toupper($1), $3}' "$LST")
echo -e "${CYAN}inhibitor codes loaded: ${#CODE[@]}${RESET}\n"

ok=0; skip=0; fail=0
for dir in "$BASE"/*/; do
    id=$(basename "${dir%/}")
    [[ "$id" =~ ^[0-9][A-Za-z0-9]{3}$ ]] || continue

    step2="$dir/step2_out.pdb"
    bak="$dir/holo_step2_out.pdb"
    code="${CODE[${id^^}]:-}"

    if [[ -z "$code" ]]; then
        echo -e "${RED}[NO CODE]  $id not in pdb_inhibitor.lst${RESET}"; ((fail++)); continue
    fi
    if [[ -f "$bak" ]]; then
        echo -e "${YELLOW}[SKIP]     $id already stripped (holo_step2_out.pdb exists)${RESET}"; ((skip++)); continue
    fi
    if [[ ! -f "$step2" ]]; then
        echo -e "${RED}[MISSING]  $id has no step2_out.pdb${RESET}"; ((fail++)); continue
    fi

    n=$(awk -v c="$code" 'substr($0,18,3)==c' "$step2" | wc -l)
    if [[ "$n" -eq 0 ]]; then
        echo -e "${RED}[EMPTY]    $id has no $code lines -- already stripped, or wrong source${RESET}"
        ((fail++)); continue
    fi
    total=$(wc -l < "$step2")

    if [[ $DRY -eq 1 ]]; then
        echo -e "${GREEN}[WOULD]${RESET}    $id  remove $n $code lines of $total"
        ((ok++)); continue
    fi

    cp -p "$step2" "$bak"
    awk -v c="$code" 'substr($0,18,3)!=c' "$bak" > "$step2"
    kept=$(wc -l < "$step2")

    if [[ $(( total - kept )) -ne $n ]]; then
        echo -e "${RED}[FATAL]    $id removed $(( total - kept )) lines, expected $n -- restoring${RESET}"
        cp -p "$bak" "$step2"; rm -f "$bak"; ((fail++)); continue
    fi
    echo -e "${GREEN}[STRIPPED]${RESET} $id  $code: $n removed, $kept kept  (backup: holo_step2_out.pdb)"
    ((ok++))
done

echo -e "\n${CYAN}=====================================${RESET}"
echo -e "${GREEN}Done.${RESET}  stripped: $ok   skipped: $skip   failed: $fail"
echo -e "Log: ${YELLOW}$LOG${RESET}"
[[ $fail -gt 0 ]] && exit 1 || exit 0
