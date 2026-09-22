#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: rm_inhib_step2_out.py
Builds apo's step2_out.pdb from holo_step2_out.pdb by deleting the kinase
inhibitor (the "stepB" hook in run_apo/submit_mcce4_s3s4.sh).  Python port of
the legacy rm_cofs_step2_out.sh.

THE FLOW:
=========
    run_holo/<ID>/step2_out.pdb          what holo's steps 1-2 produced
              |                          (holo keeps its own pristine copy as
              |                           BK_holo_step2_out.pdb, via stepB there)
              |  1-prepare_run_apo.py    seeds apo, naming it holo_step2_out.pdb
              v
    run_apo/<ID>/holo_step2_out.pdb      stepB's INPUT -- never modified
              |
              |  stepB (this script)     delete the inhibitor
              v
    run_apo/<ID>/step2_out.pdb           what apo's step3 reads
    run_apo/<ID>/inhib_step2_out.pdb     the lines that were deleted

Input and output are different files, so re-running rebuilds step2_out.pdb from
an untouched source instead of editing it in place.  Deleting the ligand is the
ONLY difference between holo's step2_out.pdb and apo's: step1 and step2 are off
in run_apo, so the pocket is never repacked and the side chains stay exactly
where the inhibitor left them.  That is what isolates its electrostatics.

FAILING CLOSED:
===============
A seeded apo directory holds NO step2_out.pdb until this script writes one, and
step3 (with step2="f") only runs when step2_out.pdb exists.  So a stepB failure
means step3 skips the directory -- it cannot quietly compute holo and label it
apo, even though driver_mcce4.sh logs a stepB failure without aborting the run.
The rule here is: never leave behind a step2_out.pdb this script wrote but could
not verify, and never delete one it did not write.

MATCHING:
=========
Anchored to the residue-name columns 18-20 (Python [17:20]), NOT an unanchored
grep for the 3-letter code.  The legacy rm_cofs_step2_out.sh grepped whole lines
against a fixed 18-code alternation that also carried "END", so it could match a
code in any column and it dropped the END record.

Only the one inhibitor belonging to this structure is removed -- its code comes
from pdb_inhibitor.lst, keyed on the run directory name.  Every other heteroatom
is deliberately KEPT: 2ITZ's _CL chloride stays in the apo structure, exactly as
the legacy script left it.

USAGE:
======
  Via submit_mcce4_s3s4.sh in run_apo (no arguments; cwd is the run directory):
      stepB="t"
      STEPB="/data/home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/rm_inhib_step2_out.py"
  The driver runs it as `$PYEX $STEPB > stepB.log`.

  Manually:
      ./rm_inhib_step2_out.py                 # rebuild ./step2_out.pdb
      ./rm_inhib_step2_out.py -d run_apo      # batch: every <PDBID>/ below run_apo
      ./rm_inhib_step2_out.py -d run_apo --dry-run
      ./rm_inhib_step2_out.py -c FMM          # force the code instead of the .lst
"""

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

RESNAME = slice(17, 20)          # PDB residue-name columns 18-20
SOURCE = "holo_step2_out.pdb"    # stepB's input: holo's step2 output, untouched
STEP2 = "step2_out.pdb"          # stepB's output: what apo's step3 reads
REMOVED = "inhib_step2_out.pdb"  # the lines that were taken out
LST_NAME = "pdb_inhibitor.lst"
PDBID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

GREEN, YELLOW, CYAN, RED, RESET = (
    ("\033[0;32m", "\033[1;33m", "\033[0;36m", "\033[0;31m", "\033[0m")
    if sys.stdout.isatty() else ("", "", "", "", "")
)


def find_list(start, explicit=None):
    """pdb_inhibitor.lst: -l, then $INHIB_LST, then upward from the run directory."""
    for cand in (explicit, os.environ.get("INHIB_LST")):
        if cand:
            p = Path(cand).expanduser()
            if p.is_file():
                return p
            sys.exit(f"ERROR: {LST_NAME} not found at {p}")
    here = Path(start).resolve()
    for d in (here, *here.parents):
        p = d / LST_NAME
        if p.is_file():
            return p
    return None


def load_codes(path):
    """PDB -> inhibitor code.  Columns: PDB  Inhibitor  Inhibitor_Code, one header row."""
    codes = {}
    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        f = raw.split()
        if lineno == 1 or len(f) < 3:
            continue
        codes[f[0].upper()] = f[2]
    if not codes:
        sys.exit(f"ERROR: {path} holds no PDB -> inhibitor rows.")
    return codes


def count_code(path, code):
    return sum(1 for ln in path.read_text().splitlines() if ln[RESNAME] == code)


def resolve_source(run_dir, pid, code):
    """
    stepB's input is holo_step2_out.pdb.  A directory copied by hand may instead
    hold holo's step2_out.pdb under its original name; promote that once, so the
    input file exists and is never edited afterwards.
    Returns (source_path, status) where status is None, 'ok' or 'failed'.
    """
    src, step2 = run_dir / SOURCE, run_dir / STEP2

    if src.is_file():
        return src, None
    if not step2.is_file():
        print(f"{RED}[MISSING]  {pid}: neither {SOURCE} nor {STEP2} -- has holo's "
              f"step2 finished, and did 1-prepare_run_apo.py seed this directory?"
              f"{RESET}")
        return None, "failed"

    if count_code(step2, code) == 0:
        # No input to rebuild from, but the existing output provably holds no
        # inhibitor, so step3 would compute apo.  Leave it alone rather than
        # destroy a finished run; say clearly that the provenance is gone.
        print(f"{YELLOW}[NO SOURCE] {pid}: {SOURCE} is missing and {STEP2} already "
              f"has no {code} lines{RESET}")
        print(f"            leaving it as is -- re-seed with 1-prepare_run_apo.py "
              f"to restore {SOURCE}")
        return None, "ok"

    shutil.copy2(step2, src)
    print(f"{CYAN}[PROMOTED] {pid}: no {SOURCE}, so this run's {STEP2} was copied to "
          f"it as the untouched input{RESET}")
    return src, None


def strip_one(run_dir, codes, forced_code=None, dry=False):
    """Rebuild one apo step2_out.pdb.  Returns 'ok' | 'failed'."""
    run_dir = Path(run_dir).resolve()
    pid = run_dir.name.upper()
    step2, removed = run_dir / STEP2, run_dir / REMOVED

    code = forced_code or codes.get(pid)
    if not code:
        print(f"{RED}[NO CODE]  {pid}: not in {LST_NAME} (use -c CODE){RESET}")
        return "failed"

    src, status = resolve_source(run_dir, pid, code)
    if status is not None:
        return status

    lines = src.read_text().splitlines()
    hit = [ln for ln in lines if ln[RESNAME] == code]
    if not hit:
        print(f"{RED}[EMPTY]    {pid}: {SOURCE} holds no {code} lines -- it is not "
              f"holo's step2 output.  Re-seed with 1-prepare_run_apo.py.{RESET}")
        return "failed"
    keep = [ln for ln in lines if ln[RESNAME] != code]

    if dry:
        verb = "rebuild" if step2.is_file() else "write"
        print(f"{GREEN}[WOULD]{RESET}    {pid}  {code}: {verb} {STEP2} from {SOURCE}, "
              f"{len(hit)} of {len(lines)} lines removed, {len(keep)} kept")
        return "ok"

    body = "\n".join(keep) + "\n"
    existed = step2.is_file()
    unchanged = existed and step2.read_text() == body

    tmp = run_dir / (STEP2 + ".tmp")
    tmp.write_text(body)
    os.replace(tmp, step2)
    removed.write_text("\n".join(hit) + "\n")

    wrote = len(step2.read_text().splitlines())
    if wrote != len(keep):
        print(f"{RED}[FATAL]    {pid}: wrote {wrote} lines, expected {len(keep)} -- "
              f"removing {STEP2} so step3 skips this directory{RESET}")
        step2.unlink(missing_ok=True)
        removed.unlink(missing_ok=True)
        return "failed"

    tag = "[UNCHANGED]" if unchanged else ("[REBUILT] " if existed else "[STRIPPED]")
    print(f"{GREEN}{tag}{RESET} {pid}  {code}: {STEP2} built from {SOURCE}, "
          f"{len(hit)} removed, {wrote} kept  (deleted lines: {REMOVED})")

    # Other heteroatoms are kept on purpose; name them so a surprise is visible.
    het = sorted({ln[RESNAME] for ln in keep
                  if ln[:6] == "HETATM" and ln[RESNAME].strip()})
    if het:
        print(f"           kept heteroatoms: {', '.join(het)}")
    return "ok"


def main():
    ap = argparse.ArgumentParser(
        description="stepB: build apo's step2_out.pdb from holo_step2_out.pdb with "
                    "the inhibitor deleted.")
    ap.add_argument("-d", "--dir",
                    help="batch mode: every <PDBID>/ below this directory "
                         "(default: the current directory only, as stepB)")
    ap.add_argument("-l", "--list", help=f"{LST_NAME} (default: searched upward)")
    ap.add_argument("-c", "--code", help="inhibitor code, overriding the .lst lookup")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    base = Path(args.dir).resolve() if args.dir else Path.cwd()
    if not base.is_dir():
        sys.exit(f"ERROR: {base} is not a directory.")

    lst = find_list(base, args.list)
    codes = load_codes(lst) if lst else {}
    if not codes and not args.code:
        sys.exit(f"ERROR: {LST_NAME} not found at or above {base} (use -l or -c).")

    print(f"{CYAN}rm_inhib_step2_out.py  (stepB){RESET}")
    print(f"  dir  : {base}")
    print(f"  build: {SOURCE} -> {STEP2}")
    print(f"  lst  : {lst if lst else '(not used, -c given)'}"
          f"{'' if not codes else f'  [{len(codes)} structures]'}")
    if args.code:
        print(f"  code : {args.code}  (forced)")
    if args.dry_run:
        print(f"  {YELLOW}*** DRY RUN - nothing written ***{RESET}")
    print()

    if args.dir:
        targets = sorted(d for d in base.iterdir()
                         if d.is_dir() and PDBID_RE.match(d.name))
        if not targets:
            sys.exit(f"ERROR: no <PDBID>/ directories under {base}.")
    else:
        targets = [base]

    tally = {"ok": 0, "failed": 0}
    for t in targets:
        tally[strip_one(t, codes, args.code, args.dry_run)] += 1

    print(f"\n{CYAN}{'='*60}{RESET}")
    print(f"built: {tally['ok']}   failed: {tally['failed']}")
    sys.exit(1 if tally["failed"] else 0)


if __name__ == "__main__":
    main()
