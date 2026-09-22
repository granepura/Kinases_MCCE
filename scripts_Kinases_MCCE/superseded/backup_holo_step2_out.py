#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: backup_holo_step2_out.py
Keeps a pristine copy of holo's step2_out.pdb as BK_holo_step2_out.pdb (the
"stepB" hook in run_holo/submit_mcce4_s1s2.sh).

PURPOSE:
========
run_holo's step2_out.pdb is the origin of the whole comparison: run_apo is that
same file with the ligand deleted, and every holo/apo difference is read against
it.  It is also the one file a re-run of step1/step2 would replace in place.
BK_holo_step2_out.pdb is the untouched record of what steps 1-2 actually produced
in THIS run, so holo's own input can always be recovered and checked.

WHERE IT RUNS:
==============
run_holo/submit_mcce4_s1s2.sh:  stepB="t", STEPB=<this script>.  step3 and step4
are off in that script, so stepB is the last thing the job does -- it sees the
final step2_out.pdb.  The driver runs it as `$PYEX $STEPB > stepB.log` with the
run directory as cwd.

Note the neighbouring names, which are different files:
    run_holo/<ID>/step2_out.pdb          what steps 1-2 produced
    run_holo/<ID>/BK_holo_step2_out.pdb  pristine copy of it (this script)
    run_apo/<ID>/holo_step2_out.pdb      the same content, seeded into apo by
                                         1-prepare_run_apo.py -- stepB's INPUT
    run_apo/<ID>/step2_out.pdb           apo: built from it, ligand deleted
    run_apo/<ID>/inhib_step2_out.pdb     apo: the deleted inhibitor lines

REFRESH RULES:
==============
The point of the file is to match the step2_out.pdb this job just produced, so a
re-run of steps 1-2 refreshes it -- a backup silently describing an older run is
the failure mode this is meant to prevent (compare head3.lst_BK, which is never
refreshed and goes stale for exactly that reason).  The log always names both
sha256 values, so a replacement is visible in stepB.log rather than silent.
Pass --no-clobber to keep an existing copy instead.

USAGE:
======
  Via submit_mcce4_s1s2.sh (no arguments; cwd is the run directory):
      stepB="t"
      STEPB="/data/home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/backup_holo_step2_out.py"

  Manually:
      ./backup_holo_step2_out.py                # back up ./step2_out.pdb
      ./backup_holo_step2_out.py -d run_holo    # batch: every <PDBID>/ below run_holo
      ./backup_holo_step2_out.py --dry-run
      ./backup_holo_step2_out.py --no-clobber   # never replace an existing copy
"""

import argparse
import hashlib
import re
import shutil
import sys
from pathlib import Path

SRC = "step2_out.pdb"
DST = "BK_holo_step2_out.pdb"
PDBID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

GREEN, YELLOW, CYAN, RED, RESET = (
    ("\033[0;32m", "\033[1;33m", "\033[0;36m", "\033[0;31m", "\033[0m")
    if sys.stdout.isatty() else ("", "", "", "", "")
)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def backup_one(run_dir, src_name, dst_name, dry=False, clobber=True):
    """Back up one run directory.  Returns 'ok' | 'skipped' | 'failed'."""
    run_dir = Path(run_dir).resolve()
    pid = run_dir.name
    src, dst = run_dir / src_name, run_dir / dst_name

    if not src.is_file():
        print(f"{RED}[MISSING]  {pid}: no {src_name} -- did step2 finish?{RESET}")
        return "failed"

    src_sum = sha256(src)
    lines = sum(1 for _ in src.open("rb"))

    if dst.is_file():
        dst_sum = sha256(dst)
        if dst_sum == src_sum:
            print(f"{GREEN}[CURRENT]{RESET}  {pid}: {dst_name} already matches "
                  f"{src_name} (sha {src_sum[:12]}, {lines} lines)")
            return "ok"
        if not clobber:
            print(f"{YELLOW}[KEPT]     {pid}: {dst_name} differs from {src_name} but "
                  f"--no-clobber is set{RESET}")
            print(f"           existing sha {dst_sum[:12]} != step2 sha {src_sum[:12]} "
                  f"-- the backup describes an older run")
            return "skipped"
        action = (f"replaced (was sha {dst_sum[:12]}, step1-2 appear to have been "
                  f"re-run)")
    else:
        action = "written"

    if dry:
        print(f"{GREEN}[WOULD]{RESET}    {pid}: {dst_name} {action}, "
              f"sha {src_sum[:12]}, {lines} lines")
        return "ok"

    tmp = run_dir / (dst_name + ".tmp")
    shutil.copy2(src, tmp)
    tmp.replace(dst)

    if sha256(dst) != src_sum:
        print(f"{RED}[FATAL]    {pid}: {dst_name} does not match {src_name} after "
              f"the copy -- removing it{RESET}")
        dst.unlink(missing_ok=True)
        return "failed"

    print(f"{GREEN}[BACKUP]{RESET}   {pid}: {dst_name} {action}, sha {src_sum[:12]}, "
          f"{lines} lines")
    return "ok"


def main():
    ap = argparse.ArgumentParser(
        description="stepB for run_holo: copy step2_out.pdb to BK_holo_step2_out.pdb, "
                    "the pristine record of what steps 1-2 produced.")
    ap.add_argument("-d", "--dir",
                    help="batch mode: every <PDBID>/ below this directory "
                         "(default: the current directory only, as stepB)")
    ap.add_argument("-f", "--file", default=SRC, help=f"source file (default: {SRC})")
    ap.add_argument("-o", "--out", default=DST, help=f"backup name (default: {DST})")
    ap.add_argument("--no-clobber", action="store_true",
                    help="keep an existing backup even when it no longer matches")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    base = Path(args.dir).resolve() if args.dir else Path.cwd()
    if not base.is_dir():
        sys.exit(f"ERROR: {base} is not a directory.")

    print(f"{CYAN}backup_holo_step2_out.py  (stepB){RESET}")
    print(f"  dir : {base}")
    print(f"  copy: {args.file} -> {args.out}")
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

    tally = {"ok": 0, "skipped": 0, "failed": 0}
    for t in targets:
        tally[backup_one(t, args.file, args.out, args.dry_run,
                         not args.no_clobber)] += 1

    print(f"\n{CYAN}{'='*60}{RESET}")
    print(f"backed up: {tally['ok']}   skipped: {tally['skipped']}   "
          f"failed: {tally['failed']}")
    sys.exit(1 if tally["failed"] else 0)


if __name__ == "__main__":
    main()
