#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: install_apo_step2_out.py
Checks apo's step2 pair and makes sure step2_out.pdb points at apo_step2_out.pdb,
right before apo's step3 (the "stepB" hook in run_apo/submit_mcce4_s3s4.sh).

THE FLOW:
=========
    run_holo/<ID>/step2_out.pdb          what holo's steps 1-2 produced
         |  make_holo_apo_step2_out.py   holo's stepB splits it in two
         +-> holo_step2_out.pdb          exact copy, inhibitor included
         +-> apo_step2_out.pdb           same file, inhibitor deleted
                    |
                    |  0-prepare_run_apo.py   copies the whole directory across
                    v
    run_apo/<ID>/apo_step2_out.pdb       the apo structure -- never modified
                    ^
                    |  step2_out.pdb -> apo_step2_out.pdb
                    |  (relative symlink, made by 0-prepare_run_apo.py;
                    |   this script verifies it and resets it if needed)
    run_apo/<ID>/step2_out.pdb           the name apo's step3 reads

Both trees therefore hold the same holo_step2_out.pdb and apo_step2_out.pdb;
only step2_out.pdb differs between them, and neither job ever edits a file it
also reads.

WHAT IT CHECKS BEFORE INSTALLING:
=================================
  * apo_step2_out.pdb holds no lines for this structure's inhibitor (columns
    18-20, from pdb_inhibitor.lst keyed on the directory name),
  * holo_step2_out.pdb, when present, DOES hold them, and apo is exactly holo
    minus those lines -- nothing else moved, nothing else dropped.
  * step2_out.pdb is a link to apo_step2_out.pdb.  If it is missing, a plain
    file, or pointing somewhere else, it is reset to the correct relative link.
That is a last-moment integrity check on the pair, immediately before the 8
minutes of step3 that depend on it.

FAILING CLOSED:
===============
step3 with step2="f" runs only when step2_out.pdb exists, and driver_mcce4.sh
logs a stepB failure without aborting the run.  So when the pair does not check
out, this script removes the link rather than leaving step3 a structure it could
not vouch for: that structure ends without a pK.out and pro_batch --check flags
it.  apo_step2_out.pdb and holo_step2_out.pdb are never modified here -- fix
them by re-splitting in run_holo and re-seeding.

USAGE:
======
  Via submit_mcce4_s3s4.sh in run_apo (no arguments; cwd is the run directory):
      stepB="t"
      STEPB="/data/home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/install_apo_step2_out.py"
  The driver runs it as `$PYEX $STEPB > stepB.log`.

  Manually:
      ./install_apo_step2_out.py              # check ./apo_step2_out.pdb + link
      ./install_apo_step2_out.py -d run_apo   # batch: every <PDBID>/ below run_apo
      ./install_apo_step2_out.py --dry-run
"""

import argparse
import os
import re
import sys
from pathlib import Path

RESNAME = slice(17, 20)        # PDB residue-name columns 18-20
APO = "apo_step2_out.pdb"      # the input: holo's step2 output minus the ligand
HOLO = "holo_step2_out.pdb"    # the cross-check: holo's step2 output
STEP2 = "step2_out.pdb"        # the output: what apo's step3 reads
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


def link_state(run_dir):
    """What step2_out.pdb currently is: 'linked' | 'wrong' | 'file' | 'absent'."""
    link = run_dir / STEP2
    if link.is_symlink():
        return "linked" if os.readlink(link) == APO else "wrong"
    return "file" if link.exists() else "absent"


def drop_link(run_dir):
    """Remove step2_out.pdb so step3 skips this directory."""
    link = run_dir / STEP2
    if link.is_symlink() or link.exists():
        link.unlink()
        print(f"           {STEP2} removed -- step3 will skip this directory")


def install_one(run_dir, codes, forced_code=None, dry=False):
    """Check the pair and (re)link step2_out.pdb.  Returns 'ok' | 'failed'."""
    run_dir = Path(run_dir).resolve()
    pid = run_dir.name.upper()
    apo, holo, step2 = run_dir / APO, run_dir / HOLO, run_dir / STEP2

    code = forced_code or codes.get(pid)
    if not code:
        print(f"{RED}[NO CODE]  {pid}: not in {LST_NAME} (use -c CODE){RESET}")
        return "failed"
    if not apo.is_file():
        print(f"{RED}[MISSING]  {pid}: no {APO} -- holo's stepB "
              f"(make_holo_apo_step2_out.py) has not run, or this directory was not "
              f"seeded by 0-prepare_run_apo.py{RESET}")
        if not dry:
            drop_link(run_dir)
        return "failed"

    apo_text = apo.read_text()
    apo_lines = apo_text.splitlines()

    left = sum(1 for ln in apo_lines if ln[RESNAME] == code)
    if left:
        print(f"{RED}[BAD APO]  {pid}: {APO} still holds {left} {code} lines -- it is "
              f"not an apo structure.  Re-split it in run_holo.{RESET}")
        if not dry:
            drop_link(run_dir)
        return "failed"

    # Cross-check the pair: apo must be holo minus exactly the inhibitor lines.
    checked = ""
    if holo.is_file():
        holo_lines = holo.read_text().splitlines()
        hit = [ln for ln in holo_lines if ln[RESNAME] == code]
        expect = [ln for ln in holo_lines if ln[RESNAME] != code]
        if not hit:
            print(f"{RED}[BAD HOLO] {pid}: {HOLO} holds no {code} lines -- the pair "
                  f"is not a holo/apo pair.  Re-split it in run_holo.{RESET}")
            if not dry:
                drop_link(run_dir)
            return "failed"
        if expect != apo_lines:
            print(f"{RED}[MISMATCH] {pid}: {APO} is not {HOLO} minus its {len(hit)} "
                  f"{code} lines ({len(apo_lines)} lines vs {len(expect)} expected) -- "
                  f"one of them has been edited.  Re-split it in run_holo.{RESET}")
            if not dry:
                drop_link(run_dir)
            return "failed"
        checked = f", verified against {HOLO} ({len(hit)} {code} lines absent)"
    else:
        print(f"{YELLOW}[NO HOLO]  {pid}: {HOLO} absent, installing {APO} without the "
              f"pair cross-check{RESET}")

    state = link_state(run_dir)
    if dry:
        verb = {"linked": "leave", "wrong": "reset", "file": "replace",
                "absent": "create"}[state]
        print(f"{GREEN}[WOULD]{RESET}    {pid}: {verb} {STEP2} -> {APO}, "
              f"{len(apo_lines)} lines{checked}")
        return "ok"

    if state != "linked":
        if state != "absent":
            (run_dir / STEP2).unlink()
        (run_dir / STEP2).symlink_to(APO)     # relative on purpose

    link = run_dir / STEP2
    if not (link.is_symlink() and os.readlink(link) == APO and link.is_file()
            and link.read_text() == apo_text):
        print(f"{RED}[FATAL]    {pid}: {STEP2} does not resolve to {APO} -- removing "
              f"it so step3 skips this directory{RESET}")
        drop_link(run_dir)
        return "failed"

    tag = {"linked": "[LINKED]  ", "wrong": "[RELINKED]",
           "file": "[RELINKED]", "absent": "[LINKED]  "}[state]
    note = {"linked": "", "wrong": " (was pointing elsewhere)",
            "file": " (was a plain file)", "absent": ""}[state]
    print(f"{GREEN}{tag}{RESET} {pid}: {STEP2} -> {APO}, {len(apo_lines)} lines"
          f"{checked}{note}")
    return "ok"


def main():
    ap = argparse.ArgumentParser(
        description="stepB for run_apo: check the step2 pair and point step2_out.pdb "
                    "at apo_step2_out.pdb.")
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

    print(f"{CYAN}install_apo_step2_out.py  (stepB){RESET}")
    print(f"  dir    : {base}")
    print(f"  link   : {STEP2} -> {APO}")
    print(f"  lst    : {lst if lst else '(not used, -c given)'}"
          f"{'' if not codes else f'  [{len(codes)} structures]'}")
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
        tally[install_one(t, codes, args.code, args.dry_run)] += 1

    print(f"\n{CYAN}{'='*60}{RESET}")
    print(f"linked: {tally['ok']}   failed: {tally['failed']}")
    sys.exit(1 if tally["failed"] else 0)


if __name__ == "__main__":
    main()
