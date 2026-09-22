#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: make_holo_apo_step2_out.py
Splits holo's finished step2_out.pdb into the two reference structures the whole
project compares (the "stepB" hook in run_holo/submit_mcce4_s1s2.sh):

    step2_out.pdb  --+--> holo_step2_out.pdb   exact copy, inhibitor included
                     |
                     +--> apo_step2_out.pdb    same file, inhibitor deleted

WHY BOTH ARE MADE HERE:
=======================
holo's step2_out.pdb is the origin of the comparison: run_apo is that same file
with the ligand deleted, sharing every conformer.  Making both variants in the
holo job, once, means the split happens where the step2 output is authoritative,
and run_apo never has to derive anything -- 0-prepare_run_apo.py copies both
files across and install_apo_step2_out.py puts the apo one in place as apo's
step2_out.pdb.  holo_step2_out.pdb also keeps a pristine record of what steps
1-2 produced, which a re-run of step1/step2 would otherwise overwrite in place.

step3 and step4 are off in submit_mcce4_s1s2.sh, so this runs last in that job
and sees the final step2_out.pdb.  The driver runs it as `$PYEX $STEPB >
stepB.log` with the run directory as cwd.

MATCHING:
=========
Anchored to the residue-name columns 18-20 (Python [17:20]), NOT an unanchored
grep for the 3-letter code.  The legacy rm_cofs_step2_out.sh grepped whole lines
against a fixed 18-code alternation that also carried "END", so it could match a
code in any column and it dropped the END record.

Only the one inhibitor belonging to this structure is removed -- its code comes
from pdb_inhibitor.lst, keyed on the run directory name.  Every other heteroatom
is deliberately KEPT: 2ITZ's _CL chloride stays in the apo structure, exactly as
the legacy script left it.  The deleted lines are not written out separately;
`diff holo_step2_out.pdb apo_step2_out.pdb` recovers them.

REFRESH RULES:
==============
Both variants are rebuilt on every run, so a re-run of steps 1-2 cannot leave
them describing an older structure (compare head3.lst_BK, which is never
refreshed and goes stale for exactly that reason).  The log names the sha256 of
each file, so a replacement is visible in stepB.log rather than silent.
Pass --no-clobber to keep existing variants instead.

FAILING CLOSED:
===============
If step2_out.pdb holds no lines for this structure's inhibitor, the two variants
would be identical and apo would silently equal holo.  That is refused, and any
half-written variant is removed, so 0-prepare_run_apo.py then refuses to seed
the structure at all.

USAGE:
======
  Via submit_mcce4_s1s2.sh (no arguments; cwd is the run directory):
      stepB="t"
      STEPB="/data/home/granepura/5-Kinases/Kinases_MCCE/scripts_Kinases_MCCE/make_holo_apo_step2_out.py"

  Manually:
      ./make_holo_apo_step2_out.py                # split ./step2_out.pdb
      ./make_holo_apo_step2_out.py -d run_holo    # batch: every <PDBID>/ below run_holo
      ./make_holo_apo_step2_out.py --dry-run
      ./make_holo_apo_step2_out.py -c FMM         # force the code instead of the .lst
"""

import argparse
import hashlib
import os
import re
import shutil
import sys
from pathlib import Path

RESNAME = slice(17, 20)        # PDB residue-name columns 18-20
STEP2 = "step2_out.pdb"        # holo's step2 output, the input here
HOLO = "holo_step2_out.pdb"    # exact copy of it
APO = "apo_step2_out.pdb"      # same, inhibitor deleted
LST_NAME = "pdb_inhibitor.lst"
PDBID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

GREEN, YELLOW, CYAN, RED, RESET = (
    ("\033[0;32m", "\033[1;33m", "\033[0;36m", "\033[0;31m", "\033[0m")
    if sys.stdout.isatty() else ("", "", "", "", "")
)


def sha256_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


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


def write_variant(path, text, label, clobber):
    """Write one variant.  Returns 'written' | 'current' | 'replaced' | 'kept'."""
    if path.is_file():
        old = path.read_text()
        if old == text:
            return "current"
        if not clobber:
            return "kept"
        state = "replaced"
    else:
        state = "written"
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)
    return state


def split_one(run_dir, codes, forced_code=None, dry=False, clobber=True):
    """Split one holo run directory.  Returns 'ok' | 'skipped' | 'failed'."""
    run_dir = Path(run_dir).resolve()
    pid = run_dir.name.upper()
    step2, holo, apo = run_dir / STEP2, run_dir / HOLO, run_dir / APO

    code = forced_code or codes.get(pid)
    if not code:
        print(f"{RED}[NO CODE]  {pid}: not in {LST_NAME} (use -c CODE){RESET}")
        return "failed"
    if not step2.is_file():
        print(f"{RED}[MISSING]  {pid}: no {STEP2} -- did step2 finish?{RESET}")
        return "failed"

    holo_text = step2.read_text()
    lines = holo_text.splitlines()
    hit = [ln for ln in lines if ln[RESNAME] == code]
    keep = [ln for ln in lines if ln[RESNAME] != code]

    if not hit:
        print(f"{RED}[EMPTY]    {pid}: {STEP2} holds no {code} lines -- {APO} would "
              f"equal {HOLO} and apo would silently equal holo.  Refusing.{RESET}")
        return "failed"

    apo_text = "\n".join(keep) + "\n"

    if dry:
        print(f"{GREEN}[WOULD]{RESET}    {pid}  {code}: {HOLO} = {len(lines)} lines, "
              f"{APO} = {len(keep)} lines ({len(hit)} removed)")
        return "ok"

    s_holo = write_variant(holo, holo_text, HOLO, clobber)
    s_apo = write_variant(apo, apo_text, APO, clobber)

    if "kept" in (s_holo, s_apo):
        print(f"{YELLOW}[KEPT]     {pid}: existing variants kept (--no-clobber); they "
              f"may describe an earlier step2 run{RESET}")
        return "skipped"

    # Verify on disk, not in memory: what step3 will eventually read.
    back_holo, back_apo = holo.read_text(), apo.read_text()
    bad = []
    if back_holo != holo_text:
        bad.append(HOLO)
    if back_apo != apo_text:
        bad.append(APO)
    if sum(1 for ln in back_apo.splitlines() if ln[RESNAME] == code):
        bad.append(f"{APO} still holds {code} lines")
    if bad:
        print(f"{RED}[FATAL]    {pid}: {', '.join(bad)} -- removing both variants so "
              f"nothing downstream can use them{RESET}")
        holo.unlink(missing_ok=True)
        apo.unlink(missing_ok=True)
        return "failed"

    tag = {"written": "[SPLIT]  ", "current": "[CURRENT]", "replaced": "[RESPLIT]"}
    state = "replaced" if "replaced" in (s_holo, s_apo) else \
            ("current" if s_holo == s_apo == "current" else "written")
    print(f"{GREEN}{tag[state]}{RESET} {pid}  {code}: "
          f"{HOLO} {len(lines)} lines (sha {sha256_text(holo_text)[:12]}), "
          f"{APO} {len(keep)} lines (sha {sha256_text(apo_text)[:12]}), "
          f"{len(hit)} removed")

    het = sorted({ln[RESNAME] for ln in keep
                  if ln[:6] == "HETATM" and ln[RESNAME].strip()})
    if het:
        print(f"           kept in apo: {', '.join(het)}")
    return "ok"


def main():
    ap = argparse.ArgumentParser(
        description="stepB for run_holo: split step2_out.pdb into holo_step2_out.pdb "
                    "and apo_step2_out.pdb.")
    ap.add_argument("-d", "--dir",
                    help="batch mode: every <PDBID>/ below this directory "
                         "(default: the current directory only, as stepB)")
    ap.add_argument("-l", "--list", help=f"{LST_NAME} (default: searched upward)")
    ap.add_argument("-c", "--code", help="inhibitor code, overriding the .lst lookup")
    ap.add_argument("--no-clobber", action="store_true",
                    help="keep existing variants even when step2_out.pdb has changed")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    base = Path(args.dir).resolve() if args.dir else Path.cwd()
    if not base.is_dir():
        sys.exit(f"ERROR: {base} is not a directory.")

    lst = find_list(base, args.list)
    codes = load_codes(lst) if lst else {}
    if not codes and not args.code:
        sys.exit(f"ERROR: {LST_NAME} not found at or above {base} (use -l or -c).")

    print(f"{CYAN}make_holo_apo_step2_out.py  (stepB){RESET}")
    print(f"  dir  : {base}")
    print(f"  split: {STEP2} -> {HOLO} + {APO}")
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

    tally = {"ok": 0, "skipped": 0, "failed": 0}
    for t in targets:
        tally[split_one(t, codes, args.code, args.dry_run,
                        not args.no_clobber)] += 1

    print(f"\n{CYAN}{'='*60}{RESET}")
    print(f"split: {tally['ok']}   skipped: {tally['skipped']}   "
          f"failed: {tally['failed']}")
    sys.exit(1 if tally["failed"] else 0)


if __name__ == "__main__":
    main()
