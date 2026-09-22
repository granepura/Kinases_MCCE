#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: 1-run_xts_corr.py   (canonical copy: scripts_Kinases_MCCE/run_xts_corr.py)
Runs xts_corr.py in every structure directory of a trial -- run_holo, run_apo
and run_inhib -- after step4 has finished.

WHY IT EXISTS:
==============
step4.py --xts only sets MONTE_TSX="t" for the Monte Carlo; it does NOT apply
the conformational entropy correction to the results.  That is a separate
post-processing pass by xts_corr.py, which reads fort.38 (plus head3.lst for the
actual charges) and writes:

    xts_fort.38            entropy-corrected conformer probabilities
    xts_sum_crg.out        the corrected charge summary -- what every published
                           figure reads
    entropy_correction.log the tool's own report of what it changed

The correction applies to non-amino-acid residues, i.e. the inhibitors, so
skipping it in one tree and not another would compare corrected numbers against
uncorrected ones.  This script exists so that cannot happen by accident: it
sweeps all three trees in one pass and reports any structure it could not do.

xts_corr.py writes its outputs into the CURRENT directory under fixed names, so
each structure is run with its own directory as cwd.

WHAT IT SKIPS:
==============
  * a structure with no fort.38 -- step4 has not finished there yet,
  * a structure whose xts_sum_crg.out is already newer than its fort.38, unless
    --force is given.  Re-running is harmless but pointless.

USAGE:
======
  cd TrialNN && ./1-run_xts_corr.py                # all three trees
                ./1-run_xts_corr.py --dry-run      # report, write nothing
                ./1-run_xts_corr.py -t run_inhib   # one tree
                ./1-run_xts_corr.py --force        # redo even if up to date
                ./1-run_xts_corr.py --all          # correct amino acids too
                                                   #   (xts_corr.py --all)
  From anywhere:  run_xts_corr.py --trial /path/to/Trial01

Each structure's console output is kept next to the results as xts_corr.log.
"""

import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

TOOL = "/home/granepura/MCCE4-Tools/mcce4_tools/xts_corr.py"
TREES = ("run_holo", "run_apo", "run_inhib")
LEGACY = {"run_holo": "run_kin", "run_apo": "run_prot2", "run_inhib": "run_cof2"}
SRC, HEAD3 = "fort.38", "head3.lst"
OUTS = ("xts_fort.38", "xts_sum_crg.out")
RUNLOG = "xts_corr.log"
PDBID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

GREEN, YELLOW, CYAN, RED, RESET = (
    ("\033[0;32m", "\033[1;33m", "\033[0;36m", "\033[0;31m", "\033[0m")
    if sys.stdout.isatty() else ("", "", "", "", "")
)


def find_trial(explicit):
    """The trial directory: --trial, else the directory holding this script."""
    if explicit:
        d = Path(explicit).resolve()
        if not d.is_dir():
            sys.exit(f"ERROR: {d} is not a directory.")
        return d
    here = Path(__file__).resolve().parent
    if any((here / t).is_dir() for t in TREES) or \
       any((here / t).is_dir() for t in LEGACY.values()):
        return here
    sys.exit(f"ERROR: {here} holds none of {', '.join(TREES)} -- run this from "
             f"inside a trial directory, or pass --trial /path/to/TrialNN.")


def cause_of(output):
    """The most informative line of a failed xts_corr.py run."""
    lines = [l.rstrip() for l in output.splitlines() if l.strip()]
    if not lines:
        return "no output from xts_corr.py"
    if any("Traceback" in l for l in lines):
        return lines[-1]                      # the exception type and message
    for pat in ("Error:", "ERROR", "No such file", "not found", "Exception"):
        for l in reversed(lines):
            if pat in l:
                return l.strip()
    return lines[-1]


def run_one(run_dir, tool, force, dry, do_all):
    """Correct one structure.  Returns (status, message, cause)."""
    src, out = run_dir / SRC, run_dir / OUTS[1]
    tag = f"{run_dir.parent.name}/{run_dir.name}"

    if not src.is_file():
        return ("waiting",
                f"{YELLOW}[WAIT]    {tag}: no {SRC} -- step4 has not finished here{RESET}",
                None)
    if not (run_dir / HEAD3).is_file():
        c = f"{HEAD3} missing -- step3 output is incomplete in this directory"
        return "failed", f"{RED}[FAILED]  {tag}: {c}{RESET}", c
    if out.is_file() and not force and out.stat().st_mtime >= src.stat().st_mtime:
        return ("current",
                f"{GREEN}[CURRENT] {tag}: {OUTS[1]} already newer than {SRC}{RESET}",
                None)
    if dry:
        verb = "redo" if out.is_file() else "correct"
        return "done", f"{GREEN}[WOULD]   {tag}: {verb}{RESET}", None

    cmd = [sys.executable, str(tool)] + (["--all"] if do_all else [])
    proc = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True)
    (run_dir / RUNLOG).write_text(proc.stdout + proc.stderr)

    if proc.returncode != 0:
        c = f"xts_corr.py exited {proc.returncode}: {cause_of(proc.stdout + proc.stderr)}"
        return ("failed",
                f"{RED}[FAILED]  {tag}{RESET}\n           cause: {c}\n"
                f"           see  : {run_dir}/{RUNLOG}", c)

    missing = [o for o in OUTS if not (run_dir / o).is_file()
               or (run_dir / o).stat().st_size == 0]
    if missing:
        c = (f"{', '.join(missing)} not written or empty despite exit 0: "
             f"{cause_of(proc.stdout + proc.stderr)}")
        return ("failed",
                f"{RED}[FAILED]  {tag}{RESET}\n           cause: {c}\n"
                f"           see  : {run_dir}/{RUNLOG}", c)

    n = sum(1 for _ in out.open())
    return ("done",
            f"{GREEN}[OK]      {tag}: {OUTS[1]} written ({n} lines), "
            f"{OUTS[0]} + entropy_correction.log{RESET}", None)


def main():
    ap = argparse.ArgumentParser(
        description="Run xts_corr.py in every structure directory of a trial.")
    ap.add_argument("-t", "--tree", action="append", dest="trees",
                    help=f"limit to one tree (repeatable); default: {', '.join(TREES)}")
    ap.add_argument("--trial", help="trial directory (default: this script's)")
    ap.add_argument("--tool", default=os.environ.get("XTS_CORR", TOOL),
                    help=f"path to xts_corr.py (default: {TOOL})")
    ap.add_argument("--all", action="store_true",
                    help="pass --all to xts_corr.py (correct amino acids too)")
    ap.add_argument("--force", action="store_true",
                    help="re-run even where xts_sum_crg.out is already up to date")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    ap.add_argument("-j", "--jobs", type=int, default=4, help="parallel workers (default: 4)")
    args = ap.parse_args()

    trial = find_trial(args.trial)
    tool = Path(args.tool).expanduser()
    if not tool.is_file():
        sys.exit(f"ERROR: xts_corr.py not found at {tool} (use --tool or $XTS_CORR).")

    wanted = args.trees or list(TREES)
    trees = []
    for name in wanted:
        if (trial / name).is_dir():
            trees.append(trial / name)
        elif name in LEGACY and (trial / LEGACY[name]).is_dir():
            trees.append(trial / LEGACY[name])       # older run_kin / run_cof2 layout
        else:
            print(f"{YELLOW}[NO TREE] {name} not in {trial}{RESET}")

    print(f"{CYAN}1-run_xts_corr.py{RESET}")
    print(f"  trial : {trial}")
    print(f"  tool  : {tool}")
    print(f"  trees : {', '.join(t.name for t in trees) or '(none)'}")
    print(f"  mode  : {'ALL residues (--all)' if args.all else 'non-amino acids only'}")
    if args.dry_run:
        print(f"  {YELLOW}*** DRY RUN - nothing written ***{RESET}")
    print()
    if not trees:
        sys.exit(1)

    totals = {"done": 0, "current": 0, "waiting": 0, "failed": 0}
    failures = []
    for tree in trees:
        dirs = sorted(d for d in tree.iterdir()
                      if d.is_dir() and PDBID_RE.match(d.name))
        if not dirs:
            print(f"{YELLOW}[EMPTY]   {tree.name}: no <PDBID>/ directories{RESET}\n")
            continue
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            results = list(pool.map(
                lambda d: run_one(d, tool, args.force, args.dry_run, args.all), dirs))
        counts = {}
        for (status, msg, cause), d in zip(results, dirs):
            totals[status] += 1
            counts[status] = counts.get(status, 0) + 1
            print(msg)
            if cause:
                failures.append((f"{d.parent.name}/{d.name}", cause))
        print(f"  {tree.name}: " + "  ".join(f"{k}={v}" for k, v in sorted(counts.items())) + "\n")

    print(f"{CYAN}{'='*60}{RESET}")
    print("  ".join(f"{k}={v}" for k, v in totals.items()))
    if failures:
        print(f"\n{RED}FAILED ({len(failures)}) -- cause per structure:{RESET}")
        for tag, cause in failures:
            print(f"  {tag:22} {cause}")
    if totals["waiting"]:
        print(f"{YELLOW}Structures are still waiting on step4 -- re-run this when they "
              f"finish.{RESET}")
    sys.exit(1 if totals["failed"] else 0)


if __name__ == "__main__":
    main()
