#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: 0-prepare_run_apo.py   (canonical copy: scripts_Kinases_MCCE/prepare_run_apo.py)
Seeds TrialNN/run_apo from TrialNN/run_holo -- step 1 of the apo workflow.

WHAT IT DOES:
=============
For every per-structure directory in run_holo that holo's stepB has split:

  1. delete run_apo/<PDBID> if it already exists, so a re-run never merges a
     new holo with the leftovers of an older one,
  2. copy the WHOLE run_holo/<PDBID> across -- every file steps 1-2 left
     (step1_out.pdb, head1.lst, head2.lst, acc.*, vdw0.lst, rot_stat, run.prm,
     param/, user_param/, the logs) plus BOTH step2 variants that holo's stepB
     made: holo_step2_out.pdb and apo_step2_out.pdb,
  3. link step2_out.pdb -> apo_step2_out.pdb, resetting the link if one is
     already there, so what step3 reads in run_apo is the apo variant and
     nothing else,
  4. write 0-prepare_run_apo_<PDBID>.log recording where the directory came
     from.

Two things are deliberately NOT copied:

  * holo's step2_out.pdb -- in run_holo that file IS holo, and apo must never
    inherit it.  In run_apo the name belongs to a relative symlink pointing at
    apo_step2_out.pdb, which apo's stepB (install_apo_step2_out.py) re-checks
    and, if necessary, resets before step3.
  * holo's step3/step4 products, if holo has already run them (EXCLUDE below).
    head3.lst_BK is never overwritten once it exists, so an inherited copy would
    leave apo holding HOLO's head3.lst as its "pristine step3 output" backup.

WHY A SYMLINK:
==============
step3 reads step2_out.pdb, so in run_apo that name has to resolve to the apo
structure.  A link makes it impossible for the two to drift: there is one apo
structure on disk, apo_step2_out.pdb, and step2_out.pdb is another name for it.
It is relative, so the directory can be moved or copied and still resolve.

NOTE: the submit scripts are NOT copied.  pro_batch symlinks them into each
structure directory (submit_mcce4_s1s2.sh -> ../submit_mcce4_s1s2.sh), so a
copied link would dangle in run_apo, which has no s1s2 script.  pro_batch
creates the one it needs when you launch the apo batch.

Only holo's step2 variants are needed, so this runs as soon as holo's steps 1-2
finish; there is no need to wait for holo's step3/4.

USAGE:
======
It sits at the trial root, next to 1-run_xts_corr.py, so the numbered scripts
read as one sequence.  It works the rest out from there: the source is that
trial's run_holo and the destination is its run_apo.  (A copy left inside
run_apo also works -- it then fills its own directory.)

  cd TrialNN && ./0-prepare_run_apo.py            # seed every structure
                ./0-prepare_run_apo.py --dry-run  # report, write nothing
                ./0-prepare_run_apo.py --keep     # skip existing dirs instead
                                                  #   of replacing them
                ./0-prepare_run_apo.py 1XKK 2ITZ  # just these structures
  From anywhere:  prepare_run_apo.py -t /path/to/Trial01

NEXT:
=====
  cd run_apo && pro_batch kin-pdb -custom submit_mcce4_s3s4.sh \\
                          -job-name apo_s3s4 -j 15
"""

import argparse
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

# Everything steps 1-2 leave is copied.  Excluded are holo's step2_out.pdb
# (apo's stepB makes apo's own) and holo's step3/step4 products, which only
# exist here if holo has already run its s3s4 job.  Globs allowed
# (shutil.ignore_patterns).  Mirrors STEP34_ARTIFACTS in trial_config.sh.
EXCLUDE = (
    "step2_out.pdb",
    "energies", "head3.lst", "head3.lst_BK", "BK_head3.lst",
    "fort.38", "sum_crg.out", "pK.out", "pK_extended.out",
    "xts_fort.38", "xts_sum_crg.out",
    "mc_out", "ms_out", "entropy.out", "entropy_correction.log", "respair.lst",
    "step3.log", "step4.log", "stepC.log", "progress_step3.log",
    # holo's stepB log: apo writes its own, and an inherited one reporting
    # success would be read as apo's
    "stepB.log",
    # pro_batch symlinks these per structure (-> ../submit_mcce4_*.sh); copied
    # into run_apo the link would dangle, and pro_batch remakes it at launch
    "submit_mcce4*.sh",
)

# What apo must and must not have once the copy is done.  step2_out.pdb is not
# copied -- it is created below as a link to apo_step2_out.pdb.
MUST_BE_ABSENT = ("energies", "head3.lst", "head3.lst_BK", "pK.out",
                  "xts_sum_crg.out")
MUST_BE_PRESENT = ("holo_step2_out.pdb", "apo_step2_out.pdb")
STEP2, APO_V = "step2_out.pdb", "apo_step2_out.pdb"
LOGTPL = "0-prepare_run_apo_{pid}.log"   # written in every seeded directory

PDBID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
SRC_NAMES = ("run_holo", "run_kin")   # run_kin is the legacy name for run_holo

GREEN, YELLOW, CYAN, RED, RESET = (
    ("\033[0;32m", "\033[1;33m", "\033[0;36m", "\033[0;31m", "\033[0m")
    if sys.stdout.isatty() else ("", "", "", "", "")
)


def find_trial(explicit):
    """
    Work out (trial directory, default destination) from where this script sits.
    Normally it sits at the trial root, so the trial is its own directory and
    the destination is run_apo.  A copy left inside run_apo also works: the
    trial is then its parent and it fills its own directory.
    """
    if explicit:
        d = Path(explicit).resolve()
        if not d.is_dir():
            sys.exit(f"ERROR: {d} is not a directory.")
        return d, "run_apo"

    here = Path(__file__).resolve().parent
    if any((here / n).is_dir() for n in SRC_NAMES):       # at the trial root
        return here, "run_apo"
    if any((here.parent / n).is_dir() for n in SRC_NAMES):  # inside run_apo
        return here.parent, here.name

    sys.exit(f"ERROR: neither {here} nor {here.parent} holds a "
             f"{' or '.join(SRC_NAMES)}/ -- run this from inside TrialNN/run_apo, "
             f"or pass -t /path/to/TrialNN.")


def link_step2(run_dir):
    """Point step2_out.pdb at apo_step2_out.pdb, resetting whatever is there."""
    link = run_dir / STEP2
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(APO_V)          # relative on purpose: survives a move
    return link


def main():
    ap = argparse.ArgumentParser(
        description="Seed TrialNN/run_apo from TrialNN/run_holo (step3/4 artifacts "
                    "purged).  The inhibitor is deleted later, by stepB.")
    ap.add_argument("pdbids", nargs="*", help="only these structures (default: all)")
    ap.add_argument("-t", "--trial", help="trial directory (default: this script's)")
    ap.add_argument("-s", "--source", help=f"source run dir (default: {SRC_NAMES[0]})")
    ap.add_argument("-o", "--out",
                    help="destination, relative to the trial "
                         "(default: this script's own directory, else run_apo)")
    ap.add_argument("--keep", action="store_true",
                    help="skip structures that already exist in run_apo instead of "
                         "replacing them")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    trial, default_out = find_trial(args.trial)
    out_name = args.out or default_out
    if args.source:
        src_root = (trial / args.source) if not Path(args.source).is_absolute() \
                   else Path(args.source)
    else:
        src_root = next((trial / n for n in SRC_NAMES if (trial / n).is_dir()), None)
    if src_root is None or not src_root.is_dir():
        sys.exit(f"ERROR: no source run directory under {trial} "
                 f"(looked for {', '.join(SRC_NAMES)}).")
    dst_root = trial / out_name

    wanted = {p.upper() for p in args.pdbids}
    sources = sorted(d for d in src_root.iterdir()
                     if d.is_dir() and PDBID_RE.match(d.name)
                     and (not wanted or d.name.upper() in wanted))

    print(f"{CYAN}0-prepare_run_apo.py{RESET}")
    print(f"  trial : {trial}")
    print(f"  from  : {src_root}")
    print(f"  to    : {dst_root}")
    if args.dry_run:
        print(f"  {YELLOW}*** DRY RUN - nothing written ***{RESET}")
    print()

    if not sources:
        sys.exit(f"{RED}ERROR: no <PDBID>/ directories in {src_root}"
                 f"{' matching ' + ', '.join(sorted(wanted)) if wanted else ''} -- "
                 f"has holo been launched?{RESET}")
    missing = wanted - {d.name.upper() for d in sources}
    for m in sorted(missing):
        print(f"{RED}[NO SUCH]  {m}: not in {src_root.name}{RESET}")

    if not args.dry_run:
        dst_root.mkdir(parents=True, exist_ok=True)

    ok = skipped = 0
    failed = len(missing)
    ignore = shutil.ignore_patterns(*EXCLUDE)

    for src in sources:
        pid = src.name
        holo_v, apo_v = src / "holo_step2_out.pdb", src / "apo_step2_out.pdb"
        dst = dst_root / pid

        if not (holo_v.is_file() and apo_v.is_file()):
            have = "step2_out.pdb only" if (src / "step2_out.pdb").is_file() \
                   else "nothing yet"
            print(f"{YELLOW}[WAIT]     {pid}: holo's stepB has not split step2_out.pdb "
                  f"({have}) -- run make_holo_apo_step2_out.py{RESET}")
            skipped += 1
            continue
        if dst.exists() and args.keep:
            print(f"{YELLOW}[SKIP]     {pid}: run_apo/{pid} exists (--keep){RESET}")
            skipped += 1
            continue

        replacing = "  (replacing existing)" if dst.exists() else ""
        route = f"{src_root.name}/{pid} -> {out_name}/{pid}"

        if args.dry_run:
            print(f"{GREEN}[WOULD]{RESET}    {route}{replacing}")
            ok += 1
            continue

        if dst.exists():
            shutil.rmtree(dst)
        try:
            shutil.copytree(src, dst, symlinks=True, ignore=ignore)
        except OSError as e:
            print(f"{RED}[FAILED]   {pid}: {e}{RESET}")
            failed += 1
            continue

        left = [a for a in MUST_BE_ABSENT if (dst / a).exists()]
        if left:
            print(f"{RED}[FATAL]    {pid}: still holds {', '.join(left)} after the "
                  f"purge -- removing {out_name}/{pid}{RESET}")
            shutil.rmtree(dst)
            failed += 1
            continue
        # apo's structure gets the name step3 reads
        link = link_step2(dst)
        if not (link.is_symlink() and link.resolve() == (dst / APO_V).resolve()
                and link.is_file()):
            print(f"{RED}[FATAL]    {pid}: {STEP2} -> {APO_V} does not resolve "
                  f"-- removing {out_name}/{pid}{RESET}")
            shutil.rmtree(dst)
            failed += 1
            continue

        gone = [a for a in MUST_BE_PRESENT if not (dst / a).is_file()]
        if gone:
            print(f"{RED}[FATAL]    {pid}: {', '.join(gone)} did not survive the copy "
                  f"-- removing {out_name}/{pid}{RESET}")
            shutil.rmtree(dst)
            failed += 1
            continue

        (dst / LOGTPL.format(pid=pid)).write_text(
            f"""prepared_by  : 0-prepare_run_apo.py
prepared_at  : {datetime.now():%Y-%m-%d %H:%M:%S}
source_dir   : {src}
step2_out.pdb: a relative symlink to apo_step2_out.pdb -- that is what step3
               reads here.  apo's stepB (install_apo_step2_out.py) re-checks the
               pair and resets the link if needed, inside the step3 job.
note         : apo reuses holo's step2 conformers.  step1/step2 are OFF in
               run_apo/submit_mcce4_s3s4.sh, so the pocket is NOT repacked --
               side chains stay in their holo positions, which is what isolates
               the ligand's electrostatic contribution.
""")

        print(f"{GREEN}[PREPARED]{RESET} {route}{replacing}")
        ok += 1

    print(f"\n{CYAN}{'='*60}{RESET}")
    print(f"prepared: {ok}   skipped: {skipped}   failed: {failed}")
    if ok and not args.dry_run:
        print(f"\nNext: cd {dst_root} && pro_batch kin-pdb -custom "
              f"submit_mcce4_s3s4.sh -job-name apo_s3s4 -j 15")
        print("      (stepB installs each apo_step2_out.pdb as step2_out.pdb, "
              "just before step3)")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
