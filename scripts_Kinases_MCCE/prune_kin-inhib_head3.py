#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: prune_kin-inhib_head3.py
Edits head3.lst between MCCE step3 and step4 (the "stepC" hook in submit_mcce4.sh).

PURPOSE:
========
Two independent edits to head3.lst:

  1. KINASE INHIBITORS (the 18 ligand codes in INHIBITORS below)
       a. vdw0 -> 0.000 for every inhibitor conformer.
       b. Only the FIRST conformer of each conformer type stays free (FL = "f").
          Every later conformer of that same conftype is turned off (FL = "t",
          occ = 0.00), so MC still chooses among ionization states but sees only
          one rotamer per state.

       Example (IRE, chain A, residue 0001):
            IRE+aA0001_001  ->  f      first of conftype "+a"  -> free
            IRE+aA0001_002  ->  t      later  of conftype "+a"  -> off
            IRE+bA0001_003  ->  f      first of conftype "+b"  -> free
            IRE+bA0001_004  ->  t      later  of conftype "+b"  -> off

       Conformers are grouped by residue + conftype + chain + resSeq, so a
       structure carrying two copies of the same ligand (e.g. 3ZOS has 0LI at
       A1000 and A1004) keeps the first conformer of each conftype per copy.

  2. ARGININE
       ARG+ stays free (FL = "f"); all neutral ARG conformers (crg == 0 and
       nH == 0, i.e. ARG01/ARG02/ARG03/...) are turned off (FL = "t", occ = 0.00).
       With every neutral state off, MC can only ever produce ARG+, while the
       ARG+ rotamers remain free to be sampled. ARG+ is deliberately NOT pinned
       to occ 1.00: a residue usually has several ARG+ rotamers and fixing each
       at occupancy 1.00 would be invalid.

Everything not matched above is left byte-for-byte untouched.

BACKUP:
=======
head3.lst is copied to head3.lst_BK before any edit. If head3.lst_BK already
exists it is NOT overwritten, so the pristine step3 output survives a rerun of
this script. Use --force-backup to refresh it deliberately.

USAGE:
======
  Via submit_mcce4.sh (no arguments; cwd is the run directory):
      stepC="t"
      STEPC="<repo>/scripts_Kinases_MCCE/prune_kin-inhib_head3.py"
  The driver runs it as `$PYEX $STEPC > stepC.log`.

  Manually:
      ./prune_kin-inhib_head3.py                    # edit ./head3.lst
      ./prune_kin-inhib_head3.py -f 3UG2/head3.lst  # edit a specific file
      ./prune_kin-inhib_head3.py --dry-run          # report only, write nothing
"""

import argparse
import shutil
import sys
from pathlib import Path

# Kinase inhibitor ligand codes (matches COFACTORS in ../rm_cofs.sh)
INHIBITORS = {
    "0LI", "0WN", "4MK", "AXI", "B49", "BAX", "DB8", "EMH", "EUI",
    "FMM", "IRE", "LEV", "LQQ", "MI1", "NIL", "STI", "VGH", "YY3",
}

# Fixed-width column layout of head3.lst, as (start, end) half-open slices.
COL = {
    "iConf":   (0, 5),
    "conf":    (6, 20),
    "fl":      (21, 22),
    "occ":     (23, 27),
    "crg":     (28, 34),
    "Em0":     (35, 40),
    "pKa0":    (41, 46),
    "ne":      (47, 49),
    "nH":      (50, 52),
    "vdw0":    (53, 60),
    "vdw1":    (61, 68),
    "tors":    (69, 76),
    "epol":    (77, 84),
    "dsolv":   (85, 92),
    "extra":   (93, 100),
    "history": (101, 111),
}

EXPECTED_HEADER_FIELDS = ["iConf", "CONFORMER", "FL", "occ", "crg", "vdw0"]


def get(line, field):
    a, b = COL[field]
    return line[a:b]


def put(line, field, value):
    """Right-justify value into a fixed-width column without disturbing the rest."""
    a, b = COL[field]
    width = b - a
    text = str(value).rjust(width)
    if len(text) > width:
        raise ValueError(f"value {value!r} too wide for column {field} ({width} chars)")
    return line[:a] + text + line[b:]


def check_header(header, path):
    missing = [f for f in EXPECTED_HEADER_FIELDS if f not in header]
    if missing:
        sys.exit(f"ERROR: {path} header is missing {missing}; is this really a head3.lst?\n"
                 f"       header: {header!r}")


def parse_conformer(line, lineno, path):
    """CONFORMER is RES[0:3] CONFTYPE[3:5] CHAIN[5] RESSEQ[6:10] '_' SEQ[11:14]."""
    conf = get(line, "conf").strip()
    if len(conf) != 14 or conf[10] != "_":
        sys.exit(f"ERROR: {path}:{lineno} unexpected CONFORMER format {conf!r} "
                 f"(expected 14 chars like 'IRE+aA0001_001')")
    return conf[0:3], conf[3:5], conf[0:10]   # resname, conftype, group key


def as_float(line, field, lineno, path):
    text = get(line, field)
    try:
        return float(text)
    except ValueError:
        sys.exit(f"ERROR: {path}:{lineno} column {field} is not a number: {text!r}")


def modify(lines, path):
    """Apply the inhibitor and ARG edits. Returns (new_lines, report)."""
    header, body = lines[0], lines[1:]
    check_header(header, path)

    seen_groups = set()
    out = [header]
    report = {"inh_vdw0": 0, "inh_free": 0, "inh_off": 0,
              "arg_free": 0, "arg_off": 0, "residues": {}}

    for i, line in enumerate(body, start=2):
        if not line.strip():
            out.append(line)
            continue

        resname, conftype, group = parse_conformer(line, i, path)
        tally = report["residues"].setdefault(resname, {"free": 0, "off": 0})

        if resname in INHIBITORS:
            # (a) zero vdw0 on every inhibitor conformer
            line = put(line, "vdw0", f"{0.0:7.3f}")
            report["inh_vdw0"] += 1

            # (b) first conformer of each conftype stays free, the rest go off
            if group in seen_groups:
                line = put(line, "fl", "t")
                line = put(line, "occ", f"{0.0:4.2f}")
                report["inh_off"] += 1
                tally["off"] += 1
            else:
                seen_groups.add(group)
                line = put(line, "fl", "f")
                report["inh_free"] += 1
                tally["free"] += 1

        elif resname == "ARG":
            crg = as_float(line, "crg", i, path)
            nH = as_float(line, "nH", i, path)
            if abs(crg) < 1e-6 and abs(nH) < 1e-6:      # neutral ARG -> off
                line = put(line, "fl", "t")
                line = put(line, "occ", f"{0.0:4.2f}")
                report["arg_off"] += 1
                tally["off"] += 1
            else:                                        # ARG+ -> free
                line = put(line, "fl", "f")
                report["arg_free"] += 1
                tally["free"] += 1

        out.append(line)

    return out, report


def main():
    ap = argparse.ArgumentParser(
        description="stepC: edit head3.lst between MCCE step3 and step4.")
    ap.add_argument("-f", "--file", default="head3.lst",
                    help="head3.lst to edit (default: ./head3.lst)")
    ap.add_argument("--backup-suffix", default="_BK",
                    help="suffix for the backup copy (default: _BK)")
    ap.add_argument("--force-backup", action="store_true",
                    help="overwrite an existing backup instead of keeping it")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change; write nothing")
    args = ap.parse_args()

    path = Path(args.file)
    if not path.is_file():
        sys.exit(f"ERROR: {path} not found (stepC runs in the MCCE run directory).")

    lines = path.read_text().splitlines()
    if len(lines) < 2:
        sys.exit(f"ERROR: {path} has no conformer lines.")

    new_lines, rep = modify(lines, path)

    print(f"prune_kin-inhib_head3.py  ->  {path.resolve()}")
    print(f"  conformers read                  : {len(lines) - 1}")
    print(f"  inhibitor conformers, vdw0 -> 0  : {rep['inh_vdw0']}")
    print(f"  inhibitor conformers free  (FL=f): {rep['inh_free']}")
    print(f"  inhibitor conformers off   (FL=t): {rep['inh_off']}")
    print(f"  ARG+ conformers free       (FL=f): {rep['arg_free']}")
    print(f"  ARG neutral conformers off (FL=t): {rep['arg_off']}")

    touched = {r: v for r, v in rep["residues"].items() if v["free"] or v["off"]}
    if touched:
        print("  per residue (free/off):")
        for r in sorted(touched):
            print(f"    {r:>4s}  {touched[r]['free']:>3d} free  {touched[r]['off']:>3d} off")
    else:
        print("  NOTE: no inhibitor or ARG conformers found; head3.lst unchanged.")

    if args.dry_run:
        print("  --dry-run: nothing written.")
        return

    backup = path.with_name(path.name + args.backup_suffix)
    if backup.exists() and not args.force_backup:
        print(f"  backup   : {backup.name} already exists, kept (original preserved)")
    else:
        shutil.copy2(path, backup)
        print(f"  backup   : {backup.name} written")

    path.write_text("\n".join(new_lines) + "\n")
    print(f"  written  : {path.name}")


if __name__ == "__main__":
    main()
