#!/usr/bin/env python3
"""
Created on Sep 22 2026
@author: Gehan Ranepura

Name: make_trials_tables_xlsx.py
Builds kinase_project-trials-tables.xlsx from whatever trials are on disk --
the Trials counterpart of kinase_project-final-tables.xlsx.

It discovers Trial*/ directories itself, so adding Trial04, Trial05 ... and
re-running is all that is needed: the per-trial columns, the AVERAGE/STDEV
ranges, n, and the seed table all widen automatically.

SHEETS:
=======
  Per-Trial Data      the raw numbers, one row per PDB, one column per trial.
                      The only hardcoded values in the workbook.
  Table 1 (Trials)    the headline table, laid out like Table 1 of
                      kinase_project-final-tables.xlsx (kinase section rows,
                      Arial 12) but with a +/- SEM column beside every value.
                      Every cell is a FORMULA over Per-Trial Data, so editing
                      or adding raw values updates the table.
  Methods             how the statistics were computed, n, the seeds, and the
                      caveats -- written from the run, not by hand.
  Outliers            the cross-trial holo-apo outlier table, if it exists.
  Reproducibility     identical / scattered / two-state breakdown, if it exists.

WHERE THE NUMBERS COME FROM:
============================
  charges      <trial>/run_{inhib,holo,apo}/<PDB>/xts_sum_crg.out  (pH 7.4,
               entropy-corrected; step4 does NOT produce this -- 1-run_xts_corr.py
               does, so run that in every trial first)
  ligand #conf <trial>/run_inhib/<PDB>/head3.lst
  kinase       the Kinase column of pdb_inhibitor.lst
  seeds        MONTE_SEED in each trial's submit script
  outliers     plots_Trials_Fig4A_holo_vs_apo_xts/ and plots_Trials_Fig3_.../
               (produced by the two plot_trials_*.py scripts)

USAGE:
======
  ./make_trials_tables_xlsx.py                  # every Trial*/ here
  ./make_trials_tables_xlsx.py --glob 'Trial0[12]'
  ./make_trials_tables_xlsx.py --outdir tables_v2
  ./make_trials_tables_xlsx.py --root /path/to/Kinases_MCCE

The workbook is written into tables_Trials/ (--outdir), next to the
plots_Trials_* directories the two plot_trials_*.py scripts produce.

NOTE: openpyxl writes formulas with no cached value, so the computed columns
read as blank until Excel opens the file and calculates them.
"""

import argparse
import glob
import math
import os
import re
import shutil
import statistics as st
import sys
from datetime import datetime

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import column_index_from_string, get_column_letter

# Excel's Number format with 3 decimal places.  Display only: the cell keeps the
# full-precision value, so a whole number shows as "1.000" and nothing is
# rounded in the data.  ("0.###" would drop the padding but Excel prints its
# decimal point regardless, giving "1." on whole numbers.)
UPTO3 = "0.000"


def r3(formula):
    """Formulas are stored unrounded; the number format does the rounding."""
    return formula if formula.startswith("=") else f"={formula}"
# Structural / tautomer annotations carried over from Table 1 of
# kinase_project-final-tables.xlsx.  They are manual assignments that MCCE does
# not output, so they are hardcoded here and the workbook no longer needs that
# file.  (PDB: (delta-Taut, DFG, kinase conformation))
ANNOTATIONS = {
    "1XKK": ("N/A", "In", "DFG-in (active)"),
    "2EUF": ("N/A", "In", "DFG-in (active)"),
    "2HYY": ("Yes", "Out", "DFG-out (inactive)"),
    "2ITO": ("No", "In", "DFG-in (active)"),
    "2ITY": ("No", "In", "DFG-in (active)"),
    "2ITZ": ("No", "In", "DFG-in (active)"),
    "2WGJ": ("No", "In", "DFG-in (active)"),
    "2XP2": ("No", "In", "DFG-in (active)"),
    "2YFX": ("No", "In", "DFG-in (active)"),
    "3AOX": ("N/A", "In", "DFG-in (active)"),
    "3CS9": ("No", "Out", "DFG-out (inactive)"),
    "3LXK": ("N/A", "In", "DFG-in (active)"),
    "3PYY": ("Yes", "Out", "DFG-out (inactive)"),
    "3UE4": ("Yes", "In", "DFG-in (active)"),
    "3UG2": ("No", "In", "DFG-in (active)"),
    "3WZD": ("No", "In", "DFG-in (active)"),
    "3WZE": ("No", "Out", "DFG-out (inactive)"),
    "3ZOS": ("Yes", "Out", "DFG-out (inactive)"),
    "4AG8": ("Yes", "Out", "DFG-out (inactive)"),
    "4AGC": ("Yes", "Out", "DFG-out (inactive)"),
    "4AGD": ("No", "Out", "DFG-out (inactive)"),
    "4AN2": ("No", "In", "DFG-in (active)"),
    "4ANQ": ("No", "In", "DFG-in (active)"),
    "4ANS": ("No", "In", "DFG-in (active)"),
    "4ASD": ("No", "Out", "DFG-out (inactive)"),
    "4G5J": ("N/A", "In", "DFG-in (active)"),
    "4G5P": ("N/A", "In", "DFG-in (active)"),
    "4I22": ("No", "In", "DFG-in (active)"),
    "4LMN": ("No", "In", "DFG-in (active)"),
    "4MKC": ("N/A", "In", "DFG-in (active)"),
    "4WKQ": ("No", "In", "DFG-in (active)"),
    "4ZAU": ("No", "In", "DFG-in (active)"),
    "5AAA": ("No", "In", "DFG-in (active)"),
    "5AAB": ("No", "In", "DFG-in (active)"),
    # 5AAC is absent from that Table 1; taken from its 5AAA/5AAB series
    "5AAC": ("No", "In", "DFG-in (active)"),
    "5L2I": ("N/A", "In", "DFG-in (active)"),
    "5MO4": ("No", "Out", "DFG-out (inactive)"),
}


# Conformer energies (kcal/mol).  MCCE's own extra.tpl is the live source --
# "EXTRA  DB8+1  0.071" -- and read_extra_energies() below prefers it, so a new
# ligand is picked up without touching this file.  The table here is the frozen
# copy from SI.3.Conf of kinase_project-final-tables.xlsx, used when extra.tpl
# cannot be found; all 63 entries were verified identical to it.
#
# Not recomputed: energy is -RT*ln(P(i) soln), and xts_fort.38 stores
# occupancies to 3 decimals, so a conformer whose population rounds to 0.000
# would have no defined energy (DB802 is one).
CONF_ENERGY = {
    "0LI+1": 0.316,
    "0LI+2": 0.573,
    "0LI01": 2.03,
    "0WN+1": 0.005,
    "0WN+a": 4.147,
    "0WN01": 2.861,
    "4MK+1": 0.0,
    "4MK+a": 5.439,
    "4MK01": 3.967,
    "AXI+1": 3.79,
    "AXI01": 0.124,
    "AXI02": 0.993,
    "B49+1": 0.008,
    "B4901": 2.573,
    "B4902": 5.676,
    "B4903": 5.692,
    "BAX+1": 3.207,
    "BAX-1": 5.84,
    "BAX01": 0.004,
    "BAX02": 3.62,
    "DB8+1": 0.071,
    "DB8+2": 1.695,
    "DB801": 1.715,
    "DB802": 5.065,
    "EMH+1": 0.342,
    "EMH01": 0.488,
    "EUI+1": 0.008,
    "EUI01": 2.536,
    "EUI02": 5.637,
    "FMM+1": 0.031,
    "FMM+a": 5.085,
    "FMM01": 1.766,
    "IRE+1": 0.717,
    "IRE+2": 3.728,
    "IRE+3": 4.395,
    "IRE+a": 4.412,
    "IRE01": 0.212,
    "LEV-1": 5.763,
    "LEV-2": 5.763,
    "LEV01": 0.0,
    "LQQ+1": 0.256,
    "LQQ01": 2.434,
    "MI1+1": 2.096,
    "MI101": 0.018,
    "NIL+1": 1.404,
    "NIL+2": 5.173,
    "NIL01": 0.058,
    "STI+1": 0.316,
    "STI+2": 0.573,
    "STI+a": 5.43,
    "STI+b": 5.687,
    "STI01": 2.031,
    "VGH+1": 0.585,
    "VGH+2": 3.587,
    "VGH+a": 0.283,
    "VGH01": 3.192,
    "YY3+1": 0.016,
    "YY3+2": 3.31,
    "YY3+3": 3.893,
    "YY3+4": 4.431,
    "YY3+a": 4.972,
    "YY3+b": 5.449,
    "YY301": 2.299,
}

RT_KCAL = 0.5925          # kcal/mol, as used in SI.3.Conf of the published workbook
TREES = {"inhib": "run_inhib", "holo": "run_holo", "apo": "run_apo"}
LEGACY = {"run_holo": "run_kin", "run_apo": "run_prot2", "run_inhib": "run_cof2"}
PDBID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

GREEN, YELLOW, CYAN, RED, RESET = (
    ("\033[0;32m", "\033[1;33m", "\033[0;36m", "\033[0;31m", "\033[0m")
    if sys.stdout.isatty() else ("", "", "", "", "")
)

ARIAL = "Arial"
F_BOLD = Font(name=ARIAL, size=12, bold=True)
F_BODY = Font(name=ARIAL, size=12)
F_NOTE = Font(name=ARIAL, size=11, italic=True)
UNDER = Border(bottom=Side(style="thin", color="999999"))
GREY = PatternFill("solid", fgColor="EEEEEE")


# ----------------------------------------------------------------- gathering
def tree_dir(trial, which):
    """run_holo etc., falling back to the legacy name (run_kin ...)."""
    name = TREES[which]
    for n in (name, LEGACY[name]):
        p = os.path.join(trial, n)
        if os.path.isdir(p):
            return p
    return os.path.join(trial, name)


def read_sum_crg(path):
    """xts_sum_crg.out -> {residue: charge}.  Net_Charge is one of the keys."""
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        f = line.split()
        if len(f) == 2 and not line.lower().startswith(" ph"):
            try:
                out[f[0]] = float(f[1])
            except ValueError:
                pass
    return out


def count_conformers(path, code):
    """Conformers of one residue type in head3.lst (CONFORMER starts at col 7)."""
    if not os.path.isfile(path):
        return None
    return sum(1 for line in open(path) if line[6:9] == code)


def load_manifest(path):
    """pdb_inhibitor.lst -> {PDB: (inhibitor, code, kinase)}."""
    out = {}
    for i, line in enumerate(open(path)):
        f = line.split()
        if i == 0 or len(f) < 3:
            continue
        out[f[0].upper()] = (f[1], f[2], f[3] if len(f) > 3 else "")
    if not out:
        sys.exit(f"ERROR: {path} holds no rows.")
    return out


def trial_seed(trial):
    """MONTE_SEED actually used, read from whichever submit script has it."""
    for rel in ("run_holo/submit_mcce4_s3s4.sh", "run_kin/submit_mcce4_s3s4.sh",
                "run_inhib/submit_mcce4.sh", "run_cof2/submit_mcce4.sh"):
        p = os.path.join(trial, rel)
        if os.path.isfile(p):
            for line in open(p):
                if line.startswith("STEP4=") and "MONTE_SEED=" in line:
                    return line.split("MONTE_SEED=")[1].split()[0].rstrip('"')
    return "?"


def gather(root, trials, manifest):
    """[{pdb, inhibitor, kinase, per: {trial: {...}}}] for every PDB in the manifest."""
    rows, incomplete = [], []
    for pdb in sorted(manifest):
        inhibitor, code, kinase = manifest[pdb]
        rec = {"pdb": pdb, "inhibitor": inhibitor, "kinase": kinase, "per": {}}
        for t in trials:
            td = os.path.join(root, t)
            inh = read_sum_crg(os.path.join(tree_dir(td, "inhib"), pdb, "xts_sum_crg.out"))
            holo = read_sum_crg(os.path.join(tree_dir(td, "holo"), pdb, "xts_sum_crg.out"))
            apo = read_sum_crg(os.path.join(tree_dir(td, "apo"), pdb, "xts_sum_crg.out"))
            lig_s = sorted(k for k in inh if k[:3] == code)
            lig_b = sorted(k for k in holo if k[:3] == code)
            rec["per"][t] = {
                "nconf": count_conformers(
                    os.path.join(tree_dir(td, "inhib"), pdb, "head3.lst"), code),
                "crg_soln": inh.get(lig_s[0]) if lig_s else None,
                "crg_bound": holo.get(lig_b[0]) if lig_b else None,
                "apo_net": apo.get("Net_Charge"),
                "holo_net": holo.get("Net_Charge"),
            }
        missing = [t for t in trials
                   if any(rec["per"][t][k] is None
                          for k in ("crg_soln", "crg_bound", "apo_net", "holo_net"))]
        if missing:
            incomplete.append((pdb, missing))
        rows.append(rec)
    return rows, incomplete


def read_fort38(path):
    """xts_fort.38 -> {conformer: occupancy at the single pH}."""
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        f = line.split()
        if len(f) >= 2 and not line.lower().startswith(" ph"):
            try:
                out[f[0]] = float(f[1])
            except ValueError:
                pass
    return out


def read_head3_charges(path):
    """head3.lst -> {conformer: charge}."""
    out = {}
    if not os.path.isfile(path):
        return out
    for i, line in enumerate(open(path)):
        if i == 0 or not line.strip():
            continue
        try:
            out[line[6:20].strip()] = float(line[28:34])
        except ValueError:
            pass
    return out


def read_extra_energies():
    """
    {conformer: energy} from MCCE4's extra.tpl ("EXTRA  DB8+1  0.071").
    Returns {} when the file cannot be found, and the caller falls back to the
    CONF_ENERGY table above.
    """
    cands = []
    if os.environ.get("MCCE_HOME"):
        cands.append(os.path.join(os.environ["MCCE_HOME"], "extra.tpl"))
    mcce = shutil.which("mcce")
    if mcce:
        cands.append(os.path.join(
            os.path.dirname(os.path.dirname(os.path.realpath(mcce))), "extra.tpl"))
    cands.append(os.path.expanduser("~/MCCE4/extra.tpl"))
    for path in cands:
        if not os.path.isfile(path):
            continue
        out = {}
        for line in open(path):
            f = line.split()
            if len(f) >= 3 and f[0] == "EXTRA":
                try:
                    out[f[1]] = float(f[2])
                except ValueError:
                    pass
        if out:
            return out, path
    return {}, None


def find_param_dir():
    """MCCE4's param/ directory, where the ligand .ftpl topologies live."""
    for cand in (os.environ.get("MCCE_HOME"), None):
        if cand:
            d = os.path.join(cand, "param")
            if os.path.isdir(d):
                return d
    mcce = shutil.which("mcce")
    if mcce:
        d = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(mcce))), "param")
        if os.path.isdir(d):
            return d
    d = os.path.expanduser("~/MCCE4/param")
    return d if os.path.isdir(d) else None


def protonated_atoms(param_dir, code):
    """
    {conformer: "NAU, NBI"} -- the polar (N/O) atoms carrying an H in each
    conformer, read from the ligand's .ftpl CONNECT records.  This is the
    "protonated atoms" column of SI.3.Conf, derived rather than transcribed.
    """
    out = {}
    if not param_dir:
        return out
    path = os.path.join(param_dir, f"{code}.ftpl")
    if not os.path.isfile(path):
        return out
    per_conf = {}
    for line in open(path):
        if not line.startswith("CONNECT"):
            continue
        m = re.match(r'CONNECT,\s*"(.{4})",\s*(\S+?):\s*\S*,?(.*)', line.strip())
        if not m:
            continue
        atom, conf, rest = m.group(1).strip(), m.group(2).rstrip(":"), m.group(3)
        if not atom or atom[0] not in "NO":
            continue
        bonded = re.findall(r'"(.{4})"', rest)
        if any(b.strip().startswith("H") for b in bonded):
            per_conf.setdefault(conf, set()).add(atom)
    for conf, atoms in per_conf.items():
        out[conf] = ", ".join(sorted(atoms))
    return out


def conformer_labels(types, charges):
    """
    {conformer: "DB8+1a"} -- the readable label of SI.3.Conf: ligand code, the
    charge, then a letter when a charge has more than one conformer (a single
    one gets no letter, as B49+1 does in the original).
    """
    by_charge = {}
    for ctype in sorted(types):
        q = charges.get(ctype)
        key = 0 if q is None else int(round(q))
        by_charge.setdefault(key, []).append(ctype)
    out = {}
    for q, members in by_charge.items():
        tag = f"+{q}" if q > 0 else (f"{q}" if q < 0 else " 0")
        for i, ctype in enumerate(members):
            suffix = "" if len(members) == 1 else chr(ord("a") + i)
            out[ctype] = f"{ctype[:3]}{tag}{suffix}"
    return out


def gather_conformers(root, trials, manifest):
    """
    Per PDB, the ligand's conformer TYPES (e.g. FMM+1, FMM01) with their
    occupancy in solution and bound, per trial.

    Types rather than individual rotamers: step2 is stochastic, so the number of
    rotamers inside a type varies between trials (FMM01 x2/x3/x2) and rotamers
    cannot be matched one-to-one.  The type -- the protonation/tautomer state --
    is stable, and it is what the original SI.3.Conf's "MCCE Conf Name" column
    lists.  Occupancies of the rotamers in a type are summed within each trial
    before averaging across trials.
    """
    out = {}
    param_dir = find_param_dir()
    for pdb in sorted(manifest):
        inhibitor, code, kinase = manifest[pdb]
        types = {}
        for t in trials:
            td = os.path.join(root, t)
            for which, tree in (("soln", "inhib"), ("bound", "holo")):
                d = os.path.join(tree_dir(td, tree), pdb)
                occ = read_fort38(os.path.join(d, "xts_fort.38"))
                crg = read_head3_charges(os.path.join(d, "head3.lst"))
                lig = sorted(k for k in occ if k[:3] == code)
                if not lig:
                    continue
                resid = lig[0][5:10]          # first copy only, as in Table 1
                for name in lig:
                    if name[5:10] != resid:
                        continue
                    ctype = name[:5]
                    rec = types.setdefault(ctype, {"charge": crg.get(name),
                                                   "soln": {}, "bound": {},
                                                   "nrot": {}})
                    rec[which][t] = rec[which].get(t, 0.0) + occ[name]
                    if which == "soln":
                        rec["nrot"][t] = rec["nrot"].get(t, 0) + 1
                    if rec["charge"] is None:
                        rec["charge"] = crg.get(name)
        if types:
            charges = {c: d["charge"] for c, d in types.items()}
            out[pdb] = {"inhibitor": inhibitor, "kinase": kinase, "types": types,
                        "labels": conformer_labels(types, charges),
                        "protons": protonated_atoms(param_dir, code)}
    return out


# ----------------------------------------------------------------- sheets
QUANTITIES = [("nconf", "Ligand #conf", "0.0"),
              ("crg_soln", "Ligand crg soln", "0.00"),
              ("crg_bound", "Ligand crg bound", "0.00"),
              ("apo_net", "Apo-kinase charge", "0.00"),
              ("holo_net", "Holo-kinase charge", "0.00")]


def sheet_raw(wb, rows, trials):
    ws = wb.create_sheet("Per-Trial Data")
    n = len(trials)
    top = ["", "", ""]
    sub = ["PDBID", "Ligand", "Kinase"]
    for _, title, _fmt in QUANTITIES:
        top += [title] + [""] * (n - 1)
        sub += [t.replace("Trial", "T") for t in trials]
    ws.append(top)
    ws.append(sub)
    for r in rows:
        line = [r["pdb"], r["inhibitor"], r["kinase"]]
        for key, _t, _f in QUANTITIES:
            line += [r["per"][t][key] for t in trials]
        ws.append(line)

    ncol = 3 + len(QUANTITIES) * n
    for c in range(1, ncol + 1):
        ws.cell(1, c).font = F_BOLD
        ws.cell(1, c).alignment = Alignment(horizontal="center")
        ws.cell(2, c).font = F_BOLD
        ws.cell(2, c).border = UNDER
        ws.cell(2, c).alignment = Alignment(horizontal="center")
    for rr in range(3, 3 + len(rows)):
        for c in range(1, ncol + 1):
            ws.cell(rr, c).font = F_BODY
            if c > 3 + n:                       # everything after #conf is a charge
                ws.cell(rr, c).number_format = "0.00"
    for col in range(1, len(QUANTITIES) + 1):   # merge each quantity's group header
        first = 4 + (col - 1) * n
        if n > 1:
            ws.merge_cells(start_row=1, start_column=first, end_row=1, end_column=first + n - 1)
    ws.freeze_panes = "D3"
    for c, w in zip("ABC", (10, 14, 9)):
        ws.column_dimensions[c].width = w
    for c in range(4, ncol + 1):
        ws.column_dimensions[get_column_letter(c)].width = 8
    note = ws.cell(4 + len(rows), 1)
    note.value = ("Source: xts_sum_crg.out (pH 7.4, entropy-corrected by xts_corr.py) in each "
                  "trial's run_inhib / run_holo / run_apo; #conf from run_inhib/head3.lst. "
                  "These are the workbook's only hardcoded values.")
    note.font = F_NOTE
    return ws, n


def sheet_table1(wb, rows, trials, n):
    ws = wb.create_sheet("Table 1 (Trials)")
    ws["A1"] = ("Table 1 (Trials): inhibitor and kinase charge at pH 7.4, "
                f"mean of {n} independent trial{'s' if n != 1 else ''}")
    ws["A1"].font = F_BOLD
    ws["A2"] = ("Each value is the mean over the trials; ± columns are the standard error of "
                "the mean (SD/√n). See the Methods sheet.")
    ws["A2"].font = F_NOTE

    # One group label per column, so each "± SEM" sits under the same group as the
    # value it belongs to.  GROUPS also drives the merges below.
    GROUPS = [("", 3), ("Ligand", 8), ("Apo-kinase", 2),
              ("Holo-kinase", 2), ("Holo-Apo", 2), ("Structure", 2)]
    group = [label for label, span in GROUPS for _ in range(span)]
    head = ["PDBID", "Ligand", "Kinase", "#conf", "crg soln", "± SEM", "crg bound", "± SEM",
            "∆crg", "± SEM", "∆Taut", "charge", "± SEM", "charge", "± SEM",
            "∆crg", "± SEM", "DFG", "Kinase conformation"]
    assert len(group) == len(head), f"group row {len(group)} != header row {len(head)}"
    ws.append([]); ws.append(group); ws.append(head)
    for c in range(1, len(head) + 1):
        for r in (4, 5):
            ws.cell(r, c).font = F_BOLD
            ws.cell(r, c).alignment = Alignment(horizontal="center")
        ws.cell(5, c).border = UNDER
    # merge each group across its columns so the spans are unambiguous
    col = 1
    for label, span in GROUPS:
        if label and span > 1:
            ws.merge_cells(start_row=4, start_column=col, end_row=4, end_column=col + span - 1)
        col += span

    # source column for each quantity on Per-Trial Data
    first_col = {key: 4 + i * n for i, (key, _t, _f) in enumerate(QUANTITIES)}
    src_row = {r["pdb"]: i + 3 for i, r in enumerate(rows)}

    def rng(key, pdb):
        a = get_column_letter(first_col[key])
        b = get_column_letter(first_col[key] + n - 1)
        return f"'Per-Trial Data'!${a}{src_row[pdb]}:${b}{src_row[pdb]}"

    by_kinase = {}
    for r in rows:
        by_kinase.setdefault(r["kinase"] or "(unassigned)", []).append(r)

    i = 6
    for kinase in sorted(by_kinase):
        for c in range(1, len(head) + 1):
            ws.cell(i, c).fill = GREY
        ws.cell(i, 2).value = kinase
        ws.cell(i, 2).font = F_BOLD
        i += 1
        # within a kinase, order by inhibitor then PDB -- as Table 1 of
        # kinase_project-final-tables.xlsx does
        for r in sorted(by_kinase[kinase], key=lambda x: (x["inhibitor"], x["pdb"])):
            p = r["pdb"]
            ws.cell(i, 1).value = p
            ws.cell(i, 2).value = r["inhibitor"]
            ws.cell(i, 3).value = r["kinase"]
            ws.cell(i, 4).value = f"=AVERAGE({rng('nconf', p)})"
            for col, key in ((5, "crg_soln"), (7, "crg_bound"),
                             (12, "apo_net"), (14, "holo_net")):
                ws.cell(i, col).value = f"=AVERAGE({rng(key, p)})"
                ws.cell(i, col + 1).value = (
                    f"=IF(COUNT({rng(key, p)})>1,"
                    f"STDEV({rng(key, p)})/SQRT(COUNT({rng(key, p)})),0)")
            ws.cell(i, 9).value = f"=G{i}-E{i}"
            ws.cell(i, 10).value = f"=SQRT(F{i}^2+H{i}^2)"
            ws.cell(i, 16).value = f"=N{i}-L{i}"
            ws.cell(i, 17).value = f"=SQRT(M{i}^2+O{i}^2)"
            taut, dfg, conf_txt = ANNOTATIONS.get(p, ("", "", ""))
            ws.cell(i, 11).value = taut
            ws.cell(i, 18).value = dfg
            ws.cell(i, 19).value = conf_txt
            for c in range(1, len(head) + 1):
                ws.cell(i, c).font = F_BODY
                if c == 4:
                    ws.cell(i, c).number_format = "0.0"
                elif c in (6, 8, 10, 13, 15, 17):
                    ws.cell(i, c).number_format = "0.000"
                elif c in (5, 7, 9, 12, 14, 16):
                    ws.cell(i, c).number_format = "0.00"
            i += 1
    note = ws.cell(i + 1, 1)
    note.value = ("∆crg columns are differences of the trial means; their ± is the propagated "
                  "SEM √(a²+b²). The residue-level holo−apo difference is instead paired "
                  "within each trial — see Methods.")
    note.font = F_NOTE
    for c, w in zip("ABC", (10, 14, 9)):
        ws.column_dimensions[c].width = w
    for c in range(4, len(head) + 1):
        ws.column_dimensions[get_column_letter(c)].width = 9
    ws.column_dimensions["K"].width = 8       # ∆Taut
    ws.column_dimensions["R"].width = 6       # DFG
    ws.column_dimensions["S"].width = 20      # Kinase conformation
    ws.freeze_panes = "A6"


def sheet_conf_raw(wb, conf, trials):
    """Per-trial occupancies per ligand conformer type -- the raw layer for SI-Table2."""
    ws = wb.create_sheet("Per-Trial Conf")
    n = len(trials)
    top = ["", "", "", ""] + ["P(i) soln"] + [""] * (n - 1) \
          + ["P(i) bound"] + [""] * (n - 1) + ["# rotamers"] + [""] * (n - 1)
    sub = ["PDBID", "Ligand", "Conf type", "charge"] + \
          [t.replace("Trial", "T") for _ in range(3) for t in trials]
    ws.append(top); ws.append(sub)
    index = {}
    r = 3
    for pdb in sorted(conf):
        rec = conf[pdb]
        for ctype in sorted(rec["types"]):
            d = rec["types"][ctype]
            ws.append([pdb, rec["inhibitor"], ctype, d["charge"]]
                      + [d["soln"].get(t) for t in trials]
                      + [d["bound"].get(t) for t in trials]
                      + [d["nrot"].get(t) for t in trials])
            index.setdefault(pdb, []).append((ctype, r))
            r += 1
    ncol = 4 + 3 * n
    for c in range(1, ncol + 1):
        ws.cell(1, c).font = F_BOLD
        ws.cell(1, c).alignment = Alignment(horizontal="center")
        ws.cell(2, c).font = F_BOLD
        ws.cell(2, c).border = UNDER
        ws.cell(2, c).alignment = Alignment(horizontal="center")
    for rr in range(3, r):
        for c in range(1, ncol + 1):
            ws.cell(rr, c).font = F_BODY
            if c == 4:
                ws.cell(rr, c).number_format = "0.000"
            elif 5 <= c <= 4 + 2 * n:
                ws.cell(rr, c).number_format = "0.000"
    for col in range(3):
        first = 5 + col * n
        if n > 1:
            ws.merge_cells(start_row=1, start_column=first, end_row=1, end_column=first + n - 1)
    ws.freeze_panes = "E3"
    for c, w in zip("ABCD", (10, 14, 11, 9)):
        ws.column_dimensions[c].width = w
    for c in range(5, ncol + 1):
        ws.column_dimensions[get_column_letter(c)].width = 8
    note = ws.cell(r + 1, 1)
    note.value = ("Occupancy at pH 7.4 from xts_fort.38, summed over the rotamers of each "
                  "conformer type, per trial. Charge from head3.lst. Types rather than "
                  "individual rotamers because step2 is stochastic: the rotamer count inside "
                  "a type varies between trials, so rotamers cannot be matched one-to-one.")
    note.font = F_NOTE
    return index, n


def sheet_si_table2(wb, conf, index, trials, n):
    """
    SI-Table2: ligand conformer populations, in the column order of SI.3.Conf of
    kinase_project-final-tables.xlsx --

        conf type | charge | energy | Boltzmann Factor | Stat Mech | P(i) soln | P(i) bound

    with the one difference that P(i) soln and P(i) bound are each a mean over
    the trials plus its SEM.  energy / Boltzmann / Stat Mech are derived from the
    mean P(i) soln, exactly as in the original.
    """
    ws = wb.create_sheet("SI-Table2")
    ws["A1"] = "SI-Table2: ligand conformer populations at pH 7.4, mean of the trials"
    ws["A1"].font = F_BOLD
    ws["A2"] = ("One row per conformer type (protonation / tautomer state); the rotamers of a "
                "type are summed within each trial, then averaged. ± columns are SEM (SD/√n).")
    ws["A2"].font = F_NOTE
    ws["A3"] = "RT (kcal/mol)"
    ws["A3"].font = F_BOLD
    ws["B3"] = RT_KCAL
    ws["B3"].font = F_BODY
    ws["B3"].number_format = "0.0000"
    ws["D3"] = ("energy: published values (SI.3.Conf) ; Boltzmann = exp(-energy/RT) ; "
                "Stat Mech = Boltzmann / ΣBoltzmann")
    ws["D3"].font = F_NOTE

    GROUPS = [("", 6), ("from P(i) soln", 3), ("P(i) soln", 2), ("P(i) bound", 2)]
    group = [label for label, span in GROUPS for _ in range(span)]
    head = ["PDBID", "Ligand", "conformer", "protonated atoms", "MCCE Conf Name",
            "charge", "energy", "Boltzmann Factor", "Stat Mech",
            "mean", "± SEM", "mean", "± SEM"]
    assert len(group) == len(head), f"{len(group)} group labels vs {len(head)} columns"
    GROUP_ROW, HEAD_ROW, FIRST_DATA = 5, 6, 7      # row 4 stays blank
    for c, (g, h) in enumerate(zip(group, head), 1):
        ws.cell(GROUP_ROW, c).value = g
        ws.cell(HEAD_ROW, c).value = h
        for rr in (GROUP_ROW, HEAD_ROW):
            ws.cell(rr, c).font = F_BOLD
            ws.cell(rr, c).alignment = Alignment(horizontal="center")
        ws.cell(HEAD_ROW, c).border = UNDER
    col = 1
    for label, span in GROUPS:
        if label and span > 1:
            ws.merge_cells(start_row=GROUP_ROW, start_column=col,
                           end_row=GROUP_ROW, end_column=col + span - 1)
        col += span

    def rng(first_col, src_row):
        a = get_column_letter(first_col)
        b = get_column_letter(first_col + n - 1)
        return f"'Per-Trial Conf'!${a}{src_row}:${b}{src_row}"

    by_kinase = {}
    for pdb in conf:
        by_kinase.setdefault(conf[pdb]["kinase"] or "(unassigned)", []).append(pdb)

    extra, extra_path = read_extra_energies()
    energies = dict(CONF_ENERGY)
    energies.update(extra)                      # extra.tpl wins where both have it
    ws["A4"] = (f"conformer energies: {extra_path}" if extra_path
                else "conformer energies: built-in table (extra.tpl not found)")
    ws["A4"].font = F_NOTE
    unpublished = []
    i = FIRST_DATA
    for kinase in sorted(by_kinase):
        for c in range(1, len(head) + 1):
            ws.cell(i, c).fill = GREY
        ws.cell(i, 2).value = kinase
        ws.cell(i, 2).font = F_BOLD
        i += 1
        # within a kinase, order by inhibitor then PDB, as Table 1 does
        for pdb in sorted(by_kinase[kinase], key=lambda x: (conf[x]["inhibitor"], x)):
            rec = conf[pdb]
            first_row = i
            for ctype, src in index[pdb]:
                ws.cell(i, 1).value = pdb
                ws.cell(i, 2).value = rec["inhibitor"]
                ws.cell(i, 3).value = rec["labels"].get(ctype, ctype)
                ws.cell(i, 4).value = rec["protons"].get(ctype, "")
                ws.cell(i, 5).value = ctype
                ws.cell(i, 6).value = r3(f"'Per-Trial Conf'!$D{src}")
                ws.cell(i, 10).value = r3(f"AVERAGE({rng(5, src)})")
                ws.cell(i, 11).value = r3(
                    f"IF(COUNT({rng(5, src)})>1,"
                    f"STDEV({rng(5, src)})/SQRT(COUNT({rng(5, src)})),0)")
                ws.cell(i, 12).value = r3(f"AVERAGE({rng(5 + n, src)})")
                ws.cell(i, 13).value = r3(
                    f"IF(COUNT({rng(5 + n, src)})>1,"
                    f"STDEV({rng(5 + n, src)})/SQRT(COUNT({rng(5 + n, src)})),0)")
                for c in range(1, len(head) + 1):
                    ws.cell(i, c).font = F_BODY
                    if c >= 6:                       # every numeric column
                        ws.cell(i, c).number_format = UPTO3
                i += 1
            last = i - 1
            # energy / Boltzmann / Stat Mech from the mean solution occupancy.
            # A conformer whose P rounds to 0.000 in xts_fort.38 has no defined
            # energy, so those cells stay blank rather than showing a fake value.
            # energy is written as a number, computed here from the mean
            # P(i) soln; Boltzmann and Stat Mech are live formulas off it, the
            # same arrangement as SI.3.Conf of the published workbook.
            for rr, (ctype, _src) in zip(range(first_row, last + 1), index[pdb]):
                if ctype in energies:
                    ws.cell(rr, 7).value = energies[ctype]
                else:
                    # unpublished conformer: fall back to -RT*ln(P(i) soln)
                    vals = [v for v in (rec["types"][ctype]["soln"].get(tr)
                                        for tr in trials) if v is not None]
                    mean_soln = (sum(vals) / len(vals)) if vals else 0.0
                    if mean_soln > 0:
                        ws.cell(rr, 7).value = -RT_KCAL * math.log(mean_soln)
                        unpublished.append(f"{pdb}/{ctype}")
                    # else: no occupancy and no tabulated energy -- left blank
                ws.cell(rr, 8).value = f'=IF(G{rr}="","",EXP(-G{rr}/$B$3))'
                ws.cell(rr, 9).value = (
                    f'=IF(H{rr}="","",H{rr}/SUM(H${first_row}:H${last}))')
            ws.cell(i, 5).value = "SUM"
            for c in (8, 9, 10, 12):
                L = get_column_letter(c)
                ws.cell(i, c).value = r3(f"SUM({L}{first_row}:{L}{last})")
            ws.cell(i + 1, 5).value = "ensemble charge"
            for c in (9, 10, 12):
                L = get_column_letter(c)
                ws.cell(i + 1, c).value = r3(
                    f"SUMPRODUCT($F{first_row}:$F{last},{L}{first_row}:{L}{last})")
            for rr in (i, i + 1):
                for c in range(1, len(head) + 1):
                    ws.cell(rr, c).font = F_BOLD
                    if c >= 6:
                        ws.cell(rr, c).number_format = UPTO3
            i += 3
    note = ws.cell(i, 1)
    note.value = ("SUM should be 1.000 within rounding. 'ensemble charge' is Σ charge × P — "
                  "computed from Stat Mech, from P(i) soln and from P(i) bound; the P(i) soln "
                  "value reproduces the ligand charge in Table 1 (crg soln) and the P(i) bound "
                  "value reproduces crg bound, so the conformer populations and the residue "
                  "charges are checked against each other.")
    note.font = F_NOTE
    for c, w in zip("ABCDE", (10, 14, 11, 24, 15)):
        ws.column_dimensions[c].width = w
    for c in range(6, len(head) + 1):
        ws.column_dimensions[get_column_letter(c)].width = 10
    ws.column_dimensions["H"].width = 15      # "Boltzmann Factor" header is the widest
    ws.freeze_panes = f"A{FIRST_DATA}"
    if unpublished:
        print(f"  {YELLOW}[NOTE] {len(unpublished)} conformer(s) have no published "
              f"energy; used -RT*ln(P): {', '.join(unpublished[:6])}"
              f"{' ...' if len(unpublished) > 6 else ''}{RESET}")


def sheet_methods(wb, trials, seeds, n, incomplete):
    ws = wb.create_sheet("Methods")

    def put(r, a, b=None, bold=False, note=False):
        ws.cell(r, 1).value = a
        ws.cell(r, 1).font = F_BOLD if bold else (F_NOTE if note else F_BODY)
        if b is not None:
            ws.cell(r, 2).value = b
            ws.cell(r, 2).font = F_BODY
            ws.cell(r, 2).alignment = Alignment(wrap_text=True, vertical="top")

    r = 1
    put(r, "Methods: how the Trials statistics were computed", bold=True); r += 1
    put(r, f"Generated {datetime.now():%Y-%m-%d %H:%M:%S} by make_trials_tables_xlsx.py",
        note=True); r += 2
    put(r, "Replicates", bold=True); r += 1
    put(r, "n", f"{n} independent trial{'s' if n != 1 else ''} of the full MCCE4 pipeline "
                f"(steps 1-4)"); r += 1
    for t in trials:
        put(r, t, f"MONTE_SEED = {seeds[t]}"); r += 1
    put(r, "What varies",
        "step2 rotamer generation does not reproduce between runs, so each trial has its own "
        "conformer set, and step4 Monte Carlo uses a per-trial seed. A trial is therefore the "
        "replicate unit: re-running step4 alone in a finished tree would hold the conformers "
        "fixed and understate the spread."); r += 2
    put(r, "Statistics", bold=True); r += 1
    for k, v in (("mean", "arithmetic mean over the n trials"),
                 ("SD", "sample standard deviation (n-1 denominator)"),
                 ("SEM", "SD / √n — every ± column, and the plotted error bar"),
                 ("∆crg (Table 1)", "difference of the trial means; ± is the propagated "
                                    "SEM √(a²+b²)"),
                 ("∆q (per residue)", "holo − apo computed WITHIN each trial, then averaged. "
                                      "apo reuses holo's conformers inside a trial, so the "
                                      "pairing cancels the conformer set; combining two "
                                      "independent SEMs would overstate the error.")):
        put(r, k, v); r += 1
    r += 1
    put(r, "Caveat", f"Error bars are ±1 SEM, NOT a confidence interval. With n={n} there are "
                     f"{max(n - 1, 0)} degrees of freedom"
                     + (", so a 95% CI would be mean ± 4.30 × SEM." if n == 3 else ".")); r += 1
    put(r, "Scope", "The spread reflects run-to-run reproducibility of the pipeline. It does "
                    "not include force-field, dielectric (ε=4) or crystal-structure "
                    "uncertainty."); r += 2
    put(r, "Provenance", bold=True); r += 1
    put(r, "Charges", "xts_sum_crg.out at pH 7.4 — the entropy-corrected output of xts_corr.py, "
                      "which step4 does not produce on its own"); r += 1
    put(r, "Ligand #conf", "conformer count for the ligand in run_inhib/head3.lst"); r += 1
    put(r, "Kinase", "the Kinase column of pdb_inhibitor.lst"); r += 1
    if incomplete:
        r += 1
        put(r, "Incomplete", bold=True); r += 1
        for pdb, miss in incomplete:
            put(r, pdb, "missing from " + ", ".join(miss)); r += 1
    ws.column_dimensions["A"].width = 20
    ws.column_dimensions["B"].width = 95


def sheet_from_table(wb, title, path, widths, numeric_from=None, note=None):
    """Drop a TSV in as a sheet; skipped silently when the file is absent."""
    if not os.path.isfile(path):
        return False
    ws = wb.create_sheet(title)
    for line in open(path):
        ws.append(line.rstrip("\n").split("\t"))
    for c in range(1, ws.max_column + 1):
        ws.cell(1, c).font = F_BOLD
        ws.cell(1, c).border = UNDER
    for rr in range(2, ws.max_row + 1):
        for c in range(1, ws.max_column + 1):
            cell = ws.cell(rr, c)
            cell.font = F_BODY
            if isinstance(cell.value, str) and cell.value.startswith(("Category", "Key", "#")):
                cell.font = F_BOLD
            try:
                cell.value = float(cell.value)
                if numeric_from and c >= numeric_from:
                    cell.number_format = "0.000"
            except (TypeError, ValueError):
                pass
    for c, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(c)].width = w
    if note:
        cell = ws.cell(ws.max_row + 2, 1)
        cell.value = note
        cell.font = F_NOTE
    ws.freeze_panes = "A2"
    return True


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        description="Build the Trials workbook from whatever Trial*/ directories exist.")
    ap.add_argument("--root", help="directory holding Trial*/ (default: this script's)")
    ap.add_argument("--glob", default="Trial*", help="which trials to include (default: Trial*)")
    ap.add_argument("--outdir", default="tables_Trials",
                    help="directory for the workbook, alongside the plots_Trials_* "
                         "directories (default: %(default)s)")
    ap.add_argument("--out", default="kinase_project-trials-tables.xlsx",
                    help="workbook filename inside --outdir (default: %(default)s)")
    args = ap.parse_args()

    root = os.path.abspath(args.root) if args.root \
        else os.path.dirname(os.path.abspath(__file__))
    trials = [os.path.basename(p) for p in sorted(glob.glob(os.path.join(root, args.glob)))
              if os.path.isdir(p) and os.path.isdir(os.path.join(p, "run_holo"))
              or os.path.isdir(os.path.join(p, "run_kin"))]
    trials = [t for t in trials if os.path.isdir(os.path.join(root, t))]
    if not trials:
        sys.exit(f"ERROR: no trials matched {os.path.join(root, args.glob)}")

    manifest_path = os.path.join(root, "pdb_inhibitor.lst")
    if not os.path.isfile(manifest_path):
        sys.exit(f"ERROR: {manifest_path} not found.")
    manifest = load_manifest(manifest_path)
    seeds = {t: trial_seed(os.path.join(root, t)) for t in trials}

    print(f"{CYAN}make_trials_tables_xlsx.py{RESET}")
    print(f"  root   : {root}")
    print(f"  trials : {', '.join(trials)}  (n = {len(trials)})")
    print(f"  seeds  : {', '.join(f'{t}={seeds[t]}' for t in trials)}")
    print(f"  pdbs   : {len(manifest)} from pdb_inhibitor.lst")
    print()

    rows, incomplete = gather(root, trials, manifest)
    for pdb, miss in incomplete:
        print(f"  {YELLOW}[INCOMPLETE] {pdb}: no xts_sum_crg.out in {', '.join(miss)}{RESET}")
    if incomplete:
        print(f"  {YELLOW}-> run 1-run_xts_corr.py in those trials; the rows are still "
              f"written, with gaps.{RESET}\n")

    wb = Workbook()
    wb.remove(wb.active)
    _, n = sheet_raw(wb, rows, trials)
    sheet_table1(wb, rows, trials, n)
    conf = gather_conformers(root, trials, manifest)
    if conf:
        idx, _ = sheet_conf_raw(wb, conf, trials)
        sheet_si_table2(wb, conf, idx, trials, n)
        print(f"  conformer types: {sum(len(v['types']) for v in conf.values())} "
              f"across {len(conf)} ligands")
    else:
        print(f"  {YELLOW}[SKIP] no xts_fort.38 found -- SI-Table2 omitted{RESET}")
    sheet_methods(wb, trials, seeds, n, incomplete)

    got = sheet_from_table(
        wb, "Outliers (holo-apo)",
        os.path.join(root, "plots_Trials_Fig4A_holo_vs_apo_xts",
                     "Trials_Fig4A_outliers_holo_vs_apo.tsv"),
        [10, 7, 13, 9, 7, 8, 12] + [11] * 7, numeric_from=8,
        note="Reproducible = every trial clears the tier threshold. 'NO' means the "
             "classification depends on which trial you look at.")
    if not got:
        print(f"  {YELLOW}[SKIP] outlier table not found -- run "
              f"plot_trials_comparison_xts_Fig4A.py{RESET}")

    rp = wb.create_sheet("Reproducibility")
    rp["A1"] = "Is the mean a fair summary of the trials?"
    rp["A1"].font = F_BOLD
    rr = 3
    for title, rel in (("Fig4A — per residue (holo vs apo)",
                        "plots_Trials_Fig4A_holo_vs_apo_xts/Trials_Fig4A_reproducibility.tsv"),
                       ("Fig3 — per inhibitor (bound vs solution)",
                        "plots_Trials_Fig3_inhibitors_xts/Trials_Fig3_reproducibility.tsv")):
        path = os.path.join(root, rel)
        if not os.path.isfile(path):
            continue
        rp.cell(rr, 1).value = title
        rp.cell(rr, 1).font = F_BOLD
        rr += 1
        for line in open(path):
            line = line.rstrip("\n")
            if not line:
                rr += 1
                continue
            for c, v in enumerate(line.split("\t"), 1):
                cell = rp.cell(rr, c)
                cell.value = v
                cell.font = F_BOLD if line.startswith(("Category", "Key", "#")) else F_BODY
                try:
                    cell.value = float(v)
                except (TypeError, ValueError):
                    pass
            rr += 1
        rr += 1
    for c, w in zip("ABCD", (34, 20, 10, 62)):
        rp.column_dimensions[c].width = w

    if os.path.isabs(args.out):
        out = args.out                      # an absolute --out wins outright
    else:
        outdir = args.outdir if os.path.isabs(args.outdir) \
            else os.path.join(root, args.outdir)
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, os.path.basename(args.out))
    wb.save(out)
    print(f"{GREEN}  wrote {os.path.relpath(out, root)}{RESET}")
    print(f"  sheets: {', '.join(wb.sheetnames)}")

if __name__ == "__main__":
    main()
