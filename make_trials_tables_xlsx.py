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
  ./make_trials_tables_xlsx.py --out mytables.xlsx
  ./make_trials_tables_xlsx.py --root /path/to/Kinases_MCCE

NOTE: openpyxl writes formulas without cached values, so the computed columns
read as blank until the file is opened in Excel (or recalculated with
LibreOffice), which fills them in.  --values-copy also writes a static copy
with the numbers baked in, for previewing.
"""

import argparse
import glob
import os
import re
import statistics as st
import sys
from datetime import datetime

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import column_index_from_string, get_column_letter

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

    group = ["", "", "", "Ligand", "Ligand", "Ligand", "Ligand", "Ligand", "Ligand",
             "Apo-kinase", "Apo-kinase", "Holo-kinase", "Holo-kinase", "Holo-Apo", "Holo-Apo"]
    head = ["PDBID", "Ligand", "Kinase", "#conf", "crg soln", "± SEM", "crg bound", "± SEM",
            "∆crg", "± SEM", "charge", "± SEM", "charge", "± SEM", "∆crg", "± SEM"]
    ws.append([]); ws.append(group); ws.append(head)
    for c in range(1, len(head) + 1):
        for r in (4, 5):
            ws.cell(r, c).font = F_BOLD
            ws.cell(r, c).alignment = Alignment(horizontal="center")
        ws.cell(5, c).border = UNDER

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
        for r in sorted(by_kinase[kinase], key=lambda x: x["pdb"]):
            p = r["pdb"]
            ws.cell(i, 1).value = p
            ws.cell(i, 2).value = r["inhibitor"]
            ws.cell(i, 3).value = r["kinase"]
            ws.cell(i, 4).value = f"=AVERAGE({rng('nconf', p)})"
            for col, key in ((5, "crg_soln"), (7, "crg_bound"),
                             (11, "apo_net"), (13, "holo_net")):
                ws.cell(i, col).value = f"=AVERAGE({rng(key, p)})"
                ws.cell(i, col + 1).value = (
                    f"=IF(COUNT({rng(key, p)})>1,"
                    f"STDEV({rng(key, p)})/SQRT(COUNT({rng(key, p)})),0)")
            ws.cell(i, 9).value = f"=G{i}-E{i}"
            ws.cell(i, 10).value = f"=SQRT(F{i}^2+H{i}^2)"
            ws.cell(i, 15).value = f"=M{i}-K{i}"
            ws.cell(i, 16).value = f"=SQRT(L{i}^2+N{i}^2)"
            for c in range(1, len(head) + 1):
                ws.cell(i, c).font = F_BODY
                if c == 4:
                    ws.cell(i, c).number_format = "0.0"
                elif c in (6, 8, 10, 12, 14, 16):
                    ws.cell(i, c).number_format = "0.000"
                elif c >= 5:
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
    ws.freeze_panes = "A6"


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
    ap.add_argument("--out", default="kinase_project-trials-tables.xlsx",
                    help="output workbook (default: %(default)s)")
    ap.add_argument("--values-copy", action="store_true",
                    help="also write a *_values.xlsx with the numbers baked in, for previewing "
                         "without Excel")
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

    out = os.path.join(root, args.out) if not os.path.isabs(args.out) else args.out
    wb.save(out)
    print(f"{GREEN}  wrote {out}{RESET}")
    print(f"  sheets: {', '.join(wb.sheetnames)}")

    if args.values_copy:
        vals = out.replace(".xlsx", "_values.xlsx")
        vb = Workbook(); vb.remove(vb.active)
        vs = vb.create_sheet("Table 1 (Trials) values")
        vs.append(["PDBID", "Ligand", "Kinase", "#conf",
                   "crg soln", "sem", "crg bound", "sem", "dcrg", "sem",
                   "apo charge", "sem", "holo charge", "sem", "holo-apo", "sem"])
        for c in range(1, 17):
            vs.cell(1, c).font = F_BOLD
        for r in rows:
            def ms(key):
                v = [r["per"][t][key] for t in trials if r["per"][t][key] is not None]
                if not v:
                    return None, None
                return st.fmean(v), (st.stdev(v) / len(v) ** 0.5 if len(v) > 1 else 0.0)
            cs, cse = ms("crg_soln"); cb, cbe = ms("crg_bound")
            an, ane = ms("apo_net"); hn, hne = ms("holo_net")
            nc, _ = ms("nconf")
            d = (cb - cs) if (cb is not None and cs is not None) else None
            de = ((cse ** 2 + cbe ** 2) ** 0.5) if (cse is not None and cbe is not None) else None
            ha = (hn - an) if (hn is not None and an is not None) else None
            hae = ((ane ** 2 + hne ** 2) ** 0.5) if (ane is not None and hne is not None) else None
            vs.append([r["pdb"], r["inhibitor"], r["kinase"], nc,
                       cs, cse, cb, cbe, d, de, an, ane, hn, hne, ha, hae])
        for rr in range(2, vs.max_row + 1):
            for c in range(1, 17):
                vs.cell(rr, c).font = F_BODY
                if c >= 5:
                    vs.cell(rr, c).number_format = "0.000" if c % 2 == 0 else "0.00"
        vb.save(vals)
        print(f"{GREEN}  wrote {vals}{RESET}  (static values, no formulas)")


if __name__ == "__main__":
    main()
