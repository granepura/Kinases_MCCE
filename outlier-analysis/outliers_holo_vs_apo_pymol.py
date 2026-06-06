"""Visualise each outlier residue from ``outliers_holo_vs_apo.tsv`` in PyMOL.

Run inside PyMOL with::

    run outliers_holo_vs_apo_pymol.py

For every row of the TSV the script:

1.  Fetches the relevant PDB entry into its own object.
2.  Shows only the cartoon backbone, removes solvent, and hides every
    other sticks / spheres representation.
3.  Colours the whole backbone ``gray70``.
4.  Aligns the structure onto the first successfully-loaded structure
    using the backbone atoms of the chain listed in the TSV.  The
    very first row becomes the alignment reference for every later one.
5.  Colours the outlier residue orange (tier ``affected``) or red
    (tier ``strong``) and shows its side-chain as sticks.
6.  Duplicates the object so the same residue is also available as
    spheres, giving two display variants per outlier.
"""

import csv
import os

from pymol import cmd


# TSV_PATH = os.path.join(os.path.dirname(__file__), "outliers_holo_vs_apo.tsv")
TSV_PATH = "./outliers_holo_vs_apo.tsv"

TIER_COLOR = {
    "affected": "orange",
    "strong": "red",
}

# Fetch types to try in order.  Some entries aren't served via the plain
# "pdb" format any more (very large assemblies, obsolete entries, etc.),
# so we fall back to mmCIF before giving up.
FETCH_TYPES = ("pdb", "cif")


def _safe(text):
    """Make a string safe to embed in a PyMOL object name."""
    return "".join(ch if ch.isalnum() else "_" for ch in str(text))


def load_outliers(tsv_path=TSV_PATH):
    rows = []
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            # Skip blank / comment lines.
            if not row.get("PDB") or not row.get("PDB").strip():
                continue
            rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def _fetch_structure(pdb, obj_name):
    """Fetch ``pdb`` into ``obj_name``.  Returns True on success."""
    for ftype in FETCH_TYPES:
        try:
            cmd.fetch(pdb, obj_name, async_=0, type=ftype)
        except Exception as exc:
            print(f"[outliers] fetch({pdb}, type={ftype}) raised: {exc}")
            continue
        if obj_name in cmd.get_object_list("(all)"):
            return True
        # Some PyMOL builds create the object under the bare PDB code when
        # the requested name is unusual; detect & rename if that happened.
        if pdb.lower() in cmd.get_object_list("(all)"):
            cmd.set_name(pdb.lower(), obj_name)
            if obj_name in cmd.get_object_list("(all)"):
                return True
    return False


def visualise_outlier(row, reference=None):
    """Render one outlier row.  Returns (obj_sticks, chain) or None."""
    pdb = row["PDB"].upper()
    tier = row["Tier"].lower()
    resname = row["ResName"].upper()
    chain = row["Chain"]
    resnum = row["ResNum"]
    inhibitor = row.get("Inhibitor", "")

    base_color = TIER_COLOR.get(tier, "yellow")

    tag = _safe(f"{pdb}_{inhibitor}_{resname}{resnum}_{tier}")
    obj_sticks = f"{tag}_sticks"
    obj_spheres = f"{tag}_spheres"

    # --- Fetch ---------------------------------------------------------
    if not _fetch_structure(pdb, obj_sticks):
        print(f"[outliers] FAILED to fetch {pdb}; skipping row.")
        return None

    # --- Strip / show cartoon only -------------------------------------
    cmd.remove(f"{obj_sticks} and (solvent or resn HOH)")
    cmd.hide("everything", obj_sticks)
    cmd.show("cartoon", obj_sticks)
    cmd.color("gray70", obj_sticks)

    # --- Align to reference on the TSV-specified chain -----------------
    if reference is not None:
        ref_obj, ref_chain = reference
        if ref_obj != obj_sticks and ref_obj in cmd.get_object_list("(all)"):
            mobile_bb = (
                f"{obj_sticks} and chain {chain} "
                f"and polymer.protein and name N+CA+C+O"
            )
            target_bb = (
                f"{ref_obj} and chain {ref_chain} "
                f"and polymer.protein and name N+CA+C+O"
            )
            try:
                cmd.align(mobile_bb, target_bb)
            except Exception as exc:
                print(f"[outliers] align failed for {obj_sticks}: {exc}")

    # --- Highlight the outlier residue ---------------------------------
    residue_sel = f"{obj_sticks} and chain {chain} and resi {resnum}"
    try:
        cmd.color(base_color, residue_sel)
        cmd.show("sticks", residue_sel)
    except Exception as exc:
        print(f"[outliers] residue highlight failed for {obj_sticks}: {exc}")

    # --- Sphere-representation duplicate -------------------------------
    cmd.create(obj_spheres, obj_sticks)
    residue_sel_spheres = f"{obj_spheres} and chain {chain} and resi {resnum}"
    cmd.hide("sticks", residue_sel_spheres)
    cmd.show("spheres", residue_sel_spheres)

    # Disable the sphere copy by default so both aren't drawn on top of
    # each other; the user can toggle it on in the object panel.
    cmd.disable(obj_spheres)

    print(
        f"[outliers] {pdb} {resname}{resnum} chain {chain} "
        f"tier={tier} colour={base_color}"
    )

    return obj_sticks, chain


def main():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)

    rows = load_outliers()
    reference = None
    for row in rows:
        try:
            result = visualise_outlier(row, reference=reference)
        except Exception as exc:
            print(f"[outliers] ERROR on row {row}: {exc}")
            continue
        if result is None:
            continue
        if reference is None:
            # The first successfully loaded structure is the alignment
            # reference for every subsequent one.
            reference = result

    cmd.orient()
    cmd.zoom()


main()
