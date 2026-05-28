#!/usr/bin/env python3
"""
fmo_fragment_prep.py
====================
Full pipeline for FMO fragment preparation:
  1. Identify close-contact residues from PIEDA CSV (R < threshold)
  2. Parse GAMESS rank*.inp (INDAT, FMOXYZ, FMOBND, FRGNAM)
  3. Extract per-fragment XYZ with hydrogen caps at bond cuts (1.09 Å)
  4. Write GAMESS MAKEFP input files
  5. Write EFP-EFP pairwise input files for R < threshold pairs

Usage
-----
Single ligand, single rank:
    python fmo_fragment_prep.py --ligand ac_67 --rank 11

All ligands, single rank:
    python fmo_fragment_prep.py --all-ligands --rank 11

Custom threshold:
    python fmo_fragment_prep.py --ligand ac_67 --rank 11 --threshold 1.5

Paths default to the project layout but can be overridden (see --help).
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------------------
# Default project paths (relative to this script's location)
# ---------------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
FMO_OUTPUT   = os.path.join(SCRIPT_DIR, "efmo_lig_frag_template")
PIEDA_DIR    = os.path.join(SCRIPT_DIR, "pieda_csv")
FRAGCHARGE   = os.path.join(SCRIPT_DIR, "fragcharge.txt")
MAKEFP_TPL   = os.path.join(SCRIPT_DIR, "MAKEFP", "makefp_input")
EFP_TPL      = os.path.join(SCRIPT_DIR, "EFP-EFP", "efp_input")
OUT_ROOT     = SCRIPT_DIR

H_BOND_LEN   = 1.09   # Å
NUCLEAR = {
    "H":1.0,"C":6.0,"N":7.0,"O":8.0,"S":16.0,
    "F":9.0,"P":15.0,"Cl":17.0,"Br":35.0
}


# ===========================================================================
# Parsing helpers
# ===========================================================================

def load_frag_charges(path):
    """Return dict {frag_index (1-based): charge} from fragcharge.txt."""
    with open(path) as f:
        content = f.read()
    vals = re.findall(r"-?\d+", content)
    return {i + 1: int(v) for i, v in enumerate(vals)}


def parse_inp(inp_file):
    """
    Parse a GAMESS FMO rank*.inp file.

    Returns
    -------
    frag_atom_ranges : list of [(start, end), ...]
        0-indexed list; frag_atom_ranges[i] = atom ranges for fragment i+1.
    atoms : list of (elem, x, y, z)
        1-indexed (atoms[0] is atom 1).
    bond_cuts : list of (atomA, atomB)
        Both positive, 1-based.
    frgnam : list of str
        Fragment names in order.
    """
    with open(inp_file) as f:
        lines = f.readlines()

    # ── INDAT: state-machine parser ─────────────────────────────────────────
    # Rule: positive = range start.
    #   next negative  → (start, |next|) range, advance 2
    #   next positive  → single atom (start, start), advance 1
    #   next 0         → single atom, let loop handle 0
    # Negative (as start) → treat as single atom |val|.
    indat_start = next(i for i, l in enumerate(lines) if "INDAT(1)" in l)
    indat_nums = []
    i = indat_start
    while i < len(lines):
        l = lines[i]
        if i > indat_start and re.search(r"\$", l):
            break
        indat_nums.extend(int(n) for n in re.findall(r"-?\d+", l))
        i += 1

    frag_atom_ranges = []
    current = []
    j = 1  # skip the leading 0 from INDAT(1)= 0
    while j < len(indat_nums):
        val = indat_nums[j]
        if val == 0:
            if current:
                frag_atom_ranges.append(current)
                current = []
            j += 1
        elif val > 0:
            s = val
            nxt = indat_nums[j + 1] if j + 1 < len(indat_nums) else 0
            if nxt < 0:
                current.append((s, abs(nxt)))
                j += 2
            else:
                current.append((s, s))
                j += 1
        else:
            current.append((abs(val), abs(val)))
            j += 1
    if current:
        frag_atom_ranges.append(current)

    # ── FMOXYZ ──────────────────────────────────────────────────────────────
    fmoxyz_start = next(i for i, l in enumerate(lines) if "$FMOXYZ" in l) + 1
    atoms = []
    i = fmoxyz_start
    while i < len(lines):
        l = lines[i].strip()
        if l.startswith("$"):
            break
        parts = l.split()
        if len(parts) == 5:
            atoms.append((parts[0], float(parts[2]), float(parts[3]), float(parts[4])))
        i += 1

    # ── FMOBND ──────────────────────────────────────────────────────────────
    fmobnd_start = next(i for i, l in enumerate(lines) if "$FMOBND" in l) + 1
    bond_cuts = []
    i = fmobnd_start
    while i < len(lines):
        l = lines[i].strip()
        if l.startswith("$"):
            break
        nums = re.findall(r"-?\d+", l)
        if len(nums) >= 2:
            bond_cuts.append((abs(int(nums[0])), int(nums[1])))
        i += 1

    # ── FRGNAM ──────────────────────────────────────────────────────────────
    frgnam = []
    in_frgnam = False
    for l in lines:
        if "FRGNAM" in l:
            in_frgnam = True
        if in_frgnam:
            frgnam.extend(re.findall(r"[A-Z][a-z]{2,3}-\d+|lig\d+|HOH\d*", l))
            if len(frgnam) >= len(frag_atom_ranges):
                break

    return frag_atom_ranges, atoms, bond_cuts, frgnam


# ===========================================================================
# Fragment extraction
# ===========================================================================

def get_close_residues(pieda_csv, frag_charges, threshold):
    """
    Return set of frag_J indices where R < threshold and
    name_I is lig[1-4] or HOH (excludes lig/HOH from frag_J).
    """
    residues = {}
    with open(pieda_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name_I = row["name_I"].strip()
            R      = float(row["R"])
            if R < threshold and (re.match(r"lig[1-4]$", name_I) or
                                   name_I.startswith("HOH")):
                fj     = int(row["frag_J"])
                name_j = row["name_J"].strip()
                if re.match(r"lig[1-4]$", name_j) or name_j.startswith("HOH"):
                    continue
                if fj not in residues:
                    residues[fj] = (name_j, frag_charges.get(fj, 0))
    return residues  # {frag_idx: (name, Z)}


def cap_H_position(inside_xyz, outside_xyz):
    v = np.array(outside_xyz) - np.array(inside_xyz)
    return np.array(inside_xyz) + H_BOND_LEN * v / np.linalg.norm(v)


def extract_fragment(frag_idx, frag_atom_ranges, atoms, bond_cuts):
    """
    Return (frag_atoms, cap_H_list).
    frag_atoms : list of (elem, x, y, z)
    cap_H_list : list of np.array [x, y, z]
    """
    ranges = frag_atom_ranges[frag_idx - 1]
    frag_set = set()
    for s, e in ranges:
        frag_set.update(range(s, e + 1))

    frag_atoms = [
        (atoms[ai - 1][0], atoms[ai - 1][1], atoms[ai - 1][2], atoms[ai - 1][3])
        for ai in sorted(frag_set)
    ]

    cap_H_list = []
    for a, b in bond_cuts:
        if a in frag_set and b not in frag_set:
            cap_H_list.append(cap_H_position(atoms[a - 1][1:], atoms[b - 1][1:]))
        elif b in frag_set and a not in frag_set:
            cap_H_list.append(cap_H_position(atoms[b - 1][1:], atoms[a - 1][1:]))

    return frag_atoms, cap_H_list


def write_xyz(path, name, frag_idx, Z, frag_atoms, cap_H_list):
    all_atoms = frag_atoms + [("H", *h) for h in cap_H_list]
    with open(path, "w") as out:
        out.write(f"{len(all_atoms)}\n")
        out.write(f"{name}  frag={frag_idx}  Z={Z}  ncap={len(cap_H_list)}\n")
        for elem, x, y, z in all_atoms:
            out.write(f"{elem:<4s}  {x:16.10f}  {y:16.10f}  {z:16.10f}\n")


# ===========================================================================
# MAKEFP input generation
# ===========================================================================

def write_makefp(path, frag_idx, Z, frag_atoms, cap_H_list, template):
    all_atoms = frag_atoms + [("H", *h) for h in cap_H_list]
    inp = template.replace("icharg=Z", f"icharg={Z}")
    inp = inp.replace("fragname=fnumber", f"fragname=f{frag_idx:03d}")
    coord_block = ""
    for elem, x, y, z in all_atoms:
        nuc = NUCLEAR.get(elem, 1.0)
        coord_block += f" {elem:<2s}  {nuc:.1f}   {x:16.10f}   {y:16.10f}   {z:16.10f}\n"
    inp += coord_block + " $end\n"
    with open(path, "w") as out:
        out.write(inp)


# ===========================================================================
# EFP-EFP input generation
# ===========================================================================

def efrag_block(frag_idx, frag_atoms, cap_H_list):
    all_atoms = frag_atoms + [("H", *h) for h in cap_H_list]
    s = f" FRAGNAME=f{frag_idx:03d}\n"
    for k, (elem, x, y, z) in enumerate(all_atoms[:3], start=1):
        s += f" A{k:02d}{elem:<2s}  {x:16.10f}  {y:16.10f}  {z:16.10f}\n"
    return s


def write_efp(path, fi, frag_i_atoms, frag_i_caps, fj, frag_j_atoms, frag_j_caps):
    inp  = " $contrl coord=fragonly $end\n"
    inp += " $system mwords=300 $end\n"
    inp += " $efrag\n"
    inp += " ISCRELEC=1\n"
    inp += efrag_block(fi, frag_i_atoms, frag_i_caps)
    inp += efrag_block(fj, frag_j_atoms, frag_j_caps)
    inp += " $end\n"
    with open(path, "w") as out:
        out.write(inp)


# ===========================================================================
# Main pipeline
# ===========================================================================

def run_pipeline(ligand, rank, threshold, fmo_output, pieda_dir,
                 fragcharge_path, makefp_tpl_path, out_root):

    # ── Locate files ────────────────────────────────────────────────────────
    ligand_dir = os.path.join(fmo_output, ligand)
    inp_files  = glob.glob(os.path.join(ligand_dir, f"*rank{rank}.inp"))
    if not inp_files:
        print(f"  [SKIP] No rank{rank}.inp found in {ligand_dir}")
        return
    inp_file = inp_files[0]

    pieda_csv = os.path.join(pieda_dir, f"{ligand}_pieda.csv")
    if not os.path.exists(pieda_csv):
        print(f"  [SKIP] No PIEDA CSV found: {pieda_csv}")
        return

    print(f"\n{'='*60}")
    print(f"  Ligand : {ligand}   Rank : {rank}   Threshold : {threshold} Å")
    print(f"  INP    : {os.path.basename(inp_file)}")
    print(f"{'='*60}")

    # ── Load charges ────────────────────────────────────────────────────────
    frag_charges = load_frag_charges(fragcharge_path)

    # ── Parse inp ───────────────────────────────────────────────────────────
    frag_atom_ranges, atoms, bond_cuts, frgnam = parse_inp(inp_file)
    nfrag = len(frag_atom_ranges)
    print(f"  Parsed {nfrag} fragments, {len(atoms)} atoms, {len(bond_cuts)} bond cuts")

    # Fragment name lookup (safe)
    def fname(idx):
        return frgnam[idx - 1] if idx <= len(frgnam) else f"frag{idx}"

    # ── Target fragments: residues (R<threshold) + HOH + lig1-4 ────────────
    close_res = get_close_residues(pieda_csv, frag_charges, threshold)

    # Detect HOH and lig frag indices from FRGNAM
    hoh_frags = [i + 1 for i, n in enumerate(frgnam) if n.startswith("HOH")]
    lig_frags  = [i + 1 for i, n in enumerate(frgnam) if re.match(r"lig\d+$", n)]

    target_set = dict(close_res)                          # residues
    for fi in hoh_frags:
        target_set[fi] = (fname(fi), frag_charges.get(fi, 0))
    for fi in lig_frags:
        target_set[fi] = (fname(fi), frag_charges.get(fi, 0))

    print(f"  Target fragments: {len(close_res)} residues + "
          f"{len(hoh_frags)} HOH + {len(lig_frags)} lig = {len(target_set)} total")

    # ── Output directories ───────────────────────────────────────────────────
    tag      = f"rank{rank}"
    xyz_dir  = os.path.join(out_root, f"residue_xyz_{tag}")
    mkfp_dir = os.path.join(out_root, "MAKEFP", f"inp_{tag}")
    efp_dir  = os.path.join(out_root, "EFP-EFP", f"inp_{tag}_R{threshold}")
    for d in (xyz_dir, mkfp_dir, efp_dir):
        os.makedirs(d, exist_ok=True)

    # ── Load MAKEFP template ────────────────────────────────────────────────
    with open(makefp_tpl_path) as f:
        makefp_template = f.read()

    # ── Extract, write XYZ + MAKEFP ────────────────────────────────────────
    frag_data = {}   # frag_idx -> (frag_atoms, cap_H_list)
    for frag_idx, (name, Z) in sorted(target_set.items()):
        fa, caps = extract_fragment(frag_idx, frag_atom_ranges, atoms, bond_cuts)
        frag_data[frag_idx] = (fa, caps)

        xyz_path  = os.path.join(xyz_dir,  f"frag{frag_idx:03d}_{name}_hcap.xyz")
        mkfp_path = os.path.join(mkfp_dir, f"frag{frag_idx:03d}_{name}_makefp.inp")

        write_xyz(xyz_path, name, frag_idx, Z, fa, caps)
        write_makefp(mkfp_path, frag_idx, Z, fa, caps, makefp_template)

    print(f"  XYZ    → {xyz_dir}")
    print(f"  MAKEFP → {mkfp_dir}")

    # ── Write ligand_full.xyz ───────────────────────────────────────────────
    if lig_frags:
        all_lig_idx = set()
        for fi in lig_frags:
            for s, e in frag_atom_ranges[fi - 1]:
                all_lig_idx.update(range(s, e + 1))
        lig_atoms = [
            (atoms[ai-1][0], atoms[ai-1][1], atoms[ai-1][2], atoms[ai-1][3])
            for ai in sorted(all_lig_idx)
        ]
        with open(os.path.join(xyz_dir, "ligand_full.xyz"), "w") as out:
            out.write(f"{len(lig_atoms)}\n")
            out.write(f"ligand (combined)  {ligand} {tag}  Z=0\n")
            for elem, x, y, z in lig_atoms:
                out.write(f"{elem:<4s}  {x:16.10f}  {y:16.10f}  {z:16.10f}\n")

    # ── EFP-EFP: pairs from PIEDA CSV with R < threshold ───────────────────
    our_set = set(frag_data.keys())
    pairs   = set()
    with open(pieda_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            fi = int(row["frag_I"])
            fj = int(row["frag_J"])
            if float(row["R"]) < threshold and fi in our_set and fj in our_set:
                pairs.add((min(fi, fj), max(fi, fj)))

    for fi, fj in sorted(pairs):
        fa_i, caps_i = frag_data[fi]
        fa_j, caps_j = frag_data[fj]
        efp_path = os.path.join(efp_dir, f"efp_f{fi:03d}_f{fj:03d}.inp")
        write_efp(efp_path, fi, fa_i, caps_i, fj, fa_j, caps_j)

    print(f"  EFP-EFP → {efp_dir}  ({len(pairs)} pairs)")


# ===========================================================================
# CLI
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="FMO fragment preparation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ligand",      help="Single ligand name, e.g. ac_67")
    group.add_argument("--all-ligands", action="store_true",
                       help="Process all ligands found in GAMESS_logfiles/")

    parser.add_argument("--rank",      type=int, required=True,
                        help="Rank number to use, e.g. 11")
    parser.add_argument("--threshold", type=float, default=2.0,
                        help="R threshold in Angstrom (default: 2.0)")

    parser.add_argument("--fmo-output",  default=FMO_OUTPUT)
    parser.add_argument("--pieda-dir",   default=PIEDA_DIR)
    parser.add_argument("--fragcharge",  default=FRAGCHARGE)
    parser.add_argument("--makefp-tpl",  default=MAKEFP_TPL)
    parser.add_argument("--out-root",    default=OUT_ROOT)

    args = parser.parse_args()

    if args.all_ligands:
        ligands = sorted(
            os.path.basename(d)
            for d in glob.glob(os.path.join(args.fmo_output, "*"))
            if os.path.isdir(d)
        )
        print(f"Found {len(ligands)} ligands: {ligands}")
    else:
        ligands = [args.ligand]

    for ligand in ligands:
        run_pipeline(
            ligand      = ligand,
            rank        = args.rank,
            threshold   = args.threshold,
            fmo_output  = args.fmo_output,
            pieda_dir   = args.pieda_dir,
            fragcharge_path = args.fragcharge,
            makefp_tpl_path = args.makefp_tpl,
            out_root    = args.out_root,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
