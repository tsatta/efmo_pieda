#!/usr/bin/env python3
"""
make_subsystem_fmo.py
=====================
Build a reduced GAMESS FMO input file containing only:
  - close-contact residues (from close_residues_R2.csv, R < threshold)
  - HOH fragment
  - lig1-4 fragments

Bond cuts that cross the subsystem boundary become H-caps (1.09 Å).
Internal bond cuts (both fragments in subsystem) are preserved in $FMOBND.

Usage
-----
    python make_subsystem_fmo.py --ligand ac_67 --rank 11
    python make_subsystem_fmo.py --all-ligands --rank 11
    python make_subsystem_fmo.py --all-ligands --all-ranks
"""

import argparse, csv, glob, os, re
import numpy as np
from collections import defaultdict

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FMO_OUTPUT  = os.path.join(SCRIPT_DIR, "efmo_lig_frag_template")
PIEDA_DIR   = os.path.join(SCRIPT_DIR, "pieda_csv")
FRAGCHARGE  = os.path.join(SCRIPT_DIR, "fragcharge.txt")
OUT_ROOT    = os.path.join(SCRIPT_DIR, "subsystem_fmo")
H_BOND_LEN  = 1.09

NUCLEAR = {"H":1.0,"C":6.0,"N":7.0,"O":8.0,"S":16.0,
           "F":9.0,"P":15.0,"Cl":17.0,"Br":35.0}

# ── Parsers (same as fmo_fragment_prep.py) ──────────────────────────────────

def load_frag_charges(path):
    with open(path) as f:
        vals = re.findall(r"-?\d+", f.read())
    return {i+1: int(v) for i,v in enumerate(vals)}

def parse_inp(inp_file):
    with open(inp_file) as f:
        lines = f.readlines()

    # INDAT state-machine
    indat_start = next(i for i,l in enumerate(lines) if "INDAT(1)" in l)
    indat_nums = []
    i = indat_start
    while i < len(lines):
        l = lines[i]
        if i > indat_start and re.search(r"\$", l): break
        indat_nums.extend(int(n) for n in re.findall(r"-?\d+", l))
        i += 1
    frag_atom_ranges = []
    current = []
    j = 1
    while j < len(indat_nums):
        val = indat_nums[j]
        if val == 0:
            if current: frag_atom_ranges.append(current); current = []
            j += 1
        elif val > 0:
            s = val
            nxt = indat_nums[j+1] if j+1 < len(indat_nums) else 0
            if nxt < 0: current.append((s, abs(nxt))); j += 2
            else: current.append((s, s)); j += 1
        else:
            current.append((abs(val), abs(val))); j += 1
    if current: frag_atom_ranges.append(current)

    # FMOXYZ
    fmoxyz_start = next(i for i,l in enumerate(lines) if "$FMOXYZ" in l) + 1
    atoms = []
    i = fmoxyz_start
    while i < len(lines):
        l = lines[i].strip()
        if l.startswith("$"): break
        parts = l.split()
        if len(parts) == 5:
            atoms.append((parts[0], float(parts[2]), float(parts[3]), float(parts[4])))
        i += 1

    # FMOBND
    fmobnd_start = next(i for i,l in enumerate(lines) if "$FMOBND" in l) + 1
    bond_cuts = []
    i = fmobnd_start
    while i < len(lines):
        l = lines[i].strip()
        if l.startswith("$"): break
        nums = re.findall(r"-?\d+", l)
        if len(nums) >= 2:
            bond_cuts.append((abs(int(nums[0])), int(nums[1])))
        i += 1

    # FRGNAM
    frgnam = []
    in_frgnam = False
    for l in lines:
        if "FRGNAM" in l: in_frgnam = True
        if in_frgnam:
            frgnam.extend(re.findall(r"[A-Z][a-z]{2,3}-\d+|lig\d+|HOH\d*", l))
            if len(frgnam) >= len(frag_atom_ranges): break

    return lines, frag_atom_ranges, atoms, bond_cuts, frgnam

def get_close_residues(pieda_csv, frag_charges, threshold):
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
                residues[fj] = (name_j, frag_charges.get(fj, 0))
    return residues

def cap_H_pos(inside_xyz, outside_xyz):
    v = np.array(outside_xyz) - np.array(inside_xyz)
    return np.array(inside_xyz) + H_BOND_LEN * v / np.linalg.norm(v)

# ── Subsystem builder ────────────────────────────────────────────────────────

def build_subsystem(inp_file, pieda_csv, frag_charges, threshold, out_path):
    lines, frag_atom_ranges, atoms, bond_cuts, frgnam = parse_inp(inp_file)

    # Target fragment set
    close_res = get_close_residues(pieda_csv, frag_charges, threshold)
    hoh_frags = [i+1 for i,n in enumerate(frgnam) if n.startswith("HOH")]
    lig_frags  = [i+1 for i,n in enumerate(frgnam) if re.match(r"lig\d+$", n)]
    target_set = set(close_res) | set(hoh_frags) | set(lig_frags)

    def fname(idx):
        return frgnam[idx-1] if idx <= len(frgnam) else f"frag{idx}"

    # For each target fragment: collect original atoms + determine H-caps
    # atom_to_frag: old atom index -> fragment index
    atom_to_frag = {}
    for fi, ranges in enumerate(frag_atom_ranges):
        for s, e in ranges:
            for ai in range(s, e+1):
                atom_to_frag[ai] = fi + 1

    frag_info = {}  # frag_idx -> {'atoms': [...], 'caps': [...]}
    for fi in sorted(target_set):
        ranges = frag_atom_ranges[fi-1]
        frag_set = set()
        for s, e in ranges:
            frag_set.update(range(s, e+1))
        orig_atoms = [(atoms[ai-1][0], atoms[ai-1][1], atoms[ai-1][2], atoms[ai-1][3])
                      for ai in sorted(frag_set)]
        caps = []
        for a, b in bond_cuts:
            if a in frag_set and b not in frag_set:
                # cap only if the other fragment is NOT in the target list
                if atom_to_frag.get(b) not in target_set:
                    caps.append(cap_H_pos(atoms[a-1][1:], atoms[b-1][1:]))
            elif b in frag_set and a not in frag_set:
                if atom_to_frag.get(a) not in target_set:
                    caps.append(cap_H_pos(atoms[b-1][1:], atoms[a-1][1:]))
        frag_info[fi] = {'atoms': orig_atoms, 'caps': caps}

    # Build new sequential atom list and index mapping
    # new_fmoxyz: list of (elem, x, y, z)
    # new_indat: list of (new_start, new_end) per kept fragment
    new_fmoxyz = []
    new_indat  = []   # [(new_s, new_e), ...] per kept frag (in kept order)
    old_to_new = {}   # old 1-based atom idx -> new 1-based atom idx
    new_idx = 1

    # First pass: map original atoms
    for fi in sorted(target_set):
        ranges = frag_atom_ranges[fi-1]
        for s, e in ranges:
            for ai in range(s, e+1):
                old_to_new[ai] = new_idx
                new_fmoxyz.append(atoms[ai-1])
                new_idx += 1

    # Second pass: build INDAT entries (original ranges + cap H appended)
    for fi in sorted(target_set):
        ranges = frag_atom_ranges[fi-1]
        caps   = frag_info[fi]['caps']
        # Convert original ranges to new indices
        new_ranges = []
        for s, e in ranges:
            new_ranges.append((old_to_new[s], old_to_new[e]))
        # Append cap H atoms
        if caps:
            cap_start = new_idx
            for h in caps:
                new_fmoxyz.append(("H", h[0], h[1], h[2]))
                new_idx += 1
            cap_end = new_idx - 1
            new_ranges.append((cap_start, cap_end))
        new_indat.append(new_ranges)

    # Internal bond cuts (both atoms in target_set) → update indices
    target_atoms = set()
    for fi in target_set:
        for s, e in frag_atom_ranges[fi-1]:
            target_atoms.update(range(s, e+1))

    new_bond_cuts = []
    for a, b in bond_cuts:
        if a in target_atoms and b in target_atoms:
            new_bond_cuts.append((old_to_new[a], old_to_new[b]))

    # ── Assemble output file ─────────────────────────────────────────────────
    kept_frags = sorted(target_set)
    nfrag_new  = len(kept_frags)
    nbnd_new   = len(new_bond_cuts)

    out = []

    # Header sections: $CONTRL through $SCF (lines 0–37 for ac_67)
    # Find end of $SCF (line before $FMO)
    # Match the standalone $FMO block (not the inline '$fmo resdim=... $end')
    fmo_start = next(i for i,l in enumerate(lines)
                     if re.match(r"\s*\$FMO\s*$", l, re.I))
    fmobnd_start = next(i for i,l in enumerate(lines) if "$FMOBND" in l)
    fmoxyz_start = next(i for i,l in enumerate(lines) if "$FMOXYZ" in l)
    fmoxyz_end   = next(i for i,l in enumerate(lines) if "$END" in l.upper() and i > fmoxyz_start)

    # 1) All lines before $FMO (header)
    out.extend(lines[:fmo_start])

    # 2) $FMO block
    out.append(" $FMO\n")
    out.append("      SCFTYP(1)=RHF\n")
    out.append("      MODMUL=0\n")
    out.append("      MAXCAO=5\n")
    out.append("      MODGRD=26\n")
    out.append("      MODMUL=0\n")
    out.append(f"      MAXBND={nbnd_new}\n")
    out.append(f"      NFRAG={nfrag_new}\n")

    MAX_COL = 72

    def wrap_line(prefix, tokens, sep="", cont_indent=None):
        """Join tokens onto lines; start a new line when col > MAX_COL."""
        if cont_indent is None:
            cont_indent = " " * len(prefix)
        lines_out = []
        cur = prefix
        for k, tok in enumerate(tokens):
            piece = ("" if k == 0 else sep) + tok
            if k > 0 and len(cur) + len(piece) > MAX_COL:
                lines_out.append(cur + "\n")
                cur = cont_indent + tok
            else:
                cur += piece
        lines_out.append(cur + "\n")
        return lines_out

    # ICHARG
    charges = [frag_charges.get(fi, 0) for fi in kept_frags]
    tokens  = [str(z) for z in charges]
    out.extend(wrap_line("      ICHARG(1)=   ", tokens, sep=",  ",
                         cont_indent="                   "))

    # FRGNAM
    names = [fname(fi) for fi in kept_frags]
    out.extend(wrap_line("      FRGNAM(1)= ", names, sep=", ",
                         cont_indent="                 "))

    # INDAT
    out.append("      INDAT(1)= 0\n")
    indent = "              "
    for ranges in new_indat:
        tokens = [f"{s:6d}  {-e:6d}" for s, e in ranges] + ["     0"]
        lines_out = wrap_line(indent, tokens, sep="  ",
                              cont_indent=indent)
        out.extend(lines_out)

    out.append(" $END\n")

    # 3) $FMOBND block
    out.append(" $FMOBND\n")
    for a, b in new_bond_cuts:
        out.append(f"      {-a:8d}  {b:8d}    6-31+Gd\n")
    out.append(" $END\n")

    # 4) $DATA block (keep as-is)
    data_start = next(i for i,l in enumerate(lines) if "$DATA" in l.upper() and "$END" not in l.upper())
    data_end   = next(i for i,l in enumerate(lines) if "$END" in l.upper() and i > data_start)
    out.extend(lines[data_start : data_end+1])

    # 5) $FMOXYZ
    out.append(" $FMOXYZ\n")
    for elem, x, y, z in new_fmoxyz:
        nuc = NUCLEAR.get(elem, 1.0)
        out.append(f" {elem:<2s}  {nuc:.1f}  {x:18.10f}  {y:18.10f}  {z:18.10f}\n")
    out.append(" $END\n")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.writelines(out)

    print(f"  Kept {nfrag_new} frags, {len(new_fmoxyz)} atoms, {nbnd_new} internal cuts")
    print(f"  Written → {out_path}")


# ── Rank discovery ───────────────────────────────────────────────────────────

def get_available_ranks(ligand_dir):
    ranks = []
    for f in glob.glob(os.path.join(ligand_dir, "*rank*.inp")):
        m = re.search(r"rank(\d+)\.inp$", os.path.basename(f))
        if m: ranks.append(int(m.group(1)))
    return sorted(set(ranks))


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Build subsystem FMO input files")
    lig_group = parser.add_mutually_exclusive_group(required=True)
    lig_group.add_argument("--ligand",      help="Single ligand, e.g. ac_67")
    lig_group.add_argument("--all-ligands", action="store_true")

    rank_group = parser.add_mutually_exclusive_group(required=True)
    rank_group.add_argument("--rank",      type=int)
    rank_group.add_argument("--all-ranks", action="store_true")

    parser.add_argument("--threshold", type=float, default=2.0)
    parser.add_argument("--fmo-output", default=FMO_OUTPUT)
    parser.add_argument("--pieda-dir",  default=PIEDA_DIR)
    parser.add_argument("--fragcharge", default=FRAGCHARGE)
    parser.add_argument("--out-root",   default=OUT_ROOT)
    args = parser.parse_args()

    ligands = (
        sorted(os.path.basename(d) for d in glob.glob(os.path.join(args.fmo_output, "*"))
               if os.path.isdir(d))
        if args.all_ligands else [args.ligand]
    )

    frag_charges = load_frag_charges(args.fragcharge)
    total = 0

    for ligand in ligands:
        ligand_dir = os.path.join(args.fmo_output, ligand)
        ranks = get_available_ranks(ligand_dir) if args.all_ranks else [args.rank]

        pieda_csv = os.path.join(args.pieda_dir, f"{ligand}_pieda.csv")
        if not os.path.exists(pieda_csv):
            print(f"[SKIP] No PIEDA CSV for {ligand}"); continue

        for rank in ranks:
            inps = glob.glob(os.path.join(ligand_dir, f"*rank{rank}.inp"))
            if not inps:
                print(f"[SKIP] {ligand} rank{rank}: no .inp found"); continue
            inp_file = inps[0]

            out_path = os.path.join(
                args.out_root, ligand,
                f"{ligand}_rank{rank}_subsystem.inp"
            )
            print(f"\n{ligand}  rank{rank}")
            build_subsystem(inp_file, pieda_csv, frag_charges,
                            args.threshold, out_path)
            total += 1

    print(f"\nDone. {total} file(s) written.")


if __name__ == "__main__":
    main()
