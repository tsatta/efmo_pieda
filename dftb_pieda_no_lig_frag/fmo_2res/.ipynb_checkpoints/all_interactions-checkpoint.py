#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import os
import sys

# Arguments
if len(sys.argv) != 3:
    print("Usage: python binding_site_info.py <log_file_path> <active_atoms_file>")
    sys.exit(1)

log_file_path = sys.argv[1]
active_atoms_path = sys.argv[2]
rank_match = re.search(r'rank(\d+)', log_file_path)
rank_label = f"Rank {rank_match.group(1)}" if rank_match else "Unknown Rank"

# Output file 
subdir_name = os.path.basename(os.path.dirname(log_file_path))
output_file_path = os.path.join(os.path.dirname(log_file_path), f"{subdir_name}_all_interactions.txt")

# Load active atoms
active_atoms = []
with open(active_atoms_path, "r") as f:
    for line in f:
        if line.strip() and not line.startswith("Index"):
            parts = line.split()
            atom = {
                "index": int(parts[0]),
                "element": parts[1],
                "x": float(parts[2]),
                "y": float(parts[3]),
                "z": float(parts[4]),
                "frag": int(parts[5]),
                "fragname": parts[6]
            }
            active_atoms.append(atom)

# Find ligand fragment ID
ligand_frag = None
with open(log_file_path, "r") as f:
    for line in f:
        if re.search(r"\bLIG\b|\bLIG\d+\b|\bLIG-\d+\b", line):
            parts = line.strip().split()
            if parts[0].isdigit():
                ligand_frag = int(parts[0])
                break

if ligand_frag is None:
    raise ValueError("Could not identify ligand fragment in the log file.")

# Parse FMO properties section
interactions = {}
fmo_section_started = False

with open(log_file_path, "r") as f:
    for line in f:
        if "Two-body FMO properties" in line:
            fmo_section_started = True
            continue
        if fmo_section_started:
            if re.match(r"\s*---+", line):
                continue
            match = re.match(r"\s*(\d+)\s+(\d+)\s+\S+\s+[-+]?\d+\s+[-+]?\d*\.\d*", line)
            if match:
                parts = line.split()
                try:
                    i_frag = int(parts[0])
                    j_frag = int(parts[1])
                    ees = float(parts[9])
                    e0 = float(parts[10])
                    edisp = float(parts[12])
                    gsol = float(parts[13])
                    interactions[(i_frag, j_frag)] = {
                        "Ees": ees, "E0": e0, "Edisp": edisp, "Gsol": gsol
                    }
                    interactions[(j_frag, i_frag)] = interactions[(i_frag, j_frag)]
                except (IndexError, ValueError):
                    continue

# Match protein fragment to ligand fragment and write output 
with open(output_file_path, "a") as out:
    out.write(f"\n========== {rank_label} ==========\n\n")
    out.write("Frag   FragName        Ees       E0    Edisp     Gsol\n")

    unique_lines = set()
    for atom in active_atoms:
        frag_id = atom["frag"]
        fragname = atom["fragname"]
        if (frag_id, ligand_frag) in interactions:
            data = interactions[(frag_id, ligand_frag)]
            line = f"{frag_id:<6} {fragname:<12} {data['Ees']:8.3f} {data['E0']:8.3f} {data['Edisp']:8.3f} {data['Gsol']:8.3f}"
            if line not in unique_lines:
                out.write(line + "\n")
                unique_lines.add(line)

print(f"Output written to: {output_file_path}")


# In[ ]:




