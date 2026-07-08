import os
import csv

# === Step 1: Load probabilities from the *_probs.txt file ===
def load_probs(filepath):
    probs = {}
    in_complex = False
    with open(filepath) as f:
        for line in f:
            if "Complex" in line:
                in_complex = True
                continue
            if in_complex and line.strip().startswith("conf"):
                parts = line.strip().split()
                rank = int(parts[1])
                prob = float(parts[2])
                probs[rank] = prob
    return probs

# === Step 2: Load all 4 energy terms per fragment per rank ===
def load_fragment_energies(filepath):
    frag_data = {}  # {frag: {rank: {'Ees': x, 'E0': y, ...}}}
    current_rank = None
    with open(filepath) as f:
        for line in f:
            if "==========" in line and "Rank" in line:
                current_rank = int(line.strip().split()[2])
                continue
            if line.strip().startswith("Frag") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) == 6:
                frag_id, frag_name = parts[0], parts[1]
                ees, e0, edisp, gsol = map(float, parts[2:])
                key = f"{frag_name}"
                if key not in frag_data:
                    frag_data[key] = {}
                frag_data[key][current_rank] = {
                    'Ees': ees,
                    'E0': e0,
                    'Edisp': edisp,
                    'Gsol': gsol
                }
    return frag_data

# === Step 3: Compute weighted averages using Boltzmann weights ===
def compute_weighted_avg(frag_data, probs):
    avg_data = {}  # final result: {frag: {Ees: x, E0: y, ...}}
    for frag, ranks in frag_data.items():
        totals = {'Ees': 0, 'E0': 0, 'Edisp': 0, 'Gsol': 0}
        for rank, energies in ranks.items():
            weight = probs.get(rank, 0) / 100  # Convert percent to fraction
            for term in energies:
                totals[term] += energies[term] * weight
        avg_data[frag] = totals
    return avg_data

# === Step 4: Save result to CSV ===
def write_csv(output_path, avg_data):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Fragment", "Ees", "E0", "Edisp", "Gsol"])
        for frag, energies in avg_data.items():
            writer.writerow([
                frag,
                f"{energies['Ees']:.3f}",
                f"{energies['E0']:.3f}",
                f"{energies['Edisp']:.3f}",
                f"{energies['Gsol']:.3f}"
            ])

# === Step 5: Loop over all ligand directories ===
root_dir = "."  # or use full path like "/mnt/lustre/koa/lab/tsqc_group/loreab/scheme1"

for ligand in sorted(os.listdir(root_dir)):
    ligand_dir = os.path.join(root_dir, ligand)
    if not os.path.isdir(ligand_dir):
        continue  # skip files

    # Construct input/output paths
    int_file = os.path.join(ligand_dir, f"{ligand}_all_interactions.txt")
    prob_file = os.path.join(ligand_dir, f"{ligand}_probs.txt")
    output_csv = os.path.join(ligand_dir, f"{ligand}_avg_all_energies.csv")

    # Check both input files exist
    if not os.path.exists(int_file) or not os.path.exists(prob_file):
        print(f"Missing files in {ligand} — skipping")
        continue

    # Run full pipeline
    probs = load_probs(prob_file)
    frag_data = load_fragment_energies(int_file)
    avg_data = compute_weighted_avg(frag_data, probs)
    write_csv(output_csv, avg_data)
    print(f"Saved: {output_csv}")
