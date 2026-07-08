# PIEDA using EFP-EFP
### May 28, 2026
Starting structures are PLQMVM2_FMO_OPT
1. fmo_fragment_prep.py Usage
2. Folder description
   
----
1. fmo_fragment_prep.py Usage

- Single ligand, single rank:
    ```python fmo_fragment_prep.py --ligand ac_67 --rank 11```

- All ligands, single rank:
    ```python fmo_fragment_prep.py --all-ligands --rank 11```

- Custom threshold:
    ```python fmo_fragment_prep.py --ligand ac_67 --rank 11 --threshold 1.5```

- Paths default to the project layout but can be overridden (see --help).

----
2. Folder description
   - PLQMVM2_FMO_OPT_PPE: contains the probability of all complexes
   - dftb_pieda_lig_frag: FMO/DFTB PIEDA log files. Ligands are fragmented.
   - dftb_pieda_no_lig_frag: FMO/DFTB PIEDA inp files. Using full ligands (no fragmentation)
   - lig_fmo_template: free ligands with fragmentation
   - efmo_lig_frag_template: efmo inputs for pieda
----
