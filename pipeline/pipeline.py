"""
Tyrosine Phosphorylation Prediction Pipeline
=============================================
Steps:
  1. Scan PDB files to find all Tyrosine (TYR) residues
  2. Extract structural features for each TYR using BioPython
     (equivalent to the original PyMOL feature extraction)
  3. Predict using the saved XGBoost structural model

Usage:
    python pipeline.py                      # processes all PDB files in pdb_files/
    python pipeline.py --pdb 1ABC.pdb       # single PDB file
    python pipeline.py --pdb_dir /path/pdbs

Output:
    output/predictions.csv
"""

import os
import argparse
import warnings
import numpy as np
import pandas as pd
import joblib
import xgboost as xgb          # must be imported before pydssp (torch conflict)
from os.path import isfile, basename

from Bio import PDB
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.vectors import calc_dihedral
import pydssp

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH  = os.path.join(SCRIPT_DIR, 'xgb_structural_model.json')
SCALER_PATH = os.path.join(SCRIPT_DIR, 'xgb_structural_scaler.pkl')
PDB_DIR     = os.path.join(SCRIPT_DIR, 'pdb_files')
OUTPUT_DIR  = os.path.join(SCRIPT_DIR, 'output')
OUTPUT_CSV  = os.path.join(OUTPUT_DIR, 'predictions.csv')

STRUCTURAL_FEATURES = [
    "Area", "Number of Residues", "Number of Atoms", "Number of COOH Atoms",
    "ss", "phi", "psi", "Length of PDB", "Mut Res/Length"
]

ACIDIC = {"GLU", "ASP"}


# ── Step 1: Find TYR residues ─────────────────────────────────────────────────
def find_tyrosines(structure):
    """Return list of Residue objects that are TYR."""
    return [
        res for model in structure
        for chain in model
        for res in chain
        if res.get_resname() == "TYR" and PDB.is_aa(res)
    ]


# ── SASA via Shrake–Rupley (built into BioPython) ─────────────────────────────
def compute_sasa(structure):
    """Returns dict mapping residue full_id → SASA (Å²)."""
    sr = PDB.SASA.ShrakeRupley()
    sr.compute(structure, level="R")
    sasa_map = {}
    for model in structure:
        for chain in model:
            for res in chain:
                sasa_map[res.full_id] = res.sasa
    return sasa_map


# ── DSSP secondary structure via pydssp (no binary needed) ───────────────────
def compute_dssp(pdb_path):
    """
    Returns dict mapping (chain_id, resnum) → ss int (0=loop, 1=helix, 2=sheet).
    Uses pydssp, a pure-Python DSSP implementation.
    """
    try:
        coords   = pydssp.read_pdbtext(open(pdb_path).read())
        ss_array = pydssp.assign(coords, out_type='onehot')
        # ss_array shape: (N_residues, 3) one-hot [loop, helix, sheet]
        # argmax: 0=loop → 0, 1=helix → 1, 2=sheet → 2  (maps directly)
        ss_indices = ss_array.argmax(axis=-1)

        # Map back to (chain, resnum) using BioPython residue order
        parser    = PDBParser(QUIET=True)
        structure = parser.get_structure("tmp", pdb_path)
        residues  = [
            res for model in structure
            for chain in model
            for res in chain
            if PDB.is_aa(res)
        ]
        dssp_map = {}
        for i, res in enumerate(residues):
            if i >= len(ss_indices):
                break
            chain_id = res.get_parent().get_id()
            resnum   = res.get_id()[1]
            dssp_map[(chain_id, resnum)] = int(ss_indices[i])  # 0=loop,1=helix,2=sheet
        return dssp_map
    except Exception as e:
        print(f"  [WARN] pydssp failed ({e}), ss will be 0 for all residues")
        return {}


# ── Step 2: Extract features for one TYR ─────────────────────────────────────
def extract_features(pdb_path, structure, tyr_res, sasa_map, dssp_map):
    pdb_name = basename(pdb_path)
    model    = structure[0]

    chain_id = tyr_res.get_parent().get_id()
    resnum   = tyr_res.get_id()[1]

    # Total residues (alpha carbons = one per residue)
    total_residues = sum(
        1 for chain in model
        for res in chain
        if PDB.is_aa(res) and "CA" in res
    )

    # SASA for this residue
    area = sasa_map.get(tyr_res.full_id, 0.0)

    # Neighbourhood within 5 Å
    all_atoms = list(model.get_atoms())
    ns = NeighborSearch(all_atoms)
    tyr_atoms = list(tyr_res.get_atoms())
    neighbour_atoms = set()
    for atom in tyr_atoms:
        for nb in ns.search(atom.get_coord(), 5.0, level='A'):
            if nb.get_parent() != tyr_res:
                neighbour_atoms.add(nb)

    natom = len(neighbour_atoms)
    neighbour_residues = {a.get_parent() for a in neighbour_atoms if PDB.is_aa(a.get_parent())}
    nres  = len(neighbour_residues)
    ncooh = sum(1 for r in neighbour_residues if r.get_resname() in ACIDIC)

    # Secondary structure
    ss = dssp_map.get((chain_id, resnum), 0)

    # Dihedral angles phi / psi
    phi, psi = None, None
    try:
        chain   = model[chain_id]
        res_ids = [r.get_id() for r in chain if PDB.is_aa(r)]
        idx     = res_ids.index(tyr_res.get_id())

        if idx > 0:
            prev_res = chain[res_ids[idx - 1]]
            if "C" in prev_res and "N" in tyr_res and "CA" in tyr_res and "C" in tyr_res:
                phi = np.degrees(calc_dihedral(
                    prev_res["C"].get_vector(),
                    tyr_res["N"].get_vector(),
                    tyr_res["CA"].get_vector(),
                    tyr_res["C"].get_vector()
                ))

        if idx < len(res_ids) - 1:
            next_res = chain[res_ids[idx + 1]]
            if "N" in tyr_res and "CA" in tyr_res and "C" in tyr_res and "N" in next_res:
                psi = np.degrees(calc_dihedral(
                    tyr_res["N"].get_vector(),
                    tyr_res["CA"].get_vector(),
                    tyr_res["C"].get_vector(),
                    next_res["N"].get_vector()
                ))
    except Exception as e:
        print(f"  [WARN] Dihedral error for chain={chain_id} resi={resnum}: {e}")

    mut_res_length_ratio = resnum / total_residues if total_residues else 0

    return {
        "pdb":                  pdb_name,
        "chain":                chain_id,
        "Residue Number":       resnum,
        "Area":                 area,
        "Number of Residues":   nres,
        "Number of Atoms":      natom,
        "Number of COOH Atoms": ncooh,
        "ss":                   ss,
        "phi":                  phi,
        "psi":                  psi,
        "Length of PDB":        total_residues,
        "Mut Res/Length":       mut_res_length_ratio,
    }


# ── Step 3: Predict ───────────────────────────────────────────────────────────
def predict(features_df, model, scaler):
    X = features_df[STRUCTURAL_FEATURES].copy()
    X_scaled = scaler.transform(X)
    return model.predict(X_scaled)


# ── Main ──────────────────────────────────────────────────────────────────────
def main(pdb_paths):
    print(f"Loading model  : {MODEL_PATH}")
    print(f"Loading scaler : {SCALER_PATH}")
    model = xgb.XGBRegressor()
    model.load_model(MODEL_PATH)
    scaler = joblib.load(SCALER_PATH)

    parser = PDBParser(QUIET=True)
    all_records = []

    for pdb_path in pdb_paths:
        pdb_name = basename(pdb_path)
        print(f"\n{'─'*50}")
        print(f"Processing: {pdb_name}")

        structure = parser.get_structure(pdb_name, pdb_path)

        # Step 1: find TYR residues
        tyrosines = find_tyrosines(structure)
        if not tyrosines:
            print(f"  No TYR residues found.")
            continue
        print(f"  Found {len(tyrosines)} TYR residue(s): "
              + ", ".join(f"chain={r.get_parent().get_id()} resi={r.get_id()[1]}" for r in tyrosines))

        # Pre-compute SASA and DSSP once per structure
        sasa_map = compute_sasa(structure)
        dssp_map = compute_dssp(pdb_path)

        # Step 2: extract features
        for tyr_res in tyrosines:
            chain_id = tyr_res.get_parent().get_id()
            resnum   = tyr_res.get_id()[1]
            print(f"  Extracting features: chain={chain_id} resi={resnum}")
            record = extract_features(pdb_path, structure, tyr_res, sasa_map, dssp_map)
            all_records.append(record)

    if not all_records:
        print("\nNo features extracted. Check PDB files in pdb_files/")
        return

    features_df = pd.DataFrame(all_records)

    # Drop rows missing phi/psi (terminal residues)
    before = len(features_df)
    features_df = features_df.dropna(subset=STRUCTURAL_FEATURES)
    dropped = before - len(features_df)
    if dropped:
        print(f"\n[INFO] Dropped {dropped} residue(s) with missing phi/psi (terminal residues).")

    if features_df.empty:
        print("No complete feature rows to predict on.")
        return

    # Step 3: predict
    features_df['Predicted_Exp'] = predict(features_df, model, scaler)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out_cols = ['pdb', 'chain', 'Residue Number'] + STRUCTURAL_FEATURES + ['Predicted_Exp']
    features_df[out_cols].to_csv(OUTPUT_CSV, index=False)

    print(f"\n{'='*50}")
    print(f"Done. Results saved to: {OUTPUT_CSV}")
    print(f"{'='*50}")
    print(features_df[['pdb', 'chain', 'Residue Number', 'Predicted_Exp']].to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tyrosine phosphorylation prediction pipeline')
    group  = parser.add_mutually_exclusive_group()
    group.add_argument('--pdb',     type=str, help='Single PDB file path')
    group.add_argument('--pdb_dir', type=str, help='Directory of PDB files (default: pdb_files/)')
    args = parser.parse_args()

    if args.pdb:
        pdb_paths = [args.pdb]
    else:
        pdb_dir = args.pdb_dir or PDB_DIR
        pdb_paths = [
            os.path.join(pdb_dir, f)
            for f in sorted(os.listdir(pdb_dir))
            if f.lower().endswith('.pdb') and isfile(os.path.join(pdb_dir, f))
        ]
        if not pdb_paths:
            print(f"No .pdb files found in {pdb_dir}")
            exit(1)

    main(pdb_paths)
