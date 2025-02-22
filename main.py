import os
import re
import requests
import subprocess
import pickle
import pandas as pd
from pymol_feature import extract_pymol_features
from EvoEF_stability import compute_evoef_stability
from Foldx_stability import compute_foldx_stability

# Define paths
PDB_DIR = "./pdb_files"
MODEL_PATH = "xgboost_model.pkl"

############TO DO##############
## 1. modify input list to include all the protein uniprot id and mutations
## 2. define the output csv file name
## 3. specify the path for EvoEF and FoldX


# List of UniProt IDs and mutations
input_data = [
    ("P04637", "YA103E"),
]

# Output csv path
OUTPUT_CSV = "predictions.csv"

# EVOEF path
evoef_path="../EvoEF_m/EvoEF"

# foldx path
foldx_path="./foldx_20251231"

######### END TODO ######### 

# AlphaFold structure retrieval
ALPHAFOLD_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb"

# Load pre-trained XGBoost model
with open(MODEL_PATH, "rb") as model_file:
    xgb_model = pickle.load(model_file)

def download_alphafold_structure(uniprot_id, pdb_dir):
    """Download AlphaFold structure using UniProt ID and save it in PDB format."""
    url = ALPHAFOLD_URL_TEMPLATE.format(uniprot_id)
    pdb_path = os.path.join(pdb_dir, f"{uniprot_id}.pdb")

    if os.path.exists(pdb_path):
        print(f"Structure for {uniprot_id} already exists. Skipping download.")
        return pdb_path

    print(f"Downloading AlphaFold structure for {uniprot_id}...")
    response = requests.get(url)

    if response.status_code == 200:
        os.makedirs(pdb_dir, exist_ok=True)
        with open(pdb_path, "wb") as f:
            f.write(response.content)
        print(f"Saved structure for {uniprot_id} at {pdb_path}")
        return pdb_path
    else:
        print(f"Failed to download structure for {uniprot_id}. HTTP Status: {response.status_code}")
        return None



# Store results
results = []

for uniprot_id, mutation in input_data:
    pdb_file = download_alphafold_structure(uniprot_id, PDB_DIR)
    match = re.search(r'^[A-Za-z]{2}(\d+)[A-Za-z]$', mutation)
    residue_num = 0
    if match:
        residue_num = int(match.group(1))
    

    if pdb_file is None:
        print(f"Skipping {uniprot_id}: Structure not available.")
        continue

    print(f"Processing {uniprot_id}, residue {residue_num}...")

    # Step 1: Extract features using PyMOL
    pymol_features = extract_pymol_features(pdb_file, residue_num)

    # Step 2: Compute EvoEF stability
    evoef_stability = compute_evoef_stability(pdb_file, mutation)

    # Step 3: Compute FoldX stability
    foldx_stability = compute_foldx_stability(pdb_file, mutation)

    # Step 4: Prepare feature vector
    feature_vector = {
        "Area": pymol_features["area"],
        "Number of Contacts": pymol_features["num_residues"],
        "Phi": pymol_features["phi"],
        "Psi": pymol_features["psi"],
        "EvoEF": evoef_stability,
        "FoldX": foldx_stability,
    }

    # Convert to DataFrame format
    feature_df = pd.DataFrame([feature_vector])

    # Step 5: Predict using XGBoost
    prediction = xgb_model.predict(feature_df)[0]

    results.append({
        "UniProt ID": uniprot_id,
        "Residue": residue_num,
        "Prediction": prediction
    })
print(results)
# Save results to CSV
result_df = pd.DataFrame(results)
result_df.to_csv(OUTPUT_CSV, index=False)

print(f"Predictions saved to {OUTPUT_CSV}")
