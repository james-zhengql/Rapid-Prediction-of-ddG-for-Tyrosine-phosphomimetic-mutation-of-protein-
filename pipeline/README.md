# Tyrosine Phosphorylation Effect Prediction Pipeline

Predicts the experimental effect of tyrosine (TYR) phosphorylation from protein structure alone. Given one or more PDB files, the pipeline automatically finds every tyrosine residue, computes structural features, and outputs a predicted effect score for each site.

---

## How it works

```
PDB file(s)
    │
    ▼
Step 1 — Find TYR residues
    Parse ATOM/HETATM records and collect every unique
    (chain, residue number) pair where residue name = TYR.
    │
    ▼
Step 2 — Extract structural features (BioPython)
    For each TYR residue:
    ├─ Area             Solvent-accessible surface area (Shrake–Rupley, Å²)
    ├─ Number of Atoms  All non-self atoms within 5 Å
    ├─ Number of Res    Unique residues within 5 Å
    ├─ COOH Atoms       Acidic neighbours (GLU / ASP) within 5 Å
    ├─ ss               Secondary structure: 0=loop, 1=helix, 2=sheet (pydssp)
    ├─ phi / psi        Backbone dihedral angles (degrees)
    ├─ Length of PDB    Total residues in the structure
    └─ Mut Res/Length   Residue number / total residues (positional ratio)
    │
    ▼
Step 3 — Predict
    Features are StandardScaler-normalised then passed to a trained
    XGBoost regressor. Outputs a continuous Predicted_Exp score.
    │
    ▼
output/predictions.csv
```

---

## Model background

The XGBoost model was trained on the Tsutsumi dataset (`output_features_tsu(all).csv`) using only structural features — no EvoEF or FoldX energy terms. 5-fold cross-validation performance on tyrosine (Mut = Y) entries:

| Metric   | Score  |
|----------|--------|
| Pearson  | 0.626  |
| Spearman | 0.611  |
| RMSE     | 1.110  |
| MSE      | 1.232  |

Feature importances (most → least):

| Feature              | Importance |
|----------------------|------------|
| Area                 | 0.434      |
| Number of Residues   | 0.109      |
| Mut Res/Length       | 0.103      |
| psi                  | 0.100      |
| Length of PDB        | 0.072      |
| Number of COOH Atoms | 0.060      |
| phi                  | 0.057      |
| Number of Atoms      | 0.046      |
| ss                   | 0.019      |

---

## Folder structure

```
pipeline/
├── pipeline.py                  # main script
├── xgb_structural_model.json    # trained XGBoost model (native JSON format)
├── xgb_structural_scaler.pkl    # fitted StandardScaler
├── pdb_files/                   # put your input PDB files here
│   └── 1A22.pdb                 # example
└── output/
    └── predictions.csv          # results written here after each run
```

---

## Setup

Create and activate a Python virtual environment, then install dependencies:

```bash
python3 -m venv venv
source venv/bin/activate

pip install biopython pandas numpy scikit-learn xgboost joblib pydssp
```

> **Note:** PyMOL is not required. All structural calculations use BioPython and pydssp.
>
> **Import order matters:** `xgboost` must be imported before `pydssp` in any script that uses both — pydssp pulls in PyTorch which can cause a segfault if XGBoost is loaded after it.

---

## Usage

Run from inside the `pipeline/` directory (or use `--pdb` / `--pdb_dir` for other paths):

```bash
# All .pdb files in pdb_files/
python pipeline.py

# Single file anywhere on disk
python pipeline.py --pdb /path/to/1XYZ.pdb

# Custom folder
python pipeline.py --pdb_dir /path/to/my_pdbs/
```

---

## Output

Results are written to `output/predictions.csv`. Each row is one TYR residue.

| Column               | Description                                        |
|----------------------|----------------------------------------------------|
| pdb                  | Source PDB filename                                |
| chain                | Chain ID                                           |
| Residue Number       | Residue sequence number                            |
| Area                 | Solvent-accessible surface area (Å²)               |
| Number of Residues   | Residues within 5 Å of the TYR                    |
| Number of Atoms      | Atoms within 5 Å of the TYR                       |
| Number of COOH Atoms | Acidic residues (GLU/ASP) within 5 Å              |
| ss                   | Secondary structure (0=loop, 1=helix, 2=sheet)     |
| phi                  | Backbone phi dihedral angle (degrees)              |
| psi                  | Backbone psi dihedral angle (degrees)              |
| Length of PDB        | Total residue count in the structure               |
| Mut Res/Length       | Positional ratio (residue number / total residues) |
| Predicted_Exp        | Predicted phosphorylation effect score             |

Terminal residues (N- or C-terminus) are dropped automatically because phi or psi cannot be computed without a flanking residue.

---

## Re-training the model

To retrain on updated data, run `xgb_structural_only.py` from the `Final/` directory, then copy the outputs into `pipeline/`:

```bash
cd ../                         # Final/
python xgb_structural_only.py  # saves xgb_structural_model.json and xgb_structural_scaler.pkl
cp xgb_structural_model.json xgb_structural_scaler.pkl pipeline/
```

The training script uses the same 5-fold CV setup and prints updated performance metrics.
