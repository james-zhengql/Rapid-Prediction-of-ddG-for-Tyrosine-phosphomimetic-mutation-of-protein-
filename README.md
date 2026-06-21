# Rapid Prediction of ddG for Tyrosine Phosphomimetic Mutation of Protein

Predicts the stability change (ddG) caused by phosphomimetic mutations at tyrosine residues. Two approaches are provided: a self-contained structural pipeline (recommended) and the original energy-function pipeline using FoldX and EvoEF.

---

## Approaches

### Structural Pipeline (Recommended)

Located in [`pipeline/`](pipeline/). Requires only a PDB file — no external tools (no FoldX, EvoEF, or PyMOL). Features are computed purely from structure using BioPython and pydssp, then passed to a trained XGBoost regressor.

**5-fold cross-validation on the Tsutsumi dataset (tyrosine entries):**

| Metric   | Score |
|----------|-------|
| Pearson  | 0.626 |
| Spearman | 0.611 |
| RMSE     | 1.110 |
| MSE      | 1.232 |

See [`pipeline/README.md`](pipeline/README.md) for setup and usage.

### Original Energy-Function Pipeline

Uses FoldX and EvoEF stability calculations combined with PyMOL-extracted structural features. Requires external installations of FoldX, EvoEF, and PyMOL.

| File                | Description                                              |
|---------------------|----------------------------------------------------------|
| `main.py`           | Top-level prediction script                              |
| `Foldx_stability.py`| Calculate ddG using FoldX given a CSV of proteins/mutations |
| `EvoEF_stability.py`| Calculate ddG using EvoEF given a CSV of proteins/mutations |
| `pymol_feature.py`  | Extract residue-level features (SASA, contacts, etc.) via PyMOL |

---

## Repository Structure

```
├── pipeline/                    # Self-contained structural pipeline (recommended)
│   ├── pipeline.py              # Main prediction script
│   ├── xgb_structural_model.json# Trained XGBoost model
│   ├── xgb_structural_scaler.pkl# Fitted StandardScaler
│   ├── pdb_files/               # Place input PDB files here
│   │   └── 1A22.pdb             # Example
│   └── README.md                # Pipeline-specific docs
│
├── main.py                      # Original pipeline entry point
├── Foldx_stability.py
├── EvoEF_stability.py
└── pymol_feature.py
```

---

## Quick Start (Structural Pipeline)

```bash
cd pipeline/
python3 -m venv venv && source venv/bin/activate
pip install biopython pandas numpy scikit-learn xgboost joblib pydssp

# Run on example PDB
python pipeline.py

# Run on your own PDB
python pipeline.py --pdb /path/to/your.pdb
```

Results are written to `pipeline/output/predictions.csv`.

---

## Prerequisites (Original Pipeline)

- [FoldX](https://foldxsuite.crg.eu/)
- [EvoEF](https://zhanggroup.org/EvoEF/)
- [PyMOL](https://www.pymol.org/)
