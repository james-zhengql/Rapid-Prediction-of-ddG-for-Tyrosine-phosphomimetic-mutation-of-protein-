# Prediction of ddG for Tyrosine phosphomimetic mutation of protein 

We used FoldX, EvoEF and structural information to predict the ddG of phosphormimetic protein. We trained the XGboost on the tsuboyama et al. dataset.

## Programs

- main.py: predicting ddG for mutations, can be run directly follow the instructions inside.
- FoldX_stability.py: Use FoldX to calculate ddG given a csv file of protein and mutations
- EvoEF_stability.py: Use EvoEF to calculate ddG given a csv file of protein and mutations
- Pymol_feature.py: Use Pymol to extract protein residue features like solvent accessible surface are, number of residue in contact etc.


## Prerequisites

Before you begin, ensure you have met the following requirements:
- FoldX https://foldxsuite.crg.es/
- EvoEF https://zhanggroup.org/EvoEF/
- Pymol https://www.pymol.org/


