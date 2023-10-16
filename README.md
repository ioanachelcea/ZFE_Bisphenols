Toxicokinetic model for bisphenols in zebrafish embryo. Knime workflow and codes were written by Ioana Chelcea, Publication Title: "Physiology-informed toxicokinetic model for the zebrafish embryo test developed for bisphenols" by Ioana Chelcea, Carolina Vogs, Timo Hamers, Jacco Koekkoek, Jessica Legradi, Maria Sapounidou, Stefan Örn,  Patrik L. Andersson. Published in Chemosphere 2023. All codes are under the Creative commons licence provided in the repository.

Repository contains 3 files: 1. Knime workflow containing models and predictions (Chelcea_2023_ZFE_git.knwf), 2. Excel file containing model parameters (ZFE_parms_KNIME.xlsx), 3. R script containing Bayesian Inference example used for model parameter calibration (ZFE_Bayes_10k_widepriors.R)

IMPORTANT: Place both KNIME workflow and Excel files in the same folder and set the local KNIME directory to that folder.

Instructions for running the KNIME workflow:

Download KKNIME workflow and Excel files in the same folder.
Open KNIME and set your local directory to that folder. (You can change your local directory under File>Switch workspaces> Browse for the folder you have the files in)
To open workflow it has to be imported. Go to File>Import KNIME workflow> Select file (Browse and select the Chelcea_2023_ZFE_git.knwf file)

Trouble shooting: Packages not found even after installing them: Check which version of R is set to use in KNIME. Using a local computer version is more suitable as packages are installed there by default. To change R version go to File>Preferences>KNIME>R> Path to R Home.
