# SIR_heterogeneity_project
# Infections are not all alike: the effect of covariation between individual-level susceptibility and transmissibility on population-level epidemic dynamics
Jeremy D Harris, Esther Gallmeier, Stephen J. Beckett, Joshua S. Weitz (2024)

 -- updated 02/01/24 by JDH --

**Code for:** "Infections are not alike: the effect of covariation between individual-level susceptibility and transmissibility on population-level epidemic dynamics." This repository contains all the MATLAB codes for running SIR models with heterogeneity in susceptibility and transmissibility using the following distributions: <br>
    (1) uncorrelated Gamma distributions with means at 1; <br>
    (2) Truncated Gaussian distributions to include correlations between transmission and susceptibility; <br>
    (3) Negative binomial distribution with postive correlations.

A preprint of the manuscript can be found on BioRxiv: [DOI]()

**Instructions:**
MATLAB was used to run model simulations and plot figures. Once the Github repository is downloaded, navigate to the subdirectory 'Code_plt' to plot the figures in the manuscript by running the appropriate function. Simulation data can be reproduced by navigating to 'Code_sims' and running appropriate functions. <u>See below for subfolder descriptions.</u>

**Folder descriptions:** <br>

- **code_sims:** All code to simulate models; subfolders organized by model with models (1)-(3) described above: <br>
  (1)  `codes_from_Esther_10102021'  <br>
  (2) `code_nb_clean_Dec2022': Jeremy cleaned up negative binomial codes




Within each of these folders, the main files simulate the model and have several user options at the top of the scripts. For instance, in the the first choice for the user is to save the simulation data using the variable 'save_ans': 0 means don't save and 1 means save. The output file will be saved to the directory 'Code_plt/sim_data/' so that the corresponding figure can be produced.

- **code_plt:**
Read in the data from 'data/' and plot the figures in the manuscript. If 'save_ans'= 1, the figures will save to the folder 'figures.' From here, they are uploaded to the Overleaf document.

- **Manuscript_forCodeReview:** not yet!


- **Data** contains simulation data from Esther to plot current results organized by the distribution used: <br>

    (1) 'GammaInd': uncorrelated Gamma distributions <br>
    (2) 'GaussianInd': uncorrelated Gaussian distributions <br>
    (3) 'forward' (`speeds' not matched): correlated Gamma, Gaussian, Negative Binomial distributions  <br>
    (4) 'gaussianTruncatedVar1': (truncated) correlated Gaussian distributions <br>
    (5) 'negbinomial': positive correlated negative binomial (`speeds' matched) <br><br>


- **figures_movies** create movies that show the evolution of the joint distributions of susceptibility and transmissibility in both susceptible and infected populations

**References:** <br>
[1]  <br>
[2]  <br>
[3] <br>
[4]  <br>
