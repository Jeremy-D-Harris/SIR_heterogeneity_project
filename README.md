# Code for: Infections are not alike: the effects of covariation between individual susceptibility and transmissibility on epidemic dynamics
Jeremy D Harris, Esther Gallmeier, Stephen J. Beckett, Joshua S. Weitz (2024)


**Code for:** "Infections are not alike: the effect of covariation between individual-level susceptibility and transmissibility on population-level epidemic dynamics." This repository contains all the MATLAB codes for both running & plotting results from SIR models. Variation is included in susceptibility ($\varepsilon$) and transmissibility ($\delta$) and compared to mdoels with variation in susceptilibility alone and the classic SIR model, i.e., homogeneous susceptibility and transmissibility. We implement the following initial joint distributions in $\varepsilon$ and $\delta$: <br>
<ol>
<li> Uncorrelated Gamma distributions;  </li> 
<li>  Low- & high-variance Gaussian distributions;  </li> 
<li>  Negative binomial distribution ([1]).  </li>
</ol>
A preprint of the manuscript can be found on BioRxiv: [DOI]()

**Instructions:**
MATLAB 2023b and 2024a was used to run model simulations and plot figures. We also use the Statistics and Machine Learning Toolbox for setting bivariate distributions. <br>

  <p>First, downloaded the Github repository. <u>To plot the figures found in the manuscript</u>, navigate to the subdirectory 'code_plt.' The name of the MATLAB script file (*.m file)  can be found in the tabulation file: 'tabulate_names_figures.docx' in the top (or main) directory. Open the appropriate script file in MATLAB. To reproduce the simulation data, navigate to the subdirectory 'code_sims' and run the appropriate function. See again, the tabulation file: <em>'tabulate_names_figures.docx'</em> in the top (or main) directory.</p> <br>

<u>See below for subfolder descriptions.</u>

**Folder descriptions:** <br>

- **Manuscript_forCodeReview:** Downloaded August 18, 2024


- **code_plt:** All code to plot model simulation results. <br>

  <p>Read in the data from 'data/' and plot the figures in the manuscript. At the top of the file, there is an option to save the figure or not. If 'save_ans'= 1, the figures will save to the folder 'figures,' a subdirectory of the top (or main) directory. These are the figures included in the manuscript.</p>


- **code_sims:** All code to simulate models. <br>

  <p> At the top of the file, there is an option to save the simulation data or not. If 'save_ans'= 1, the data will save to the folder 'data,' a subdirectory of the top (or main) directory. To run and save additional results, set each of the options equal to 1. <br> 

  <u>To simulate (nonzero) correlations, first set the options, 'save_distributions' = 1 and 'readin_init_joint = 0.'</u> This will produce a transient at the beginning of the epidemic dynamics and save the distributions at time point 'index_day_distribution' = 40. <u>Next, run the simulation again, but now setting the options, 'save_distributions' = 0 and 'readin_init_joint = 1.'</u> This will remove the transient, because the distributions have already converged to the eigendistributions of the epidemic dynamics.</p>


- **data:** Simulation data is saved here and used to plot current results. <br>


- **figures:** Figure files are saved to this folder.


**References:** <br>
[1]  Famoye, F. (2010). On the bivariate negative binomial regression model. Journal of Applied Statistics, 37(6), 969-981. <br><br>


