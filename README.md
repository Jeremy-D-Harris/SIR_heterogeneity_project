# *Code for:* "Infections are not alike: the effects of covariation between individual susceptibility and transmissibility on epidemic dynamics"
Jeremy D Harris, Esther Gallmeier, Jonathan Dushoff, Stephen J. Beckett, Joshua S. Weitz (2024)


This repository contains all the MATLAB codes for both running & plotting results from SIR-like models incorporating both susceptibility and transmissibility variation. Variation is included in susceptibility ($\varepsilon$) and transmissibility ($\delta$) and compared to models with variation in susceptilibility alone and the standard SIR model, i.e., homogeneous susceptibility and transmissibility. We implement the following initial joint distributions in $\varepsilon$ and $\delta$: <br>
<ol>
<li> Uncorrelated gamma distributions;  </li> 
<li>  Low- & high-variance truncated Gaussian distributions.  </li> 
</ol>
A preprint of the manuscript can be found on medRxiv. Please cite our work as:


> Harris J.D., Gallmeier E., Dushoff J., Beckett S.J., Weitz J.S. 2024. Infections are not alike: the effects of covariation between individual susceptibility and transmissibility on epidemic dynamics. medRxiv aaaa doi: https://doi.org/10.1101/2024.bbbbbbbb 

<br>

**Instructions:** <br>
MATLAB 2023b and 2024a was used to run model simulations and plot figures. We also use the Statistics and Machine Learning Toolbox for setting bivariate distributions. <br>

  <p>First, downloaded the Github repository. <u>To plot the figures found in the manuscript</u>, navigate to the subdirectory 'code_plt.' The name of the MATLAB script file (*.m file)  can be found in the tabulation file: 'tabulatefigures.md' in the main directory. Open the appropriate script file in MATLAB. To reproduce the simulation data, navigate to the subdirectory 'code_sims' and run the appropriate function. See again, the tabulation file: <em>'tabulatefigures.md'</em> in the main directory.</p> <br>

<u>See below for subfolder descriptions.</u>

**Folder descriptions:** <br>

- **code_plt:** All code to plot model simulation results. <br>

  <p>Read in the data from 'data/' and plot the figures in the manuscript. At the top of these files, there is an option to save the figure or not. If 'save_ans'= 1, the figures will save to the folder 'figures,' a subdirectory of the main directory. These are the figures included in the manuscript.</p>


- **code_sims:** All code to simulate models. <br>

  <p> At the top of these files, there is an option to save the simulation data or not. If 'save_ans'= 1, the data will save to the folder 'data,' a subdirectory of the main directory. To run and save additional results, set each of the options equal to 1. <br> 

  <u>To simulate (nonzero) correlations, first set the options, 'save_distributions' = 1 and 'readin_init_joint = 0.'</u> This will produce a transient at the beginning of the epidemic dynamics and save the distributions at time point 'index_day_distribution' = 40. <u>Next, run the simulation again, but now setting the options, 'save_distributions' = 0 and 'readin_init_joint = 1.'</u> This will remove the transient, because the distributions have already converged to the eigendistributions of the epidemic dynamics.</p>


- **data:** Simulation data is saved here and used by plotting codes. <br>


- **figures:** Figure files are saved to this folder.