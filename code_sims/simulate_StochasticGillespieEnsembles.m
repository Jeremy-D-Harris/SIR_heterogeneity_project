%code to load intial joint distributions, perform ensemble stochastic
%simulations and save data.

clear all; close all; clc;

%set a random seed for rng
seed = 1321910;
rng(seed)
% set population and number of outbreaks to simulate for each type
PopulationSize = 10000;
Trajectories = 1000;
% set up experiments
stochSIR_sir.description = "SIR model, no variation";
stochSIR_sir.Population_Size = PopulationSize;
stochSIR_uncorr.description = "uncorrelated susceptibility and transmissibility";
stochSIR_uncorr.Population_Size = PopulationSize;
stochSIR_poscorr.description = "positive correlation between susceptibility and transmissibility";
stochSIR_poscorr.Population_Size = PopulationSize;
stochSIR_negcorr.description = "negative correlation between susceptibility and transmissibility";
stochSIR_negcorr.Population_Size = PopulationSize;

%set epi parameters for the SIR  (others are loaded from file)
stochSIR_sir.beta = 0.2;
stochSIR_sir.gamma = 0.1;
stochSIR_sir=RunStochasticGillepsieEnsemble(stochSIR_sir,[],Trajectories,false);
disp("stochSIR_sir complete")

%load and run other cases
file_location = '../data/';
stochSIR_uncorr=RunStochasticGillepsieEnsemble(stochSIR_uncorr,append(file_location,'GaussianNoCorrelation.mat'),Trajectories,false);
disp("stochSIR_uncorr complete")
stochSIR_poscorr=RunStochasticGillepsieEnsemble(stochSIR_poscorr,append(file_location,'GaussianPositiveCorrelation_0pt6.mat'),Trajectories,false);
disp("stochSIR_poscorr complete")
stochSIR_negcorr=RunStochasticGillepsieEnsemble(stochSIR_negcorr,append(file_location,'GaussianNegativeCorrelation_0pt6.mat'),Trajectories,false);
disp("stochSIR_negcorr complete")

%save output to file
file_location = '../data/';
timenow = replace(string(datetime()),{' ',':','-'},'');
save(append(file_location,"stochastic_output_",timenow,".mat"),"stochSIR_sir","stochSIR_uncorr","stochSIR_poscorr","stochSIR_negcorr");