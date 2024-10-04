
clear all; close all; clc;

save_fig_ans = 1;
% save figure:
% 0 = no, 1 = yes

figure_name = 'SuppFigure_correlations_initjoints_073124';

%create color gradiets
c2 = [133,192,249]/255; % light blue
c1 = [15,32,128]/255; % dark blue

depth = 9;
grad1=colorGradient(c2,c1,depth);
colors_rgb = grad1;


%% load results from file
file_location = '../data/';


%% (1) positive
this_infile = 'GaussianPositiveCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));

init_joint_S(1,:,:) = params.init_joint_S;
init_joint_I(1,:,:) = params.init_joint_I;
corr_coeff(1) = params.corr_coeff;


%% (2) positive
this_infile = 'GaussianPositiveCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

init_joint_S(2,:,:) = params.init_joint_S;
init_joint_I(2,:,:) = params.init_joint_I;
corr_coeff(2) = params.corr_coeff;


%% (3) independent
this_infile = 'GaussianNoCorrelation.mat';
load(strcat(file_location,this_infile));

% should update?
init_joint_S(3,:,:) = params.init_joint_S;
init_joint_I(3,:,:) = params.init_joint_I;
corr_coeff(3) = params.corr_coeff;


%% (4) negative
this_infile = 'GaussianNegativeCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

init_joint_S(4,:,:) = params.init_joint_S;
init_joint_I(4,:,:) = params.init_joint_I;
corr_coeff(4) = params.corr_coeff;


%% (5) negative
this_infile = 'GaussianNegativeCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));

init_joint_S(5,:,:) = params.init_joint_S;
init_joint_I(5,:,:) = params.init_joint_I;
corr_coeff(5) = params.corr_coeff;


%% Plotting
X = get(0,'ScreenPixelsPerInch'); %determine screen pixels per inch (96 on windows, 72 on mac os)
factor = X/72;
f1 = figure(1); set(f1, 'Position', [100 500 factor*1050 factor*450]);


%% plot disstributions for:
% rho = -0.6, -0.3, 0, 0.3, 0.6
for count = 1:5

    this_joint_S(:,:) = init_joint_S(count,:,:);
    this_joint_I(:,:) = init_joint_I(count,:,:);

    plt_initdistributions(params.eps,params.del,this_joint_S,this_joint_I, corr_coeff(count), count)
end


%% save figure
if save_fig_ans==1

    figures_location = '../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));

else

    fprintf('Figure not saved...\n'); % want to be close to 25 days in
end