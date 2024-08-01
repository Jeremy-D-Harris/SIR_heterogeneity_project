% function void = main_SIRvariation_Gamma(void)

% simulate SIR model with transmissibility & susceptibility variation
% using Gamma distribution

%%
clear all; close all; clc;


%% want to save results?
save_results = 0;
% 0: don't save
% 1: save

filename_results = 'Gaussian_lowvariance_update072924.mat';

%% options
% run & save Classic SIR?
run_classic_SIR = 1;
% save_results_classic = 1;

% run & save variation in susceptibility SIR?
run_variation_susc_SIR = 1;
% save_results_variation_susc = 1;

% run & save reduced model?
run_reduced_SIR = 1;
% save_results_reduced = 1;

% note! - could get a save error, if didn't run but wanted to save
save_additional_results = 1;


% save distributions during exponential growth
% want to plot distributions at certain time?
want_to_plt_distributions = 1;
save_distributions = 1; % save distribution at certain time?
index_day_distribution = 50; % what time? (days)

% want to read in distribution from a file?
readin_init_joint = 0;
filename_distributions_load = 'GaussianIndependent_joint_expgrowth.mat';

%  = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];
my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;


%% parameters
% recovery rate
gam=1/10; % assume mean infectious period is 10 days


intended_R0 = 2;
% match basic reproduction numbers: R0=2
bet=intended_R0*gam; %Poisson-like, rho=0


% starting population
N = 1; % population size

% discrete mesh
n = 100; % number transmissibility classes
m = n; % susceptibilty classes


eps_end = 6;
eps = linspace(0,eps_end,m);
del = linspace(0,eps_end,n);
dx = eps_end/n;

eps_plt = eps;
del_plt = del;
eps(1)=dx/2;
del(1)=dx/2;

[Eps,Del] = meshgrid(eps,del);

% create parameter structure
params.bet = bet;
params.gam = gam;
params.N = N;
params.eps=eps;
params.del=del;
params.E=Eps;
params.D=Del;

% time is in days; 200 days is about 6-7 months; 250 is 8-9 months
t_start = 0; t_end = 200;
dt = 1;
t_span = t_start:dt:t_end;
params.t_span = t_span;


%% Initialize Joint Distributions & calculate Marginals

if readin_init_joint

    % read in joint distribution from file, e.g., during exponential growth

    folder_location = './sim_results/';
    filename_distributions_load = 'Gamma_joint_expgrowth_nocorr.mat';

    load(strcat(folder_location,filename_distributions_load));
    init_joint_S  = data.init_joint_S;
    init_joint_I  = data.init_joint_I;
    init_joint_R  = data.init_joint_R;
    index_day_distribution = data.index_day_distribution;

    fprintf('Initial Distributions: \n');
    fprintf('Load From: \n');
    fprintf(strcat(filename_distributions_load,'\n\n'));

else

    fprintf('Initial Distributions: \n');
    fprintf('Paremtrized Gaussian Distribution \n\n');

    %% variation in susceptibility/transmissibility
    % multivariate Normal Distribution:
    % https://www.mathworks.com/help/stats/mvnpdf.html


    %% Initialize Distributions

    % could better name 'intended'
    set_mean_eps_S = 1;
    set_mean_delta_S = 1;
    set_mean_eps_I = set_mean_eps_S;
    set_mean_delta_I = set_mean_delta_S;
    
    params.set_mean_eps_S = set_mean_eps_S;
    params.set_mean_delta_S = set_mean_delta_S;
    % params.set_mean_eps_I = set_mean_eps_I;
    % params.set_mean_delta_I = set_mean_delta_I;

    [X1, X2] = meshgrid(eps,del);
    X = [X1(:) X2(:)];
    mean_S = [set_mean_eps_S set_mean_delta_S]; %means for S
    mean_I = [set_mean_eps_I set_mean_delta_I]; %means for I

    corr_coeff = 0; % -0.8, 0, 0.8;
    set_variance_eps = 0.15;
    set_variance_delta = 0.05;

    offdiagonal_Sigma = corr_coeff*sqrt(set_variance_eps)*sqrt(set_variance_delta);

    Sigma_S = [set_variance_eps, offdiagonal_Sigma; offdiagonal_Sigma, set_variance_delta]; %change off-diagonals for correlation

    init_joint_S = mvnpdf(X, mean_S, Sigma_S);
    init_joint_S = reshape(init_joint_S,length(del),length(eps))/sum(sum(init_joint_S))/dx/dx;
    init_joint_I = init_joint_S;
    
    % marginals
    marginal_eps_S = dx*sum(init_joint_S);
    marginal_delta_S = dx*sum(init_joint_S,2)';

    marginal_eps_I = dx*sum(init_joint_I);
    marginal_delta_I = dx*sum(init_joint_I,2)';

    % Calculated Means
    mean_eps_S = dx*sum(eps.*marginal_eps_S);
    mean_delta_S = dx*sum(del.*marginal_delta_S);
    mean_eps_I = dx*sum(eps.*marginal_eps_I);
    mean_delta_I = dx*sum(del.*marginal_delta_I);

    % Calculated variance in S
    variance_eps_S = dx*sum((eps- mean_eps_S*ones(size(eps))).^2.*marginal_eps_S);
    variance_delta_S = dx*sum((del- mean_delta_S*ones(size(del))).^2.*marginal_delta_S);

    params.mean_eps_S = mean_eps_S
    params.mean_delta_S = mean_delta_S
    params.mean_delta_I = mean_delta_I;
    params.mean_eps_I = mean_eps_I;
    params.variance_eps_S = variance_eps_S;
    params.variance_delta_S = variance_delta_S;
    params.init_joint_S = init_joint_S;
    params.init_joint_I = init_joint_I;


end



init_marginal_eps_S = dx*sum(init_joint_S);
init_marginal_delta_S = dx*sum(init_joint_S,2)';

init_marginal_eps_I = dx*sum(init_joint_I);
init_marginal_delta_I = dx*sum(init_joint_I,2)';


%% Initializing Eigendirections
eps_perturb = 1.05e-4; % SIR peaks at 100 days

[eigen_direction_SIR_ed] = get_eigendirection_mean_eps_delta(params);

S_init = N + eps_perturb*eigen_direction_SIR_ed(1);
I_init = eps_perturb*eigen_direction_SIR_ed(2);
R_init = eps_perturb*eigen_direction_SIR_ed(3);

init_values_S_eps_delta = (dx)^2*init_joint_S*S_init;
params.init_values_S_eps_delta = init_values_S_eps_delta;

init_values_I_eps_delta = (dx)^2*init_joint_S*I_init;
params.init_values_I_eps_delta = init_values_I_eps_delta;

init_values_R_eps_delta = (dx)^2*init_joint_S*R_init;
params.init_values_R_eps_delta = init_values_R_eps_delta;

init_conds = [reshape(init_values_S_eps_delta, m*n,1); reshape(init_values_I_eps_delta,m*n,1); reshape(init_values_R_eps_delta, m*n,1)];

%check should equal to population size
% sum(sum((init_S_eps_delta_values + init_I_eps_delta_values + init_R_eps_delta_values),2));

% plot initial distributions - if you want!
if 1

    plt_distributions(eps_plt, del_plt, init_joint_S, init_joint_I, init_marginal_eps_S, init_marginal_delta_S, init_marginal_eps_I, init_marginal_delta_I, my_rgb_colors)

end


%% Simulate model

tic;
fprintf('Simulating... \n');

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

S_traj_eps_delta_array =zeros(length(t_span),m,n);
I_traj_eps_delta_array = zeros(length(t_span),m,n);
R_traj_eps_delta_array = zeros(length(t_span),m,n);

[t,y_traj] = ode45(@(t,y)simulate_SIR_eps_delta(t,y,params), params.t_span, init_conds, options);

% reshape from vector to matrix form
% S,I,R Trajectories - variation in eps, delta
S_traj_eps_delta_array(:,:,:) = reshape(y_traj(:,1:(m*n)), length(t_span),m,n);
I_traj_eps_delta_array(:,:,:) = reshape(y_traj(:,(m*n+1):(2*m*n)),length(t_span),m,n);
R_traj_eps_delta_array(:,:,:) = reshape(y_traj(:,(2*m*n+1):(3*m*n)),length(t_span),m,n);

% S,I,R Trajectories
S_traj = sum(S_traj_eps_delta_array,[2,3]); % population sizes
I_traj = sum(I_traj_eps_delta_array,[2,3]);
R_traj = sum(R_traj_eps_delta_array,[2,3]);

% final outbreak size
FOS_traj_eps_delta = I_traj + R_traj;

toc;
fprintf('All done with simulating! \n\n');


%% Now calculate over time
% joints, marginals, and means
for kk=1:length(t_span)

    % S(t,eps,delta)
    this_S_eps_delta_array(:,:) = S_traj_eps_delta_array(kk,:,:);

    % joint in S
    this_joint_S = this_S_eps_delta_array/S_traj(kk)/dx/dx;
    joint_S_traj(kk,:,:) = this_joint_S;

    % marginals in S
    marginal_eps_S_traj(kk,:) = dx*reshape(sum(this_joint_S,1),1,m);
    marginal_delta_S_traj(kk,:) = dx*reshape(sum(this_joint_S,2),1,n);

    % mean susceptibility
    mean_eps_S_traj(kk,1) = dx*sum(eps.*marginal_eps_S_traj(kk,:));
    mean_delta_S_traj(kk,1) = dx*sum(del.*marginal_delta_S_traj(kk,:));

    % variance in susceptibility
    variance_eps_S_traj(kk,1) = dx*sum((eps- mean_eps_S_traj(kk)*ones(size(eps))).^2.*marginal_eps_S_traj(kk,:));
    variance_delta_S_traj(kk,1) = dx*sum((del- mean_delta_S_traj(kk)*ones(size(del))).^2.*marginal_delta_S_traj(kk,:));

    % I(t,eps,delta)
    this_I_eps_delta_array(:,:) = I_traj_eps_delta_array(kk,:,:);

    % joint in I
    this_joint_I = this_I_eps_delta_array/I_traj(kk)/dx/dx;
    joint_I_traj(kk,:,:) = this_joint_I;

    % marginals in I
    marginal_eps_I_traj(kk,:) = dx*reshape(sum(this_joint_I,1),1,m);
    marginal_delta_I_traj(kk,:) = dx*reshape(sum(this_joint_I,2),1,n);

    % mean transmissibility
    mean_delta_I_traj(kk,1) = dx*sum(del.*marginal_delta_I_traj(kk,:));

    % variance transmissibility
    variance_delta_I_traj(kk,1) = dx*sum(((del- mean_delta_I_traj(kk)*ones(size(del))).^2).*marginal_delta_I_traj(kk,:));

    % I(t,eps,delta)
    this_R_eps_delta_array(:,:) = R_traj_eps_delta_array(kk,:,:);

    % joint in R
    this_joint_R(:,:) = reshape(R_traj_eps_delta_array(kk,:,:),m,n)/R_traj(kk)/dx/dx;
    joint_R_traj(kk,:,:) = this_joint_R;


end



%% get Rt - variation in eps & delta
SIR_traj_eps_delta = zeros(length(t_span),5);
SIR_traj_eps_delta(:,1) = S_traj;
SIR_traj_eps_delta(:,2) = I_traj;
SIR_traj_eps_delta(:,3) = R_traj;
SIR_traj_eps_delta(:,4) = mean_eps_S_traj;
SIR_traj_eps_delta(:,5) = mean_delta_I_traj;

Rt_traj_eps_delta = transpose(get_Rt_SIR_eps_delta(params,SIR_traj_eps_delta));
results.Rt_traj = Rt_traj_eps_delta;

% total incidence
total_incidence = bet*mean_delta_I_traj.*I_traj.*mean_eps_S_traj.*S_traj;
results.total_incidence = total_incidence;


% coefficient of variation (squared)

% CV^2 susceptibility
CV2_eps_S_traj = variance_eps_S_traj./(mean_eps_S_traj.^2);
results.CV2_eps_S_traj = CV2_eps_S_traj;

% CV^2 potential transmissibility
CV2_delta_S_traj = variance_delta_S_traj./(mean_delta_S_traj.^2);
results.CV2_delta_S_traj = CV2_delta_S_traj;

% CV^2 transmissibility
CV2_delta_I_traj = variance_delta_I_traj./(mean_delta_I_traj.^2);
results.CV2_delta_I_traj = CV2_delta_I_traj;


%% optional to run: save & run classic SIR

if run_classic_SIR

    S_traj_SIR_classic =zeros(length(t_span),1);
    I_traj_SIR_classic = zeros(length(t_span),1);
    R_traj_SIR_classic = zeros(length(t_span),1);

    init_conds_SIR_classic = [S_traj(1);I_traj(1);R_traj(1)];

    [t,y_traj_classic] = ode45(@(t,y)simulate_SIR(t,y,params), params.t_span, init_conds_SIR_classic, options);

    S_traj_SIR_classic = y_traj_classic(:,1);
    I_traj_SIR_classic = y_traj_classic(:,2);
    R_traj_SIR_classic = y_traj_classic(:,3);

    results_classic.S_traj = S_traj_SIR_classic;
    results_classic.I_traj = I_traj_SIR_classic;
    results_classic.R_traj = R_traj_SIR_classic;


    %% get Rt - classic SIR!
    SIR_traj_classic = zeros(length(t_span),3);
    SIR_traj_classic(:,1) = S_traj_SIR_classic;
    SIR_traj_classic(:,2) = I_traj_SIR_classic;
    SIR_traj_classic(:,3) = R_traj_SIR_classic;

    Rt_traj_classic = transpose(get_Rt_SIR_classic(params,SIR_traj_classic));
    results_classic.Rt_traj = Rt_traj_classic;

    % total incidence
    total_incidence_classic = bet*I_traj_SIR_classic.*S_traj_SIR_classic;
    results_classic.total_incidence = total_incidence_classic;

end

%% optional to run: save & run Susceptible Variation SIR
% % eps_pertub = 1e-4;
if run_variation_susc_SIR

    init_S_eps_values = transpose(sum(init_values_S_eps_delta,1));
    init_I_eps_values = transpose(sum(init_values_I_eps_delta,1));
    init_R_eps_values = transpose(sum(init_values_R_eps_delta,1));
    % init_R_eps_values = zeros(m,1);

    S_traj_eps_values = zeros(length(t_span),m);
    I_traj_eps_values = zeros(length(t_span),m);
    R_traj_eps_values = zeros(length(t_span),m);

    % to pass to ODE function
    init_conds_var_susc = [init_S_eps_values; init_I_eps_values; init_R_eps_values];

    % simulate variation in susceptibility
    [t,y_traj_eps] = ode45(@(t,y)simulate_SIR_eps(t,y,params), params.t_span, init_conds_var_susc, options);

    % S-I-R values: time goes down, eps goes across
    S_traj_eps_values = y_traj_eps(:,1:m);
    I_traj_eps_values = y_traj_eps(:,(m+1):(2*m));
    R_traj_eps_values = y_traj_eps(:,(2*m+1):end);

    results_var_susc.S_traj_eps_values = S_traj_eps_values;
    results_var_susc.I_traj_eps_values = I_traj_eps_values;
    results_var_susc.R_traj_eps_values = R_traj_eps_values;

    % sum values across eps variation
    S_traj_var_susc = sum(S_traj_eps_values,2);
    I_traj_var_susc = sum(I_traj_eps_values,2);
    R_traj_var_susc = sum(R_traj_eps_values,2);

    results_var_susc.S_traj = S_traj_var_susc;
    results_var_susc.I_traj = I_traj_var_susc;
    results_var_susc.R_traj = R_traj_var_susc;


    for kk=1:length(t_span)

        marginal_eps_S_traj_var_susc(kk,:) = S_traj_eps_values(kk,:)/S_traj_var_susc(kk)/dx;
        mean_epsilon_S_traj_var_susc(kk,1) = dx*sum(eps.*marginal_eps_S_traj_var_susc(kk,:));
        variance_eps_S_traj_var_susc(kk,1) = dx*sum((eps- mean_epsilon_S_traj_var_susc(kk)*ones(size(eps))).^2.*marginal_eps_S_traj_var_susc(kk,:));

    end
    
    results_var_susc.marginal_eps_S_traj = marginal_eps_S_traj_var_susc;
    results_var_susc.mean_epsilon_S_traj = mean_epsilon_S_traj_var_susc;
    results_var_susc.variance_eps_S_traj = variance_eps_S_traj_var_susc;
    
    % CV^2 susceptibility
    CV2_eps_S_traj_var_susc = variance_eps_S_traj_var_susc./(mean_epsilon_S_traj_var_susc.^2);
    results_var_susc.CV2_eps_S_traj = CV2_eps_S_traj_var_susc;



    %% get Rt - variation in eps - only!!
    SIR_traj_var_susc = zeros(length(t_span),4);
    SIR_traj_var_susc(:,1) = S_traj_var_susc;
    SIR_traj_var_susc(:,2) = I_traj_var_susc;
    SIR_traj_var_susc(:,3) = R_traj_var_susc;
    SIR_traj_var_susc(:,4) = mean_epsilon_S_traj_var_susc;

    Rt_traj_var_susc = transpose(get_Rt_SIR_eps(params,SIR_traj_var_susc));
    results_var_susc.Rt_traj = Rt_traj_var_susc;

    % total incidence
    total_incidence_var_susc = bet*I_traj_var_susc.*mean_epsilon_S_traj_var_susc.*S_traj_var_susc;
    results_var_susc.total_incidence = total_incidence_var_susc;

end


if run_reduced_SIR

    S_traj_SIR_reduced =zeros(length(t_span),1);
    I_traj_SIR_reduced = zeros(length(t_span),1);
    R_traj_SIR_reduced = zeros(length(t_span),1);
    mean_eps_traj_SIR_reduced = zeros(length(t_span),1);

    init_conds_SIR_reduced = [S_traj(1);I_traj(1);R_traj(1);mean_eps_S_traj(1)];

    [t,y_traj_reduced] = ode45(@(t,y)simulate_SIR_reducedgaussian(t,y,params), params.t_span, init_conds_SIR_reduced, options);

    S_traj_SIR_reduced = y_traj_reduced(:,1);
    I_traj_SIR_reduced = y_traj_reduced(:,2);
    R_traj_SIR_reduced = y_traj_reduced(:,3);
    mean_eps_traj_SIR_reduced = y_traj_reduced(:,4);

    results_reduced.S_traj = S_traj_SIR_reduced;
    results_reduced.I_traj = I_traj_SIR_reduced;
    results_reduced.R_traj = R_traj_SIR_reduced;


    %% get Rt - reduced model
    SIR_traj_reduced = zeros(length(t_span),3);
    SIR_traj_reduced(:,1) = S_traj_SIR_reduced;
    SIR_traj_reduced(:,2) = I_traj_SIR_reduced;
    SIR_traj_reduced(:,3) = R_traj_SIR_reduced;
    SIR_traj_reduced(:,4) = mean_eps_traj_SIR_reduced;

    Rt_traj_reduced = transpose(get_Rt_SIR_reducedgaussian(params,SIR_traj_reduced));
    results_reduced.Rt_traj = Rt_traj_reduced;

    % total incidence
    total_incidence_reduced = bet*I_traj_SIR_reduced.*S_traj_SIR_reduced.*mean_eps_traj_SIR_reduced;
    results_reduced.total_incidence = total_incidence_reduced;

end



%% Plotting
f2=figure(2); set(f2, 'Position', [900   50   400   930]);

subplot(3,1,1);

if run_classic_SIR

    plot(params.t_span, S_traj_SIR_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;
    plot(params.t_span, I_traj_SIR_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;
    plot(params.t_span, R_traj_SIR_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;

end


q(1)=plot(params.t_span, S_traj,'Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
q(2)=plot(params.t_span, I_traj,'Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
q(3)=plot(params.t_span, R_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;

% if run_variation_susc_SIR
% 
%     plot(params.t_span, S_traj_var_susc,'k--','LineWidth',2); hold on;
%     plot(params.t_span, I_traj_var_susc,'k--','LineWidth',2); hold on;
%     plot(params.t_span, R_traj_var_susc,'k--','LineWidth',2); hold on;
% 
% end

if run_reduced_SIR

    plot(params.t_span, S_traj_SIR_reduced,'k--','LineWidth',2); hold on;
    plot(params.t_span, I_traj_SIR_reduced,'k--','LineWidth',2); hold on;
    plot(params.t_span, R_traj_SIR_reduced,'k--','LineWidth',2); hold on;

end

axis([0 t_end 0 1.1]);
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
title('Dynamics')
legend(q,{'S','I','R'},'Location','SouthWest');
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

subplot(3,1,2);
r(2)=plot(params.t_span, mean_delta_I_traj,'-','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
r(1)=plot(params.t_span, mean_eps_S_traj,'-','Color',my_rgb_colors(3,:),'LineWidth',2); hold on;


axis([0 t_end 0.2 1.2]);
% ylim([0 2])
title('Mean Trajectories');
xlabel('Time (days)'); %ylabel({'Population'; 'Fraction'});
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

legend(r,{'$\bar{\varepsilon}$','$\bar{\delta}_I$'},'Interpreter','Latex','Location','SouthWest');

subplot(3,1,3);

if run_classic_SIR

    plot(params.t_span, Rt_traj_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;

end

plot(params.t_span, Rt_traj_eps_delta,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;

axis([0 t_end 0.2 2.2]);

xlabel('Time (days)'); ylabel([{'Effective'; 'Reproduction'; 'Number'}]);
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');


%% plot joint in S and I at certain time point
if want_to_plt_distributions

    this_ind = index_day_distribution;

    this_joint_S(:,:) = joint_S_traj(this_ind,:,:);
    this_joint_I(:,:) = joint_I_traj(this_ind,:,:);
    this_joint_R(:,:) = joint_R_traj(this_ind,:,:);


    this_marg_eps_S = dx*sum(this_joint_S);
    this_marg_delta_S = dx*sum(this_joint_S,2)';

    this_marg_eps_I = dx*sum(this_joint_I);
    this_marg_delta_I = dx*sum(this_joint_I,2)';

    plt_distributions(eps_plt, del_plt, this_joint_S, this_joint_I, this_marg_eps_S, this_marg_delta_S, this_marg_eps_I, this_marg_delta_I, my_rgb_colors)

    % want to save this distribution?
    if save_distributions

        % save the joints to a data structure
        data.init_joint_S = this_joint_S;
        data.init_joint_I = this_joint_I;
        data.init_joint_R = this_joint_R;
        data.index_day_distribution = index_day_distribution;

        folder_location = './sim_results/';

        save(strcat(folder_location,filename_distributions_load),'data');

        fprintf('Saved Distribution to File: \n');
        fprintf(strcat(filename_distributions_load,'\n\n'));

    else

        fprintf('Distribution Not Saved. \n\n');

    end




end




%% collect results
% SIR - eps & delta

% initial distributions
results.init_joint_S = init_joint_S;
results.init_joint_I = init_joint_I;

results.init_marginal_eps_S = init_marginal_eps_S;
results.init_marginal_delta_S = init_marginal_delta_S;
results.init_marginal_eps_I = init_marginal_eps_I;
results.init_marginal_delta_I = init_marginal_delta_I;

% trajectories
results.S_traj = S_traj;
results.S_traj_eps_delta_array = S_traj_eps_delta_array;

results.I_traj = I_traj;
results.I_traj_eps_delta_array = I_traj_eps_delta_array;

results.R_traj = R_traj;
results.R_traj_eps_delta_array = R_traj_eps_delta_array;

% final outbreak size
results.FOS_traj = FOS_traj_eps_delta;

% marginals over time
results.marginal_eps_S_traj = marginal_eps_S_traj; % susceptibility in S
results.marginal_delta_S_traj = marginal_delta_S_traj;
results.marginal_eps_I_traj = marginal_eps_I_traj;
results.marginal_delta_I_traj = marginal_delta_I_traj; % transmissibility in I

% mean value trajectories
results.mean_eps_S_traj = mean_eps_S_traj;
results.mean_delta_S_traj = mean_delta_S_traj;
results.mean_delta_I_traj = mean_delta_I_traj;

% variance trajectories
results.variance_eps_S_traj = variance_eps_S_traj;
results.variance_delta_I_traj = variance_delta_I_traj;

% CV^2 susceptibility
results.CV2_eps_S_traj = CV2_eps_S_traj;

% CV^2 potential transmissibility
results.CV2_delta_S_traj = CV2_delta_S_traj;

% CV^2 transmissibility
results.CV2_delta_I_traj = CV2_delta_I_traj;



%%
% save simulated results
if save_results==1

    folder_location = './../data/';

    if save_additional_results

        % + classic + variation susceptibility + reduced model
        save(strcat(folder_location,filename_results),'params','results','results_classic','results_var_susc','results_reduced');


    else

        % save variation in eps & delta, exclusively
        save(strcat(folder_location,filename_results),'params','results');

    end

    fprintf('Saved Results to File: \n');
    fprintf(strcat(filename_results,'\n'));

else

    fprintf('Results Not Saved. \n');

end




