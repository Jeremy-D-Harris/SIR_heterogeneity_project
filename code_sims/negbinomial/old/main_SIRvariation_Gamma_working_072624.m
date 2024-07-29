% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation
% using negative binomial distributions

%%
clear all; close all; clc;


%% want to save results?
save_results = 0;
% 0: don't save
% 1: save

% filename = 'negbinomial_Poissonlike.mat';
% filename_results = 'negbinomial_overdispersed_nocorr_update072524.mat';
filename_results = 'Gamma_update072524.mat';

%% options
% run & save Classic SIR
run_classic_SIR = 1;
save_results_classic = 1;

% run & save variation in susceptibility SIR
run_variation_susc_SIR = 1;

% save distributions during exponential growth
% want to plot distributions at certain time?
want_to_plt_distributions = 1;
save_distributions = 0; % save distribution at certain time?
index_day_distribution = 50; % what time? (days)
filename_distributions_load = 'Gamma_joint_expgrowth.mat';

readin_init_joint = 0;

%  = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];
my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;


%% parameters
% recovery rate
gam=1/10; % assume mean infectious period is 10 days



intended_R0 = 2;
% match basic reproduction numbers: R0=3
bet=intended_R0*gam; %Poisson-like, rho=0

% want means of marginals to be 1
% mean_joint = 1;

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
eps(1)=dx/2;
del(1)=dx/2;

[Eps,Del] = meshgrid(eps,del);

% create parameter structure
params.bet = bet;
params.gam = gam;
% params.mean = mean_joint;
params.N = N;
params.eps=eps;
params.del=del;
params.E=Eps;
params.D=Del;

% time is in days; 200 days is about 6-7 months; 250 is 8-9 months
t_start = 0; t_end = 300;
dt = 1;
t_span = t_start:dt:t_end;
params.t_span = t_span;


%% Initialize Joint Distributions & calculate Marginals

if readin_init_joint

    % read in joint distribution (e.g., during exponential growth)

    folder_location = './sim_results/';
    % filename_joint_expgrowth = 'negbinomial_joint_expgrowth_poscorr.mat';
    % filename_load = 'negbinomial_joint_expgrowth_nocorr.mat';

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
    fprintf('Paremtrized Gamma Distribution \n\n');

    %% variation in susceptibility/transmissibility
    % using Gamma distributions
    % see: https://www.mathworks.com/help/stats/gampdf.html
    % y = gampdf(x,a,b) = (1/b)^a x^(a-1)/(n-1)! exp(-x/b)
    % mean = a*b
    % variance = a*b^2
    % Let a be the shape parameter, and then:
    % substitute: b = eps_bar/a
    % -> mean = eps_bar; variance = eps_bar^2/a
    % -> CV^2 = 1/a -> CV = 1/sqrt(a)

    % shape parameters
    shape_eps = 3;
    shape_delta = 10;

    params.shape_delta = shape_delta;
    params.shape_eps = shape_eps;

    init_mean_eps_S = 1;
    init_mean_delta_S = 1;
    init_mean_eps_I = 1;
    init_mean_delta_I = 1;


    %% Initialize Distributions
    init_marg_eps_S = gampdf(eps,shape_eps,init_mean_eps_S/shape_eps);
    % init_marg_eps_S = init_marg_eps_S/sum(init_marg_eps_S)/dx;
    
    init_marg_delta_S = gampdf(del,shape_delta,init_mean_delta_S/shape_delta);
    % init_marg_delta_S = init_marg_delta_S/sum(init_marg_delta_S)/dx;


    init_marg_eps_I = gampdf(eps,shape_eps,init_mean_eps_I/shape_eps);
    % init_marg_eps_I = init_marg_eps_I/sum(init_marg_eps_I)/dx;
    

    init_marg_delta_I = gampdf(del,shape_delta,init_mean_delta_I/shape_delta);
    % init_marg_delta_I = init_marg_delta_I/sum(init_marg_delta_I)/dx;


    %Independent Joint Distributions - Product of Marginals
    init_joint_S = init_marg_eps_S'*init_marg_delta_S; 
    init_joint_I = init_marg_eps_I'*init_marg_delta_I;
    init_joint_R = init_marg_eps_S'*init_marg_delta_S;

    % %Independent Joint Distributions - Product of Marginals
    % init_joint_S = init_marg_eps_S'*init_marg_delta_S;
    % init_joint_S = init_joint_S/sum(sum(init_joint_S));
    % init_joint_I = init_joint_S/sum(sum(init_joint_S));
    % init_joint_R = init_joint_S/sum(sum(init_joint_S));

end


init_marg_eps_S = dx*sum(init_joint_S);
init_marg_delta_S = dx*sum(init_joint_S,2)';

init_marg_eps_I = dx*sum(init_joint_I);
init_marg_delta_I = dx*sum(init_joint_I,2)';


% calculated mean values
init_mean_eps_S = dx*sum(eps.*init_marg_eps_S);
init_mean_delta_S = dx*sum(del.*init_marg_delta_S);
init_mean_eps_I = dx*sum(eps.*init_marg_eps_I);
init_mean_delta_I = dx*sum(del.*init_marg_delta_I);

%calculated values
init_covariance = (dx^2)*sum(sum(((eps' - init_mean_eps_S*ones(size(eps')))*(del- init_mean_delta_S*ones(size(del)))).*init_joint_S));

init_variance_eps = dx*sum((eps - init_mean_eps_S*ones(size(eps))).^2.*init_marg_eps_S);
init_sd_eps = sqrt(init_variance_eps);

init_variance_delta = dx*sum((del- init_mean_delta_S*ones(size(del))).^2.*init_marg_delta_S);
init_sd_delta = sqrt(init_variance_delta);

corr_coef = init_covariance/(init_sd_eps*init_sd_delta);


%% Initial conditions
% small perturbation in direction of eigenvector
eps_perturb = 1.25e-6;

params.eps_perturb = eps_perturb;

params.mu_eps_S = init_mean_eps_S;
params.mu_delta_I = init_mean_delta_I;

%% Initializing Eigendirections
% eps = 1.25e-6; %change as needed

% eigen_direction_SIR = get_eigendirection_SIRdelta_eps(params);

% eigen_direction_SIR = get_eigendirection_SIRdelta_eps(params);

% S_init = N - eps_perturb*abs(eigen_direction_SIR(1));
% I_init = eps_perturb*abs(eigen_direction_SIR(2));
% R_init = eps_perturb*abs(eigen_direction_SIR(3));

S_init = N - eps_perturb;
I_init = eps_perturb;
R_init = 0;

% pause;

%
init_S_eps_delta_values = init_joint_S.*(S_init); %joint probability matrix S_{i,j}
params.init_S_eps_delta_values = init_S_eps_delta_values;

init_I_eps_delta_values = init_joint_I.*(I_init);
params.init_I_eps_delta_values = init_I_eps_delta_values;

init_R_eps_delta_values = init_joint_R.*(R_init);
params.init_R_eps_delta_values = init_R_eps_delta_values;

%check should equal to population size
% sum(sum((init_S_eps_delta_values + init_I_eps_delta_values + init_R_eps_delta_values),2));

% plot initial distributions - if you want!
if 1

    plt_distributions(eps_plt, del, init_joint_S, init_joint_I, init_marg_eps_S, init_marg_delta_S, init_marg_eps_I, init_marg_delta_I, my_rgb_colors)

end


% need to reshape from matrix to vector from in order to pass to ODE function
init_conds = [reshape(init_S_eps_delta_values, m*n,1); reshape(init_I_eps_delta_values,m*n,1); reshape(init_R_eps_delta_values,m*n,1)];


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
FOS_traj_eps_delta = I_traj + R_traj;

marginal_eps_S_traj = dx*reshape(sum(S_traj_eps_delta_array,2),length(t_span),m)./S_traj;
marginal_delta_S_traj = dx*reshape(sum(S_traj_eps_delta_array,3),length(t_span),n)./S_traj;
marginal_eps_I_traj = dx*reshape(sum(I_traj_eps_delta_array,2),length(t_span),m)./I_traj;
marginal_delta_I_traj = dx*reshape(sum(I_traj_eps_delta_array,3),length(t_span),n)./I_traj;
marginal_eps_R_traj = dx*reshape(sum(R_traj_eps_delta_array,2),length(t_span),m)./R_traj;
marginal_delta_R_traj = dx*reshape(sum(R_traj_eps_delta_array,3),length(t_span),n)./R_traj;

for kk=1:length(t_span)
    % mean epsilon_S(t), delta_I(t)
    mean_epsilon_S_traj(kk) = dx*sum(eps.*marginal_eps_S_traj(kk,:));
    mean_delta_I_traj(kk) = dx*sum(del.*marginal_delta_I_traj(kk,:));

end

%% get Rt - variation in eps & delta
SIR_traj_eps_delta = zeros(length(t_span),5);
SIR_traj_eps_delta(:,1) = S_traj;
SIR_traj_eps_delta(:,2) = I_traj;
SIR_traj_eps_delta(:,3) = R_traj;
SIR_traj_eps_delta(:,4) = mean_epsilon_S_traj;
SIR_traj_eps_delta(:,5) = mean_delta_I_traj;

Rt_traj_eps_delta = get_Rt_SIR_eps_delta(params,SIR_traj_eps_delta);

%%

for kk = 1:length(t_span)

    %     marginal_eps_S_traj(kk,:) = marginal_eps_S_traj(kk,:)/sum(marginal_eps_S_traj(kk,:));
    %     marginal_delta_S_traj(kk,:) = marginal_delta_S_traj(kk,:)/sum(marginal_delta_S_traj(kk,:));
    %
    %     marginal_eps_I_traj(kk,:) = marginal_eps_I_traj(kk,:)/sum(marginal_eps_I_traj(kk,:));
    %     marginal_delta_I_traj(kk,:) = marginal_delta_I_traj(kk,:)/sum(marginal_delta_I_traj(kk,:));

    variance_eps_S_traj(kk) = sum((eps- mean_epsilon_S_traj(kk)*ones(size(eps))).^2.*marginal_eps_S_traj(kk,:));
    variance_delta_I_traj(kk) = sum(((del- mean_delta_I_traj(kk)*ones(size(del))).^2).*marginal_delta_I_traj(kk,:));

    sd_S_eps_traj(kk) = sqrt(variance_eps_S_traj(kk));
    sd_I_delta_traj(kk) = sqrt(variance_delta_I_traj(kk));

end


CV2_Seps_traj = variance_eps_S_traj./(mean_epsilon_S_traj.^2);
CV2_Idelta_traj = variance_delta_I_traj./(mean_delta_I_traj.^2);

toc;
fprintf('All done with simulating! \n\n');



%% optional to run: save & run classic SIR

if run_classic_SIR

    params_classic = params;

    S_traj_SIR_classic =zeros(length(t_span),1);
    I_traj_SIR_classic = zeros(length(t_span),1);
    R_traj_SIR_classic = zeros(length(t_span),1);

    init_conds_SIR_classic = [S_traj(1);I_traj(1);R_traj(1)];

    [t,y_traj_classic] = ode45(@(t,y)simulate_SIR(t,y,params_classic), params_classic.t_span, init_conds_SIR_classic, options);

    S_traj_SIR_classic = y_traj_classic(:,1);
    I_traj_SIR_classic = y_traj_classic(:,2);
    R_traj_SIR_classic = y_traj_classic(:,3);

    results_classic.S_traj_SIR_classic = S_traj_SIR_classic;
    results_classic.I_traj_SIR_classic = I_traj_SIR_classic;
    results_classic.R_traj_SIR_classic = R_traj_SIR_classic;


    %% get Rt - variation in eps & delta
    SIR_traj_classic = zeros(length(t_span),3);
    SIR_traj_classic(:,1) = S_traj;
    SIR_traj_classic(:,2) = I_traj;
    SIR_traj_classic(:,3) = R_traj;

    Rt_traj_classic = get_Rt_SIR_classic(params_classic,SIR_traj_classic);

    results_classic.Rt_traj_classic = Rt_traj_classic;

end

%% optional to run: save & run Susceptible Variation SIR

if run_variation_susc_SIR

    init_S_eps_values = marginal_eps_S_traj(1,:).*S_init;
    params.init_S_eps_values = init_S_eps_values;

    init_I_eps_values = marginal_eps_I_traj(1,:).*I_init;
    params.init_I_eps_values = init_I_eps_values;

    init_R_eps_values = marginal_eps_R_traj(1,:).*R_init;
    params.init_R_eps_values = init_R_eps_values;


    S_traj_eps_values =zeros(length(t_span),(m+1));
    I_traj_eps_values = zeros(length(t_span),(m+1));
    R_traj_eps_values = zeros(length(t_span),(m+1));

    % to pass to ODE function
    init_conds_var_susc = [transpose(init_S_eps_values); transpose(init_I_eps_values); transpose(init_R_eps_values)];

    % simulate variation in susceptibility
    [t,y_traj_eps] = ode45(@(t,y)simulate_SIR_eps(t,y,params), params.t_span, init_conds_var_susc, options);

    % S-I-R values: time goes down, eps goes across
    S_traj_eps_values = y_traj_eps(:,1:(m+1));
    I_traj_eps_values = y_traj_eps(:,(m+2):(2*(m+1)));
    R_traj_eps_values = y_traj_eps(:,(2*(m+1)+1):end);

    results_var_susc.S_traj_SIR_eps = S_traj_eps_values;
    results_var_susc.I_traj_SIR_eps = I_traj_eps_values;
    results_var_susc.R_traj_SIR_eps = R_traj_eps_values;

    % sum values across eps variation
    S_traj_var_susc = sum(S_traj_eps_values,2);
    I_traj_var_susc = sum(I_traj_eps_values,2);
    R_traj_var_susc = sum(R_traj_eps_values,2);

    results_var_susc.S_traj_var_susc = S_traj_var_susc;
    results_var_susc.I_traj_var_susc = I_traj_var_susc;
    results_var_susc.R_traj_var_susc = R_traj_var_susc;

    marginal_eps_S_traj_var_susc = reshape(sum(S_traj_eps_delta_array,2),length(t_span),m)./S_traj;

    results_var_susc.marginal_eps_S_traj_var_susc = marginal_eps_S_traj_var_susc;

    for kk=1:length(t_span)
        % mean epsilon_S(t)
        mean_epsilon_S_traj_var_susc(kk) = sum(eps.*marginal_eps_S_traj_var_susc(kk,:));


    end

    results_var_susc.mean_epsilon_S_traj_var_susc = mean_epsilon_S_traj_var_susc;

    %% get Rt - variation in eps - only!!
    SIR_traj_var_susc = zeros(length(t_span),4);
    SIR_traj_var_susc(:,1) = S_traj_var_susc;
    SIR_traj_var_susc(:,2) = I_traj_var_susc;
    SIR_traj_var_susc(:,3) = R_traj_var_susc;
    SIR_traj_var_susc(:,4) = mean_epsilon_S_traj_var_susc;

    Rt_traj_var_susc = get_Rt_SIR_eps(params,SIR_traj_var_susc);
    results_var_susc.Rt_traj_var_susc = Rt_traj_var_susc;

end

%% Plotting
f2=figure(2); set(f2, 'Position', [900   50   400   930]);

subplot(3,1,1);

if run_classic_SIR

    plot(params.t_span, S_traj_SIR_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;
    plot(params.t_span, I_traj_SIR_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;
    plot(params.t_span, R_traj_SIR_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;

end

q(1)=plot(params.t_span, S_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
q(2)=plot(params.t_span, I_traj,'Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
q(3)=plot(params.t_span, R_traj,'Color',my_rgb_colors(3,:),'LineWidth',2); hold on;



% axis([0 t_end 10^-5 0.1*10]);
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
title('Dynamics')
legend(q,{'S','I','R'},'Location','NorthEast');
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

subplot(3,1,2);
r(1)=plot(params.t_span, mean_epsilon_S_traj,'-','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
r(2)=plot(params.t_span, mean_delta_I_traj,'-','Color',my_rgb_colors(3,:),'LineWidth',2); hold on;

% axis([0 t_end 0 1.1]);
% ylim([0 2])
title('Mean Trajectories');
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

legend(r,{'$\bar{\varepsilon}$','$\bar{\delta}_I$'},'Interpreter','Latex','Location','NorthEast');

subplot(3,1,3);

if run_classic_SIR

    plot(params.t_span, Rt_traj_classic,'Color',[0.65 0.65 0.65],'LineWidth',2); hold on;

end

plot(params.t_span, Rt_traj_eps_delta,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
% q(2) = plot(params.t_span, Rt_SIR_iv,'--','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
% q(3) = plot(params.t_span, Rt_SIR_v,':','Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
% axis([0 t_end 0 4]);
% ylim([0 4]);

xlabel('Time (days)'); ylabel([{'Effective'; 'Reproduction'; 'Number'}]);
% legend(q,{'Correlated','independent','classic'},'Location','SouthEast');
%'correlated',
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');


%% plot joint in S and I at certain time point
if want_to_plt_distributions

    this_ind = index_day_distribution;

    this_joint_S(:,:) = reshape(S_traj_eps_delta_array(this_ind,:,:),m,n)/S_traj(this_ind);
    this_joint_I(:,:) = reshape(I_traj_eps_delta_array(this_ind,:,:),m,n)/I_traj(this_ind);
    this_joint_R(:,:) = reshape(R_traj_eps_delta_array(this_ind,:,:),m,n)/R_traj(this_ind);

    this_marg_eps_S = dx*sum(this_joint_S);
    this_marg_delta_S = dx*sum(this_joint_S,2)';

    this_marg_eps_I = dx*sum(this_joint_I);
    this_marg_delta_I = dx*sum(this_joint_I,2)';

    plt_distributions(eps_plt, del, this_joint_S, this_joint_I, this_marg_eps_S, this_marg_delta_S, this_marg_eps_I, this_marg_delta_I, my_rgb_colors)

    % want to save this distribution?
    if save_distributions

        data.init_joint_S = this_joint_S;
        data.init_joint_I = this_joint_I;
        data.init_joint_R = this_joint_R;
        data.index_day_distribution = index_day_distribution;

        folder_location = './sim_results/';

        save(strcat(folder_location,filename_distributions_save),'data');

        fprintf('Saved Distribution to File: \n');
        fprintf(strcat(filename_distributions_save,'\n\n'));

    else

        fprintf('Distribution Not Saved. \n\n');

    end




end




%% collect results
% SIR - eps & delta
% trajectories

results.S_traj = S_traj;
results.S_traj_eps_delta_array = S_traj_eps_delta_array;
% results.S_traj_v = S_traj_v;

results.I_traj = I_traj;
results.I_traj_eps_delta_array = I_traj_eps_delta_array;

results.R_traj = R_traj;
results.R_traj_eps_delta_array = R_traj_eps_delta_array;

% effective reproduction number
results.Rt_traj_eps_delta = Rt_traj_eps_delta;

% final outbreak size
results.FOS_traj_eps_delta = FOS_traj_eps_delta;

% marginals over time
results.marginal_eps_S_traj = marginal_eps_S_traj; % susceptibility in S
results.marginal_delta_S_traj = marginal_delta_S_traj;
results.marginal_eps_I_traj = marginal_eps_I_traj;
results.marginal_delta_I_traj = marginal_delta_I_traj; % transmissibility in I



% mean value trajectories
results.mean_epsilon_S_traj = mean_epsilon_S_traj;
results.mean_delta_I_traj = mean_delta_I_traj;

% results.mu_delta_S = mu_delta_S;
% results.mu_eps_I = mu_eps_I;
% results.mu_delta_I = mu_delta_I;

% variance
results.var_eps_S_traj = variance_eps_S_traj;
results.var_delta_I_traj = variance_delta_I_traj;
% results.dispSe_traj = dispSe_traj;
% results.dispId_traj = dispId_traj;

results.joint_s = init_joint_S;
results.joint_i = init_joint_I;

results.marg_eps_s = init_marg_eps_S;
results.marg_delta_S=init_marg_delta_S;
results.marg_eps_I=init_marg_eps_I;
results.marg_delta_I=init_marg_delta_I;




%%
% save simulated results
if save_results==1

    folder_location = './sim_results/';

    if save_results_classic

        save(strcat(folder_location,filename_results),'params','results','results_classic','results_var_susc');
    else
        save(strcat(folder_location,filename_results),'params','results');
    end
    fprintf('Saved Results to File: \n');
    fprintf(strcat(filename_results,'\n'));

else

    fprintf('Results Not Saved. \n');

end

