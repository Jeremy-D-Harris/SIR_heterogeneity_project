% function void = main_SIRvariation_Gamma(void)

% simulate SIR model with transmissibility & susceptibility variation
% using Gamma distribution

%%
clear all; close all; clc;


%% want to save results?
save_results = 0;
% 0: don't save
% 1: save

filename_results = 'Gamma_update072824.mat';

%% options
% run & save Classic SIR
run_classic_SIR = 1;
save_results_classic = 1;

% run & save variation in susceptibility SIR
run_variation_susc_SIR = 1;
save_results_variation_susc = 1;


% save distributions during exponential growth
% want to plot distributions at certain time?
want_to_plt_distributions = 1;
save_distributions = 0; % save distribution at certain time?
index_day_distribution = 50; % what time? (days)

% want to read in distribution from a file?
readin_init_joint = 0;
% filename_distributions_load = 'Gamma_joint_expgrowth.mat';

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
t_start = 0; t_end = 300;
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


    %% Initialize Distributions

    % % shape parameters
    shape_eps = 10; % larger shape, smaller variance
    shape_delta = 3;

    params.shape_delta = shape_delta;
    params.shape_eps = shape_eps;

    % could better name 'intended'
    %
    set_mean_eps_S = 1;
    set_mean_delta_S = 1;
    set_mean_eps_I = 1;
    set_mean_delta_I = 1;
    %
    params.mean_eps_S = set_mean_eps_S;
    params.mean_delta_I = set_mean_delta_I;


    % susceptibility distribution in susceptible pop - g_S(eps)
    marg_eps_S = gampdf(eps,shape_eps,set_mean_eps_S/shape_eps);
    marg_eps_S = marg_eps_S./sum(marg_eps_S)/dx; %normalize

    %potential transmissibility distribution in susceptible pop - h_S(delta)
    marg_delta_S = gampdf(del,shape_delta, set_mean_delta_S/shape_delta);
    marg_delta_S = marg_delta_S./sum(marg_delta_S)/dx; %normalize

    %inactive susceptibility distribution in infectious pop
    marg_eps_I = gampdf(eps,shape_eps,set_mean_eps_I/shape_eps);
    marg_eps_I = marg_eps_I/sum(marg_eps_I)/dx; % normalize

    %effective transmissibility distribution in infectious pop - h_I(delta)
    marg_delta_I = gampdf(del,shape_delta, set_mean_delta_I/shape_delta);
    marg_delta_I = marg_delta_I./sum(marg_delta_I)/dx; %normalize


    %Independent Joint Distributions = Product of Marginals
    init_joint_S = marg_eps_S'*marg_delta_S;
    init_joint_I = marg_eps_I'*marg_delta_I;

    %Calculated Means
    mean_delta_S = sum(del.*marg_delta_S);
    mean_eps_S = sum(eps.*marg_eps_S);
    mean_delta_I = sum(del.*marg_delta_I);

    params.mean_eps_S = mean_eps_S;
    params.mean_delta_I = mean_delta_I;
    params.init_joint_S = init_joint_S;
    params.init_joint_I = init_joint_I;


end



init_marginal_eps_S = dx*sum(init_joint_S);
init_marginal_delta_S = dx*sum(init_joint_S,2)';

init_marginal_eps_I = dx*sum(init_joint_I);
init_marginal_delta_I = dx*sum(init_joint_I,2)';


% %% plotting initial distributions
% plotInitDists(marg_e_s, marg_d_s, marg_e_i, marg_d_i, eps,del, joint_s, joint_i,E,D, default_colors);

%% Initializing Eigendirections
eps_perturb = 1.25e-6; %change as needed to adjust onset

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
    mean_epsilon_S_traj(kk) = dx*sum(eps.*marginal_eps_S_traj(kk,:));

    % variance in susceptibility
    variance_eps_S_traj(kk) = sum((eps- mean_epsilon_S_traj(kk)*ones(size(eps))).^2.*marginal_eps_S_traj(kk,:));

    % I(t,eps,delta)
    this_I_eps_delta_array(:,:) = I_traj_eps_delta_array(kk,:,:);

    % joint in I
    this_joint_I = this_I_eps_delta_array/I_traj(kk)/dx/dx;
    joint_I_traj(kk,:,:) = this_joint_I;

    % marginals in I
    marginal_eps_I_traj(kk,:) = dx*reshape(sum(this_joint_I,1),1,m);
    marginal_delta_I_traj(kk,:) = dx*reshape(sum(this_joint_I,2),1,n);

    % mean transmissibility
    mean_delta_I_traj(kk) = dx*sum(del.*marginal_delta_I_traj(kk,:));

    % variance transmissibility
    variance_delta_I_traj(kk) = sum(((del- mean_delta_I_traj(kk)*ones(size(del))).^2).*marginal_delta_I_traj(kk,:));

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
SIR_traj_eps_delta(:,4) = mean_epsilon_S_traj;
SIR_traj_eps_delta(:,5) = mean_delta_I_traj;

Rt_traj_eps_delta = get_Rt_SIR_eps_delta(params,SIR_traj_eps_delta);


CV2_eps_traj = variance_eps_S_traj./(mean_epsilon_S_traj.^2);
CV2_delta_traj = variance_delta_I_traj./(mean_delta_I_traj.^2);



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

    init_S_eps_values = sum(init_values_S_eps_delta,2)';

    init_I_eps_values = sum(init_values_I_eps_delta,2)';

    init_R_eps_values = sum(init_values_R_eps_delta,2)';

    S_traj_eps_values =zeros(length(t_span),(m+1));
    I_traj_eps_values = zeros(length(t_span),(m+1));
    R_traj_eps_values = zeros(length(t_span),(m+1));

    % to pass to ODE function
    init_conds_var_susc = [init_S_eps_values; init_I_eps_values; init_R_eps_values];

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


    for kk=1:length(t_span)

        marginal_eps_S_traj_var_susc(kk,:) = S_traj_var_susc(kk,:)/S_traj(kk)/dx;
        mean_epsilon_S_traj_var_susc(kk) = dx*sum(eps.*marginal_eps_S_traj_var_susc(kk,:));
        variance_eps_S_traj_var_susc(kk) = sum((eps- mean_epsilon_S_traj_var_susc(kk)*ones(size(eps))).^2.*marginal_eps_S_traj_var_susc(kk,:));

    end

    results_var_susc.marginal_eps_S_traj_var_susc = marginal_eps_S_traj_var_susc;
    results_var_susc.mean_epsilon_S_traj_var_susc = mean_epsilon_S_traj_var_susc;
    results_var_susc.variance_eps_S_traj_var_susc = variance_eps_S_traj_var_susc;

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

q(1)=plot(params.t_span, S_traj,'Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
q(2)=plot(params.t_span, I_traj,'Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
q(3)=plot(params.t_span, R_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;



axis([0 t_end 0 1.1]);
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
title('Dynamics')
legend(q,{'S','I','R'},'Location','SouthWest');
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

subplot(3,1,2);
r(2)=plot(params.t_span, mean_delta_I_traj,'-','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
r(1)=plot(params.t_span, mean_epsilon_S_traj,'-','Color',my_rgb_colors(3,:),'LineWidth',2); hold on;


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

axis([0 t_end 0 2.1]);
% ylim([0 4]);

xlabel('Time (days)'); ylabel([{'Effective'; 'Reproduction'; 'Number'}]);
% legend(q,{'Correlated','independent','classic'},'Location','SouthEast');
%'correlated',
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

        save(strcat(folder_location,filename_distributions_save),'data');

        fprintf('Saved Distribution to File: \n');
        fprintf(strcat(filename_distributions_save,'\n\n'));

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

% variance trajectories
results.variance_eps_S_traj = variance_eps_S_traj;
results.variance_delta_I_traj = variance_delta_I_traj;



%%
% save simulated results
if save_results==1

    folder_location = './sim_results/';

    if (save_results_classic*save_results_variation_susc)

        % + classic + variation susceptibility
        save(strcat(folder_location,filename_results),'params','results','results_classic','results_var_susc');

    elseif save_results_classic

        % + classic
        save(strcat(folder_location,filename_results),'params','results','results_classic');

    elseif save_results_variation_susc

        % + variation susceptibility
        save(strcat(folder_location,filename_results),'params','results','results_var_susc');

    else

        % save variation in eps & delta, exclusively
        save(strcat(folder_location,filename_results),'params','results');

    end

    fprintf('Saved Results to File: \n');
    fprintf(strcat(filename_results,'\n'));

else

    fprintf('Results Not Saved. \n');

end




