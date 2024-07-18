% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation
%using negative binomial distributions

%%

clear all; close all; clc;

%% want to save?
save_ans = 1;
% 0: don't save
% 1: save

% filename = 'negbinomial_Poissonlike.mat';
% filename = 'negbinomial_overdispersed_uncorr_update07172024.mat';
filename = 'negbinomial_overdispersed_poscorr_update07172024.mat';

%  = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];
my_rgb_colors = [78 132 193; 209 109 106; 236 180 118]/255;

%% parameters
% recovery rate
gam=1/10; % assume mean infectious period is 10 days

%% variation in susceptibility/transmissibility
% using joint negative binomial from Felix Famoye:
% https://doi.org/10.1080/02664760902984618

% dispersion parameters: kappa->infty = Poisson; kappa=1 = Geometric;
% kappa<1 overdispersed
kappa1 = 50;  % dispersion (susceptibility)
kappa2 = 0.5; % dispersion (transmissibility)

% For the Singapore outbreak (Lloyd-smith, 2005), the maximum-likelihood estimate k̂ is 0.16
% (90% confidence interval 0.11–0.64), indicating highly overdispersed

%% correlation coefficient for initial joint in S population
rho_S = 0.6; % 0, 0.6

intended_R0 = 3;
% match basic reproduction numbers: R0=3
bet=intended_R0*gam; %Poisson-like, rho=0
% bet = intended_R0*gam*.5; % overdispersed, rho=0.6

% want means of marginals to be 1
mean_joint = 1;

% starting population
N = 1; % population size

% discrete mesh
n = 100; % number transmissibility classes
m = n; % susceptibilty classes

dx = 1;
eps = 0:dx:m;
del = 0:dx:n;

eps_plt = eps;
eps(1)=0.05;
% del(1)=0.05;

[Eps,Del] = meshgrid(eps,del);

% create parameter structure
params.beta = bet;
params.gamma = gam;
params.mean = mean_joint;
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


%% Initialize Joint Distributions & calculate Marginal
% can use other bivariation distributions, including nonparametric !
% x1 = eps, x2 = delta

m1 = 1/kappa1; % reciprocal of dispersion (susceptibility)
m2 = 1/kappa2; % reciprocal of dispersion (transmissibility)


%theta values - theta_t = mt*mean_joint/(mt*mean_joint +1)
% want means of 1
mu1 = 1;
mu2 = mu1;
theta1 = m1*mu1/(m1*mu2+1);
theta2 = m2*mu2/(m2*mu2+1);
% rearranging: means of marginals
% mu1 = theta1/(1-theta1)/m1;
% mu2 = theta2/(1-theta2)/m2;

% variances of marginals
var1 = theta1/(1-theta1)^2/m1;
var2 = theta2/(1-theta2)^2/m2;

c1 = ((1-theta1)/(1-theta1*exp(-1)))^(1/m1);
c2 = ((1-theta2)/(1-theta2*exp(-1)))^(1/m2);

A1 = m1^(-1)*theta1*exp(-1)/(1-theta1*exp(-1))-m1^(-1)*theta1/(1-theta1);
A2 = m2^(-1)*theta2*exp(-1)/(1-theta2*exp(-1))-m2^(-1)*theta2/(1-theta2);

% calculate lambda - multiplicative factor
lambda_S = rho_S*(sqrt(var1)*sqrt(var2))/(c1*c2*A1*A2);

%bivariate probability function for S compartment
init_joint_S = bivariate_negativebinomial_pdf(Eps,Del,theta1,theta2,m1,m2,lambda_S);
init_joint_I = init_joint_S;

%normalize here
init_joint_S = init_joint_S./sum(sum(init_joint_S));
init_joint_I = init_joint_I./sum(sum(init_joint_I));

init_marg_eps_S = sum(init_joint_S);
init_marg_delta_S = sum(init_joint_S,2)';

init_marg_eps_I = sum(init_joint_I);
init_marg_delta_I = sum(init_joint_I,2)';

% calculated mean values
init_mu_eps_S = sum(eps.*init_marg_eps_S);
init_mu_delta_S = sum(del.*init_marg_delta_S);
init_mu_eps_I = sum(eps.*init_marg_eps_I);
init_mu_delta_I = sum(del.*init_marg_delta_I);

%calculated values
init_covariance = sum(sum(((eps' - init_mu_eps_S*ones(size(eps')))*(del- init_mu_delta_S*ones(size(del)))).*init_joint_S));

init_variance_eps = sum((eps - init_mu_eps_S*ones(size(eps))).^2.*init_marg_eps_S);
init_sd_eps = sqrt(init_variance_eps);

init_variance_delta = sum((del- init_mu_delta_S*ones(size(del))).^2.*init_marg_delta_S);
init_sd_delta = sqrt(init_variance_delta);

corrcoef = init_covariance/(init_sd_eps*init_sd_delta);

%% plot initial distributions - if you want!
if 1
    
    
    f1=figure(1); set(f1, 'Position', [10   50   840   350]);
    subplot(1,2,1);
    
    imagesc(eps_plt,del,init_joint_S);
    % pcolor(eps,del,joint_S);
    %axis xy;
    set(gca,'YDir','normal');
    colorbar;
    xlim([0 4.5]); ylim([0 4.5]);
    
    xlabel('susceptibility $\varepsilon$','interpreter','latex');
    ylabel({'transmissibility $\delta$'},'interpreter','latex');
    zlabel('Population Fraction')
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Joint Distribution in S');
    
    % cbh = colorbar();
    % set color range
    % caxis([0,.2])
    % set ticks
    % set(cbh, 'YTick', [0.001, 0.01, 0.05, .1], ...
    %     'YTickLabel', {'p=0.001', 'p=0.01', 'p=0.05', 'p=.1'})
    % cbh.Label.String = 'Population Density';
    
    
    % f1=figure(1);
    subplot(1,2,2);
    q(1)=plot(eps_plt,init_marg_eps_S,'.-','Color',my_rgb_colors(3,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    q(2)=plot(del, init_marg_delta_S,'.-','Color',my_rgb_colors(2,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    % cat_margs = [marg_eps_S;marg_delta_S];
    % q=bar(eps,cat_margs);
    axis([0 4.5 0 1]);
    
    axis([0 4 0 1]);
    title('Marginal Distributions');
    xlabel('transmissibility d');
    
    xlabel('susceptibility $\varepsilon$ or transmissibility $\delta$','interpreter','latex');
    ylabel({'Population Fraction'});
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    legend(q,{'Susceptibility','Potential Transmissibility'},'Location','NorthEast');
    
end

%%
% pause;
% small perturbation in direction of eigenvector
eps_perturb = 5e-10;

params.eps_perturb = eps_perturb;

params.mu_eps_S = init_mu_eps_S;
params.mu_delta_I = init_mu_delta_I;

eigen_direction_SIR = get_eigendirection_SIRdelta_eps(params);

S_init = 1 - eps_perturb*abs(eigen_direction_SIR(1));
I_init = eps_perturb*abs(eigen_direction_SIR(2));
R_init = eps_perturb*abs(eigen_direction_SIR(3));

% pause;

%
init_dist_S = init_joint_S.*(S_init); %joint probability matrix S_{i,j}
params.init_dist_S = init_dist_S;


init_dist_I = init_joint_I.*I_init;
params.init_dist_I = init_dist_I;

init_dist_R = init_joint_S.*R_init;

%check should equal 1
%sum(sum((init_dist_S + init_dist_I + init_dist_R),2));

% need to reshape from matrix to vector from in order to pass to ODE function
init_conds = [reshape(init_dist_S, (m+1)*(n+1),1); reshape(init_dist_I,(m+1)*(n+1),1); reshape(init_dist_R, (m+1)*(n+1),1)];

%% Simulate model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

S_traj_discrete =zeros(length(t_span),(m+1),(n+1));
I_traj_discrete = zeros(length(t_span),(m+1),(n+1));
R_traj_discrete = zeros(length(t_span),(m+1),(n+1));

[t,y_traj] = ode45(@(t,y)simulate_SIRdelta_eps(t,y,params), params.t_span, init_conds, options);

% reshape from vector to matrix form
S_traj_discrete = reshape(y_traj(:,1:(m+1)*(n+1)), length(t_span),(m+1),(n+1));
I_traj_discrete = reshape(y_traj(:,((m+1)*(n+1)+1):(2*(m+1)*(n+1))),length(t_span),(m+1),(n+1));
R_traj_discrete = reshape(y_traj(:,(2*(m+1)*(n+1)+1):(3*(m+1)*(n+1))),length(t_span),(m+1),(n+1));


% S,I,R Trajectories
S_traj = sum(S_traj_discrete,[2,3]); % population sizes
I_traj = sum(I_traj_discrete,[2,3]);
R_traj = sum(R_traj_discrete,[2,3]);
FOS_traj = I_traj + R_traj;

marginal_eps_S_traj = reshape(sum(S_traj_discrete,2),length(t_span),m+1);
marginal_delta_S_traj = reshape(sum(S_traj_discrete,3),length(t_span),n+1);
marginal_eps_I_traj = reshape(sum(I_traj_discrete,2),length(t_span),m+1);
marginal_delta_I_traj = reshape(sum(I_traj_discrete,3),length(t_span),n+1);

for kk=1:length(t_span)
    % mean epsilon_S(t), delta_I(t)
    mu_epsilon_S_traj(kk) = sum(eps.*marginal_eps_S_traj(kk,:))./S_traj(kk);
    mu_delta_I_traj(kk) = sum(del.*marginal_delta_I_traj(kk,:))./I_traj(kk);
    
end


SIRed_traj = zeros(length(t_span),5);
SIRed_traj(:,1) = S_traj;
SIRed_traj(:,2) = I_traj;
SIRed_traj(:,3) = R_traj;
SIRed_traj(:,4) = mu_epsilon_S_traj;
SIRed_traj(:,5) = mu_delta_I_traj;

Rt_SIR_traj = get_Rt_SIR_delta_eps(params,SIRed_traj);


for kk = 1:length(t_span)
    
    %     marginal_eps_S_traj(kk,:) = marginal_eps_S_traj(kk,:)/sum(marginal_eps_S_traj(kk,:));
    %     marginal_delta_S_traj(kk,:) = marginal_delta_S_traj(kk,:)/sum(marginal_delta_S_traj(kk,:));
    %
    %     marginal_eps_I_traj(kk,:) = marginal_eps_I_traj(kk,:)/sum(marginal_eps_I_traj(kk,:));
    %     marginal_delta_I_traj(kk,:) = marginal_delta_I_traj(kk,:)/sum(marginal_delta_I_traj(kk,:));
    
    variance_eps_S_traj(kk) = sum((eps- mu_epsilon_S_traj(kk)*ones(size(eps))).^2.*marginal_eps_S_traj(kk,:));
    variance_delta_I_traj(kk) = sum(((del- mu_delta_I_traj(kk)*ones(size(del))).^2).*marginal_delta_I_traj(kk,:));
    
    sd_S_eps_traj(kk) = sqrt(variance_eps_S_traj(kk));
    sd_I_delta_traj(kk) = sqrt(variance_delta_I_traj(kk));
    
end

%dispersion parameters
% dispSe_traj = std_Se_traj./mu_epsilon_S_traj;
% dispId_traj = std_Id_traj./mu_delta_I_traj;

CV_Se_traj = variance_eps_S_traj./(mu_epsilon_S_traj.^2);
CV_Id_traj = variance_delta_I_traj./(mu_delta_I_traj.^2);


%% Plotting
f2=figure(2); set(f2, 'Position', [900   50   400   930]);

subplot(3,1,1);
q(1)=plot(params.t_span, S_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
q(2)=plot(params.t_span, I_traj,'Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
q(3)=plot(params.t_span, R_traj,'Color',my_rgb_colors(3,:),'LineWidth',2); hold on;


axis([0 t_end 10^-5 0.1*10]);
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
title('Dynamics')
legend(q,{'S','I','R'},'Location','NorthEast');
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

subplot(3,1,2);
r(1)=plot(params.t_span, mu_epsilon_S_traj,'-','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
r(2)=plot(params.t_span, mu_delta_I_traj,'-','Color',my_rgb_colors(3,:),'LineWidth',2); hold on;

% axis([0 t_end 0 1.1]);
% ylim([0 2])
title('Mean Trajectories');
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

legend(r,{'$\bar{\varepsilon}$','$\bar{\delta}_I$'},'Interpreter','Latex','Location','NorthEast');

subplot(3,1,3);
plot(params.t_span, Rt_SIR_traj,'Color',my_rgb_colors(1,:),'LineWidth',2); hold on;
% q(2) = plot(params.t_span, Rt_SIR_iv,'--','Color',my_rgb_colors(2,:),'LineWidth',2); hold on;
% q(3) = plot(params.t_span, Rt_SIR_v,':','Color',my_rgb_colors(3,:),'LineWidth',2); hold on;
% axis([0 t_end 0 4]);
% ylim([0 4]);

xlabel('Time (days)'); ylabel([{'Effective'; 'Reproduction'; 'Number'}]);
% legend(q,{'Correlated','independent','classic'},'Location','SouthEast');
%'correlated',
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');


%% plot joint in S and I at certain time point
if 0
    
    f3=figure(3); set(f3, 'Position', [100   450   840   350]);
    
    this_ind = 50;
    subplot(1,2,1);
    this_joint_S = reshape(S_traj_discrete(this_ind,:,:),(m+1),(n+1))/S_traj(this_ind);
    imagesc(eps_plt,del,this_joint_S);
    % pcolor(eps,del,joint_S);
    
    set(gca,'YDir','normal');
    colorbar;
    xlim([-0.5 6.5]); ylim([-0.5 6.5]);
    caxis([0 0.35]);
    
    xlabel('susceptibility $\varepsilon$','interpreter','latex');
    ylabel({'transmissibility $\delta$'},'interpreter','latex');
    
    f2=gca;
    f2.LineWidth = 1;
    f2.FontSize = 14;
    f2.FontWeight = 'normal';
    f2.FontName = 'Times';
    
    subplot(1,2,2);
    this_joint_I = reshape(I_traj_discrete(this_ind,:,:),(m+1),(n+1))/I_traj(this_ind);
    imagesc(eps_plt,del,this_joint_I);
    % pcolor(eps,del,joint_S);
    %axis xy;
    set(gca,'YDir','normal');
    colorbar;
    xlim([-0.5 6.5]); ylim([-0.5 6.5]);
    % clim([]);
    caxis([0 0.35]);
    
    xlabel('susceptibility $\varepsilon$','interpreter','latex');
    ylabel({'transmissibility $\delta$'},'interpreter','latex');
    
    f2=gca;
    f2.LineWidth = 1;
    f2.FontSize = 14;
    f2.FontWeight = 'normal';
    f2.FontName = 'Times';
    
end


%% Save results

results.S_traj = S_traj;
results.S_traj_discrete = S_traj_discrete;
% results.S_traj_v = S_traj_v;

results.I_traj = I_traj;
results.I_traj_discrete = I_traj_discrete;
% results.I_traj_v = I_traj_v;

results.R_traj = R_traj;
results.R_traj_discrete = R_traj_discrete;
% results.R_traj_v = R_traj_v(1:s:end);

results.Rt_SIR_traj = Rt_SIR_traj;
% results.Rt_SIR_v = Rt_SIR_v(1:s:end);

results.FOS_traj = FOS_traj;
% results.FOS_SIR_traj_v = FOS_SIR_traj_v;

% marginals over time
results.marg_eps_I_traj = marginal_eps_I_traj;
results.marg_delta_I_traj = marginal_delta_I_traj;
results.marg_delta_S_traj = marginal_delta_S_traj;
results.marg_eps_S_traj = marginal_eps_S_traj;

% mean value trajectories
results.mu_epsilon_S_traj = mu_epsilon_S_traj;
results.mu_delta_I_traj = mu_delta_I_traj;

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


%% homogeneous
% results.FOS_traj_v = FOS_SIR_traj_v;

%%
% save simulated results
if save_ans==1
    
    folder_location = './sim_results/';
    save(strcat(folder_location,filename),'params','results');
    
    fprintf('Saved to file: \n');
    fprintf(strcat(filename,'\n'));
    
else
    
    fprintf('Not Saved. \n');
    
end

