% simulate SIR model with transmissibility & susceptibility variation
% Version: Independent Gaussian (Not Truncated)

%% 

clear;close all; clc;

default_colors = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];  

%% parameters
gam=1/10;

bet=2*gam; %change bet accordingly if want to match exponential growth rates

delta_S = 1; % starting mean transmissibility of susceptible population
delta_I = 1; % starting mean transmissibility of infectious population
epsilon_S =1;% starting mean susceptibility of susceptible population
epsilon_I =1;
N = 1; %starting population; 

n = 200 ; % n different transmissibility classes
m = n; % m different susceptibilty classes

end_val = 3;
ifVar = 1; % 1 - True, 0 - False
d =getTransValues(n,ifVar,end_val);
dd = end_val/n; %like dx
e = linspace(0,end_val,m);
de = end_val/m;
d = d + dd/2;
e = e + de/2;
D = d.*ones(n,n); % same across columns
E = e'.*ones(m,m); % same across rows


params.beta = bet;
params.gamma = gam;
params.mean = 1; 
params.delta_S = delta_S;
params.epsilon_S = epsilon_S;
params.delta_I = delta_I;
params.epsilon_I = epsilon_I;
params.N = N;
params.D = D;
params.E = E;
params.end_val = end_val;

t_start = 0; t_end = 300; % 200 days is about 6-7 months; 250 is 8-9 months
t_span = t_start:1:t_end;
params.t_span = t_span;


%% Initialize Joint Distributions & calculate Marginal
% joint Guassian distribution the correlation matrix (Sigma)
[X1, X2] = meshgrid(e,d);
X = [X1(:) X2(:)];
mus = [delta_S epsilon_S]; %means for S
mui = [delta_I epsilon_I]; %means for I


Sigmas = [0.04 0.05; 0.05 0.16]; %change off-diagonals for correlation

Sigmai = [0.04 0.05; 0.05 0.16]; %change diagonals for variance

joint_s = mvnpdf(X, mus,Sigmas); 
joint_s = reshape(joint_s, length(d),length(e));

joint_i = mvnpdf(X, mui,Sigmai);
joint_i = reshape(joint_i, length(d), length(e));

joint_s = joint_s./sum(sum(joint_s),2); %normalize here 
joint_i = joint_i./sum(sum(joint_i),2);


marg_d_s = sum(joint_s); 
marg_e_s = sum(joint_s,2)'; 

marg_d_i = sum(joint_i);
marg_e_i = sum(joint_i,2)'; 

%Calculated Means
mu_delta_S = sum(d.*marg_d_s); 
mu_eps_S = sum(e.*marg_e_s);
mu_delta_I = sum(d.*marg_d_i);
mu_eps_I = sum(e.*marg_e_i);



params.mu_eps_S = mu_eps_S;
params.mu_delta_S = mu_delta_S;
params.mu_eps_I = mu_eps_I;
params.mu_delta_I = mu_delta_I;
params.joint_s = joint_s;
params.joint_i = joint_i;
params.Sigmas = Sigmas;
params.Sigmai = Sigmai;

%Computed Covariance, Variance, Standard Deviation and Correlation Coefficient
cov = sum(sum(((E(:,1)- mu_eps_S*ones(size(E(:,1))))*(D(1,:)- mu_delta_S*ones(size(D(1,:))))).*joint_s),2);
var_e = sum((E(:,1)'- mu_eps_S*ones(size(E(:,1)'))).^2.*marg_e_s);
sd_e = sqrt(var_e);
var_d = sum((D(1,:)- mu_delta_S*ones(size(D(1,:)))).^2.*marg_d_s,2);
sd_d = sqrt(var_d);
corrcoef = cov/(sd_e*sd_d); %'actual' correlation coefficient
params.corrcoef = corrcoef;

%% plotting initial distributions

plotInitDists(marg_e_s, marg_d_s, marg_e_i, marg_d_i, e,d, joint_s, joint_i,E,D, default_colors);

%% Initializing Eigendirections
eps = 1e-6; % change as needed

params.eps = eps;

[eigen_direction_SIR_ed] = get_eigendirection_SIRdelta_eps(params);

S_init = 1 + eps*eigen_direction_SIR_ed(1);
I_init = eps*eigen_direction_SIR_ed(2);
R_init = eps*eigen_direction_SIR_ed(3);

init_dist_S = joint_s*(S_init); %joint probability matrix S_{i,j}
params.init_dist_S = init_dist_S;

init_dist_I = joint_i*I_init; 
params.init_dist_I = init_dist_I;

init_dist_R = joint_s*R_init;

init_conds = [reshape(init_dist_S, m*n,1); reshape(init_dist_I,m*n,1); reshape(init_dist_R, m*n,1)];

%% Simulate model

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

S_traj_discrete =zeros(length(t_span),m,n);
I_traj_discrete = zeros(length(t_span),m,n);
R_traj_discrete = zeros(length(t_span),m,n);


[t,y_traj] = ode45(@(t,y)simulate_SIRdelta_eps(t,y,params), params.t_span, init_conds, options);

S_traj_discrete(:,:,:) = abs(reshape(y_traj(:,1:m*n), length(t_span),m,n));
I_traj_discrete(:,:,:) = abs(reshape(y_traj(:,m*n+1:2*m*n),length(t_span),m,n));
R_traj_discrete(:,:,:) = abs(reshape(y_traj(:,2*m*n+1:3*m*n),length(t_span),m,n));

%% S,I,R & Marginal eps & delta Distribution Trajectories

S_traj = sum(S_traj_discrete,[2,3]); % population size
I_traj = sum(I_traj_discrete,[2,3]);
R_traj = sum(R_traj_discrete,[2,3]);

marg_S_e_traj = reshape(sum(S_traj_discrete,3),length(t_span),m); 
marg_S_d_traj = reshape(sum(S_traj_discrete,2),length(t_span),n);
marg_I_e_traj = reshape(sum(I_traj_discrete,3),length(t_span),m);
marg_I_d_traj = reshape(sum(I_traj_discrete,2),length(t_span),n);

delta_S_dist_traj = sum(d.*marg_S_d_traj,2); 
mu_delta_S_traj = sum(delta_S_dist_traj,2)./S_traj;  

epsilon_S_dist_traj = sum(e.*marg_S_e_traj,2);
mu_epsilon_S_traj = sum(epsilon_S_dist_traj,2)./S_traj;

delta_I_dist_traj = sum(d.*marg_I_d_traj,2);
mu_delta_I_traj = sum(delta_I_dist_traj,2)./I_traj;

epsilon_I_dist_traj = sum(e.*marg_I_e_traj,2);
mu_epsilon_I_traj = sum(epsilon_I_dist_traj,2)./I_traj;

SIRed_traj = zeros(length(t_span),5);
SIRed_traj(:,1) = S_traj;
SIRed_traj(:,2) = I_traj;
SIRed_traj(:,3) = R_traj;
SIRed_traj(:,4) = mu_epsilon_S_traj;
SIRed_traj(:,5) = mu_delta_I_traj;

Rt_SIRed = get_Rt_SIR_delta_eps(params,SIRed_traj);



%% Joint Distribution Trajectories for Sus. & Infected

joint_s_traj = zeros(length(t_span),m,n);
joint_i_traj = zeros(length(t_span),m,n);
for t = 1:length(t_span)
    joint_s_traj(t,:,:) = marg_S_e_traj(t,:,:)'*marg_S_d_traj(t,:,:); 
    joint_i_traj(t,:,:) = marg_I_e_traj(t,:,:)'*marg_I_d_traj(t,:,:);
end

%% Calculate final size of outbreak

FOSed_traj = getFinalOutbreakSize(I_traj, R_traj);

%% Simulate classic SIR model

[S_traj_v, I_traj_v, R_traj_v, Rt_SIR_v, FOS_SIR_traj, params_v] = simulateClassic();


%% Plotting
f6=figure(6); set(f6, 'Position', [900   50   400   620]);
set(gca,'LineWidth',1,'FontSize',14, 'FontWeight','normal','FontName','Times');

subplot(3,1,1);
q(1)=semilogy(params.t_span, S_traj,'Color',default_colors(1,:),'LineWidth',2); hold on;
q(2)=semilogy(params.t_span, I_traj,'Color',default_colors(2,:),'LineWidth',2); hold on;
q(3)=semilogy(params.t_span, R_traj,'Color',default_colors(3,:),'LineWidth',2); hold on;
q(4)=semilogy(params.t_span, S_traj_v,':','Color',default_colors(1,:),'LineWidth',2); hold on;
q(5)=semilogy(params.t_span, I_traj_v,':','Color',default_colors(2,:),'LineWidth',2); hold on;
q(6)=semilogy(params.t_span, R_traj_v,':','Color',default_colors(3,:),'LineWidth',2); hold on;

axis([0 t_end 10^-5 0.1*10]);
xlabel('Time (days)'); ylabel({'Population'; 'Fraction'});

legend(q,{'S','I','R','S-classic','I-classic','R-classic'},'Location','SouthEast');


subplot(3,1,2);
q(1)=plot(params.t_span, mu_epsilon_S_traj,'-','Color',default_colors(2,:),'LineWidth',2); hold on;
q(2)=plot(params.t_span, mu_delta_I_traj,'-','Color',default_colors(3,:),'LineWidth',2); hold on;

axis([0 t_end 0 3]);
ylim([0 2])
xlabel('Time (days)'); ylabel({'Mean'; 'Trajectories'});


legend(q,{'$\mu_{\epsilon,S}$','$\mu_{\delta,I}$'},'Interpreter','Latex','Location','SouthEast');

subplot(3,1,3)
q(1) = plot(params.t_span, Rt_SIRed,'Color',default_colors(1,:),'LineWidth',2); hold on;
q(2) = plot(params.t_span, Rt_SIR_v,':','Color',default_colors(3,:),'LineWidth',2); hold on;
axis([0 t_end 0 5]);
ylim([0 4]);
xlabel('Time (days)'); ylabel([{'Effective'; 'Reproduction'; 'Number'}]);
legend(q,{'Independent','Classic'},'Location','SouthEast');


f7 = figure(7);
set(gca,'FontSize', 14, 'LineWidth',1,'FontWeight','normal','FontName','Times');
q(1)=plot(params.t_span, 100.*FOSed_traj,'Color',default_colors(1,:),'LineWidth',2); hold on;
q(2)=plot(params.t_span, 100.*FOS_SIR_traj,':','Color',default_colors(3,:),'LineWidth',2); hold on;
axis([0 t_end 10^(-6) 10]);
ylim([0,105]);
xlabel('Time (days)'); ylabel({'Cumulated Infected (%)'});
legend(q,{'Independent', 'Classic'},'Location', 'SouthEast');
title('Cumulative Infected (I + R)');


%% Save Variables
s = 1;
params.s = s;
params.t_span = t_span(1:s:end);
data.FOS_SIR_traj = FOS_SIR_traj(1:s:end);
data.FOSed_traj = FOSed_traj(1:s:end);
data.I_traj = I_traj(1:s:end);
data.I_traj_discrete = I_traj_discrete(1:s:end,:,:);
data.I_traj_v = I_traj_v(1:s:end);
data.marg_I_e_traj = marg_I_e_traj(1:s:end,:);
data.marg_I_d_traj = marg_I_d_traj(1:s:end,:);
data.marg_S_d_traj = marg_S_d_traj(1:s:end,:);
data.marg_S_e_traj = marg_S_e_traj(1:s:end,:);
data.mu_delta_I_traj = mu_delta_I_traj(1:s:end);
data.mu_delta_S_traj = mu_delta_S_traj(1:s:end);
data.mu_epsilon_I_traj = mu_epsilon_I_traj(1:s:end);
data.mu_epsilon_S_traj = mu_epsilon_S_traj(1:s:end);
data.R_traj = R_traj(1:s:end);
data.R_traj_discrete = R_traj_discrete(1:s:end,:,:);
data.R_traj_v = R_traj_v(1:s:end);
data.Rt_SIR_v = Rt_SIR_v(1:s:end);
data.Rt_SIRed = Rt_SIRed(1:s:end);
data.S_traj = S_traj(1:s:end);
data.S_traj_discrete = S_traj_discrete(1:s:end,:,:);
data.S_traj_v = S_traj_v(1:s:end);

%save('dataFigure2LowVariance.mat','params','data','params_iv', 'params_v');
