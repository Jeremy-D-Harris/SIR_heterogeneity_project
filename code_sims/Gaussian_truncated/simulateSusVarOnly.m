function [S_traj, I_traj, R_traj, Rt_SIRed, FOSed_traj, params] = simulateSusVarOnly()

default_colors = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];  

%% parameters

gam=1/10;
bet=2*gam;
delta_S = 1.0; % starting mean transmissibility of susceptible population
delta_I =1; % starting mean transmissibility of infectious population
epsilon_S =1; % starting mean susceptibility of susceptible population
N = 1; %starting population; 

n = 300 ; % n different transmissibility classes
m = n; % m different susceptibilty classes


end_val = 6; 
d =linspace(0, end_val,n);
dd = end_val/n; %like dx
e = linspace(0,end_val,m);
de = end_val/m;
d = d + dd/2;
e = e + de/2;
D = d.*ones(n,n); % same across columns
E = e'.*ones(m,m); % same across rows

%NOTE: idx needs to be manually changed if end_val or n are changed
%idx = index where d is equal or closest to 1
idx = 50;


params.beta = bet;
params.gamma = gam;
params.delta_S = delta_S;
params.epsilon_S = epsilon_S;
params.N = N;
params.D = D;
params.E = E;

t_start = 0; t_end = 400; % 200 days is about 6-7 months; 250 is 8-9 months
t_span = t_start:0.2:t_end;
params.t_span = t_span;


%% Initialize Distributions
%potential transmissibility distribution in susceptible pop - g
marg_d_s = zeros(1,n);
marg_d_s(idx) = 1; 
marg_d_s = marg_d_s./sum(marg_d_s); %normalize

% susceptibility distribution in susceptible pop - h
marg_e_s =  normpdf(e, epsilon_S,sqrt(0.16));
marg_e_s = marg_e_s./sum(marg_e_s); %normalize

%effective transmissibility distribution in infectious pop - q
marg_d_i = zeros(1,n);
marg_d_i(idx) = 1; 
marg_d_i = marg_d_i./sum(marg_d_i); %normalize

%inactive susceptibility distribution in infectious pop - p
marg_e_i = normpdf(e,epsilon_S,sqrt(0.16));
marg_e_i = marg_e_i/sum(marg_e_i); % normalize

%Joint Distributions
joint_s = marg_e_s'*marg_d_s; 
joint_i = marg_e_i'*marg_d_i;

%Calculated Means
mu_delta_S = sum(d.*marg_d_s); 
mu_eps_S = sum(e.*marg_e_s);
mu_delta_I = sum(d.*marg_d_i);

params.mu_eps_S = mu_eps_S;
params.mu_delta_I = mu_delta_I;
params.joint_s = joint_s;
params.joint_i = joint_i;
params.mean = mu_eps_S;

%% plotting initial distributions
plotInitDists(marg_e_s, marg_d_s, marg_e_i, marg_d_i, e,d, joint_s, joint_i,E,D, default_colors);

%% Initializing Eigendirections
eps = 1.35e-6; %change as needed

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

S_traj = sum(S_traj_discrete,[2,3]);
I_traj = sum(I_traj_discrete,[2,3]);
R_traj = sum(R_traj_discrete,[2,3]);

marg_S_e_traj = reshape(sum(S_traj_discrete,3),length(t_span),m); 
marg_S_d_traj = reshape(sum(S_traj_discrete,2),length(t_span),n);
marg_I_e_traj = reshape(sum(I_traj_discrete,3),length(t_span),m);
marg_I_d_traj = reshape(sum(I_traj_discrete,2),length(t_span),n);

delta_S_dist_traj = sum(d.*marg_S_d_traj,2)./S_traj;
mu_delta_S_traj = sum(delta_S_dist_traj,2);

epsilon_S_dist_traj = sum(e.*marg_S_e_traj,2)./S_traj;
mu_epsilon_S_traj = sum(epsilon_S_dist_traj,2);

delta_I_dist_traj = sum(d.*marg_I_d_traj,2)./I_traj;
mu_delta_I_traj = sum(delta_I_dist_traj,2);

epsilon_I_dist_traj = sum(e.*marg_I_e_traj,2)./I_traj;
mu_epsilon_I_traj = sum(epsilon_I_dist_traj,2);

SIRed_traj = zeros(length(t_span),5);
SIRed_traj(:,1) = S_traj;
SIRed_traj(:,2) = I_traj;
SIRed_traj(:,3) = R_traj;
SIRed_traj(:,4) = mu_epsilon_S_traj;
SIRed_traj(:,5) = mu_delta_I_traj;

Rt_SIRed = get_Rt_SIR_delta_eps(params,SIRed_traj);

FOSed_traj = getFinalOutbreakSize(I_traj, R_traj);

%% Joint Distribution Trajectories for Sus. & Infected

joint_s_traj = zeros(length(t_span),m,n);
joint_i_traj = zeros(length(t_span),m,n);
for t = 1:length(t_span)
    joint_s_traj(t,:,:) = marg_S_e_traj(t,:,:)'*marg_S_d_traj(t,:,:); 
    joint_i_traj(t,:,:) = marg_I_e_traj(t,:,:)'*marg_I_d_traj(t,:,:);
end

%% Simulate classic SIR model

[S_traj_v, I_traj_v, R_traj_v, Rt_SIR_v, FOS_SIR_traj, params_v] = simulateClassic();


%% Calculate Variance
for kk=1:length(params.t_span)
    this_S_marg(1,:) = sum(S_traj_discrete(kk,:,:),3);
    marg_S_e_traj(kk,:) = (1/S_traj(kk))*this_S_marg;
end

var = zeros(length(params.t_span),1);
for tt = 1:length(params.t_span)
    var(tt) = sum((params.E(:,1)'- mu_epsilon_S_traj(tt)*ones(size(params.E(:,1)'))).^2.*marg_S_e_traj(tt,:));
    
end

%%
s = 20;
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
data.Svar = var(1:s:end);


%save('dataIVGaussianSusVarOnly.mat','params','data');


end