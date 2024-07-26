function [S_traj_v, I_traj_v, R_traj_v, Rt_SIR_v, FOS_SIR_traj, params] = simulateClassic()
    default_colors = [0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250];  

    gam=1/10;
    bet=2*gam; % careful, beta is a built-in function
    delta_S = 1; % starting mean transmissibility of susceptible population
    delta_I = 1; % starting mean transmissibility of infectious population
    epsilon_I =1;
    epsilon_S =1; % starting mean susceptibility of susceptible population
    N = 1; %starting population; 
    
    n = 300 ; % n different transmissibility classes
    m = n; % m different susceptibilty classes
 
    end_val = 6;
    d = delta_S*ones(1,n); % homogeneity
    dd = end_val/n; %like dx
    e = epsilon_S*ones(1,n); % homogeneity
    de = end_val/m;
    D = d.*ones(n,n); % same across columns
    E = e'.*ones(m,m); % same across rows

    
    params.mean = 1; 
    params.beta = bet;
    params.gamma = gam;
    params.delta_S = delta_S;
    params.delta_I = delta_I;
    params.epsilon_S = epsilon_S;
    params.epsilon_I = epsilon_I;
    params.N = N;
    params.D = D;
    params.E = E;
    params.e = e;
    params.n = n;
    params.m = m;

    t_start = 0; t_end = 400; % 200 days is about 6-7 months; 250 is 8-9 months
    t_span = t_start:0.2:t_end;
    params.t_span = t_span;


    %% Initialize Distributions
    
    sigma = sqrt(1); %standard deviation

    marg_d_s = normpdf(d,delta_S,sigma); %still homogeneous because e is one value
    marg_d_s = marg_d_s./sum(marg_d_s); %normalize

    marg_e_s = normpdf(e,epsilon_S,sigma);
    marg_e_s = marg_e_s./sum(marg_e_s); %normalize

    marg_d_i = normpdf(d,delta_I,sigma);
    marg_d_i = marg_d_i./sum(marg_d_i); %normalize

    marg_e_i = normpdf(e,epsilon_I,sigma);
    marg_e_i = marg_e_i/sum(marg_e_i); % normalize

    %Joint Distributions 
    joint_s = marg_e_s'*marg_d_s; 
    joint_i = marg_e_i'*marg_d_i;

    %Calculated Means
    mu_delta_S = sum(d.*marg_d_s)/sum(marg_d_s);
    mu_eps_S = sum(e.*marg_e_s)/sum(marg_e_s);
    mu_delta_I = sum(d.*marg_d_i)/sum(marg_d_i);
    mu_eps_I = sum(e.*marg_e_i)/sum(marg_e_i);

    params.mu_eps_S = mu_eps_S;
    params.mu_delta_I = mu_delta_I;
    params.joint_s = joint_s;
    params.joint_i = joint_i;


    %% plotting initial distributions
    
    %plotInitDists(marg_e_s, marg_d_s, marg_e_i, marg_d_i, e,d, joint_s, joint_i,E,D, default_colors);

    %% Initialize Eigendirections

    eps = 1.1e-6; %change as needed
    
    params.eps = eps;
    params.end_val = end_val;

    [eigen_direction_SIR_ed] = get_eigendirection_SIRdelta_eps(params);


    S_init = 1 + eps*eigen_direction_SIR_ed(1);
    I_init = eps*eigen_direction_SIR_ed(2);
    R_init = eps*eigen_direction_SIR_ed(3);

    init_dist_S = joint_s*(S_init); %joint probability matrix S_{i,j}
    params.init_dist_S = init_dist_S;

    init_dist_I = joint_i*I_init; 
    params.init_dist_I = init_dist_I;

    init_dist_R = joint_s*R_init;



    init_conds = [reshape(init_dist_S, (m)*(n),1); reshape(init_dist_I,(m)*(n),1); reshape(init_dist_R, (m)*(n),1)];



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

    delta_I_dist_traj = sum(d.*marg_I_d_traj,2)./I_traj; %fishy
    mu_delta_I_traj = sum(delta_I_dist_traj,2);

    epsilon_I_dist_traj = sum(e.*marg_I_e_traj,2)./I_traj; %fishy
    mu_epsilon_I_traj = sum(epsilon_I_dist_traj,2);

    SIRed_traj = zeros(length(t_span),5);
    SIRed_traj(:,1) = S_traj;
    SIRed_traj(:,2) = I_traj;
    SIRed_traj(:,3) = R_traj;
    SIRed_traj(:,4) = mu_epsilon_S_traj;
    SIRed_traj(:,5) = mu_delta_I_traj;

    Rt_SIR_v = get_Rt_SIR_delta_eps(params,SIRed_traj);
    FOS_SIR_traj = getFinalOutbreakSize(I_traj, R_traj);


    %% Joint Distribution Trajectories for Sus. & Infected

    joint_s_traj = zeros(length(t_span),m,n);
    joint_i_traj = zeros(length(t_span),m,n);
    for t = 1:length(t_span)
        joint_s_traj(t,:,:) = marg_S_e_traj(t,:,:)'*marg_S_d_traj(t,:,:); 
        joint_i_traj(t,:,:) = marg_I_e_traj(t,:,:)'*marg_I_d_traj(t,:,:);
    end


    S_traj_v = S_traj;
    I_traj_v = I_traj;
    R_traj_v = R_traj;

end