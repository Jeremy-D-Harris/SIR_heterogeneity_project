function [time_array,SIR_traj,dist_S_traj,dist_I_traj,dist_R_traj]= simulate_SIRdelta_eps_stochastic(params,init_dist_S,init_dist_I,init_dist_R)

%% parameters to local variables
N = params.N; % number of individuals
% num_steps = params.num_steps; % reaction per trajectory
% cycles = 100; % number of trajectories iterated over

bet=params.beta;
gam = params.gamma;
E = params.E;
D = params.D;

% m = length(E);
n = length(D);


%% read in init distribution
% init_dist_S = y(1:(m*n));
% init_dist_I = y((m*n+1):(2*m*n));
% init_dist_R = y((2*m*n+1):(3*m*n));

% initial populations
this_S = N-1;
this_I = 1;
this_R = 0;

% initial values: sum to N, population size
this_dist_S = init_dist_S;
this_dist_I = init_dist_I;
this_dist_R = init_dist_R;

% initial distributions f_S(eps,delta)*deps*ddelta = 1, 
% where f_S(eps,delta) = S(eps,delta)/S(t)

this_S_truedistribution = init_dist_S/this_S;
this_I_truedistribution = init_dist_I/this_I;
this_R_truedistribution = init_dist_R/this_R;



% %reshaping
% this_dist_S = reshape(this_dist_S,m,n);
% this_dist_I = reshape(this_dist_I,m,n);
% this_dist_R = reshape(this_dist_R,m,n);


this_t = 0;

time_array = [];
SIR_traj = [];
dist_traj = [];


%% main code loop: Gillespie alogrithm
% number of steps in time
% for k=1:num_steps

k = 1;
while this_I > 0


    % rates from differential equations
    transmission_rates_matrix = -bet*sum(sum(D.*this_dist_I)).*E.*this_S_truedistribution;
    recovery_rates_matrix = gam.*this_I_values;

    %reshaping
    transmission_rates_vector = reshape(transmission_rates_matrix,[m*n,1]);
    recovery_rates_vector = reshape(recovery_rates_matrix,[m*n,1]);

    %% create distribution from reaction rates
    % array of rates
    reaction_rates_array = [transmission_rates_vector,recovery_rates_vector];
    % total reaction rate
    total_rate = sum(reaction_rates_array);

    % get the next time based on the total rate, assuming exponential
    % distributed events
    u1=rand; %uniform random number
    tau = -(1/total_rate)*log(u1); % time elapsed
    this_t = this_t + tau;

    % cumul_rates = cumsum(dydt);
    tmp_weights = dydt/total_rate;
    tmp_order=cumsum(tmp_weights);

    % gives the index of the bin, corresponding to the reaction event drawn
    u2 = rand; % uniform random number
    x=find(u2<=tmp_order,1); % which event happens

    % get m (col index) & n (row index) indexes, which gives S(m,n), I(m,n)
    if mod(x,n) == 0
        % at the end
        mx = x/n; % integer value
    else

        mx = mod(x,n);

    end

    if mod(x,m) == 0
        % at the bottom
        nx = x/m; % integer value

    else

        nx = mod(x,m);

    end


    partition_rates = m*n+1;
    % m*n+1 = 1st rate of dIdt
    % 2*m*n+1 = 1st rate of dRdt

    % strictly less than lands on the left the last rate
    % e.g., number < m*n+1 includes up to the last rate in dSdt

    %% update populations based on which event happened
    if x < partition_rates(1)

        % transmission happened

        % these are truly distribution, meaning 
        % sum(this_S_truedistribution)*deps*ddelta = 1
        this_S_truedistribution(mx,nx) = this_S_truedistribution(mx,nx)-1/N;
        this_I_truedistribution(mx,nx) = this_I_truedistribution(mx,nx)+1/N;
        % this_dist_R =

        this_S = this_S - 1;
        this_I = this_S + 1;
        
        % these sume to N, the total population
        this_dist_S = this_S*this_S_truedistribution;
        this_dist_I = this_I*this_I_truedistribution;
        this_dist_R = this_R*this_R_truedistribution;

    else
        % recovery happened
        % this_dist_S =
        this_I_truedistribution(mx,nx) = this_I_truedistribution(mx,nx)-1/N;
        this_R_truedistribution(mx,nx) = this_R_truedistribution(mx,nx)+1/N;

        this_I = this_I - 1;
        this_R = this_R + 1;

        this_dist_S = this_S*this_S_truedistribution;
        this_dist_I = this_I*this_I_truedistribution;
        this_dist_R = this_R*this_R_truedistribution;



    end
    % end

    %% store stuff

    time_array(k) =  this_t;
    SIR_traj(k,:) = [this_S,this_I,this_R];
    dist_S_traj(k,:,:) = this_dist_S;
    dist_I_traj(k,:,:) = this_dist_I;
    dist_R_traj(k,:,:) = this_dist_R;

    k = k+1;

end

end
