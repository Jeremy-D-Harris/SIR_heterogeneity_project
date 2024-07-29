function dydt = simulate_SIR_eps_delta(t,y,params)

%% for SIR variation in both:
% susceptibility (eps) & transmissibility (delta)

% note! 'S_eps_delta_array' is an array of values that sum to 1
% If you divide this by S, you'll get the distribution, f_S(eps,delta)

% parameters to local variables
bet=params.bet; 
gam = params.gam;
E = params.E; % epsilon array
D = params.D; % delta array

m = length(E);
n = length(D);
dx = 1/m;

S_eps_delta_values = y(1:(m*n)); 
I_eps_delta_values = y((m*n+1):(2*m*n)); 
R_eps_delta_values = y((2*m*n+1):(3*m*n));

%reshapeing
S_eps_delta_array = reshape(S_eps_delta_values,m,n);
I_eps_delta_array = reshape(I_eps_delta_values,m,n);
R_eps_delta_array = reshape(R_eps_delta_values,m,n);


% incident infections
incident_infections = bet*sum(sum(D.*I_eps_delta_array)).*E.*S_eps_delta_array;

% recovery
recovery = gam.*I_eps_delta_array;


%differential equations
dSdt = -incident_infections;

dIdt = incident_infections - recovery; 

dRdt = recovery;

%reshaping
dSdt = reshape(dSdt,[m*n,1]);
dIdt = reshape(dIdt,[m*n,1]);
dRdt = reshape(dRdt,[m*n,1]);


dydt = [dSdt; dIdt; dRdt];

