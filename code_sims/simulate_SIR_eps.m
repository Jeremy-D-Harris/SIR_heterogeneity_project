function dydt = simulate_SIR_eps(t,y,params)

%% for SIR variation in susceptibility alone:
% susceptibility (eps) 

% note! 'S_eps_array' is an array of values that sum to 1
% If you divide this by S, you'll get the distribution, f_S(eps)

% parameters to local variables
bet = params.bet;
gam = params.gam;
eps_values = params.eps;
N = params.N;
% D = params.D;

m = length(eps_values);
% n = length(D);

% S(t,eps),I(t,eps),R(t,eps)
S_eps_values(1,:) = y(1:m); 
I_eps_values(1,:) = y((m+1):(2*m)); 
R_eps_values(1,:) = y((2*m+1):(3*m));

% force of infection
incident_infections = bet*sum(I_eps_values)*eps_values.*S_eps_values/N;

% recovery
recovery = gam.*I_eps_values;

%differential equations
dSdt = -incident_infections;

dIdt = incident_infections - recovery; 

dRdt = recovery;

% %reshaping
dSdt = transpose(dSdt);
dIdt = transpose(dIdt);
dRdt = transpose(dRdt);


dydt = [dSdt; dIdt; dRdt];
