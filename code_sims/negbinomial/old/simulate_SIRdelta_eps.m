function dydt = simulate_SIRdelta_eps(t,y,params)

% note! 'S_eps_delta_array' is an array of values that sum to 1
% If you divide this by S, you'll get the distribution, f_S(eps,delta)

% parameters to local variables
bet=params.beta; 
gam = params.gamma;
E = params.E;
D = params.D;

m = length(E);
n = length(D);


S_eps_delta_values = y(1:(m*n)); 
I_eps_delta_values = y((m*n+1):(2*m*n)); 
R_eps_delta_values = y((2*m*n+1):(3*m*n));

%reshapeing
S_eps_delta_array = reshape(S_eps_delta_values,m,n);
I_eps_delta_array = reshape(I_eps_delta_values,m,n);
R_eps_delta_array = reshape(R_eps_delta_values,m,n);


%differential equations
dSdt = -bet*sum(sum(D.*I_eps_delta_array)).*E.*S_eps_delta_array;

dIdt = bet*sum(sum(D.*I_eps_delta_array)).*E.*S_eps_delta_array - gam.*I_eps_delta_array; 

dRdt = gam.*I_eps_delta_array;

%reshaping
dSdt = reshape(dSdt,[m*n,1]);
dIdt = reshape(dIdt,[m*n,1]);
dRdt = reshape(dRdt,[m*n,1]);


dydt = [dSdt; dIdt; dRdt];

