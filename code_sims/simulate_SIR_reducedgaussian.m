function dydt = simulate_SIR_reducedgaussian(t,y,params)

%% for SIR model reduced model
% assumes delta_bar is constant = 1
% variance is set equal to the initial

% parameters to local variables
bet= params.bet; 
gam = params.gam;
init_variance_eps = params.variance_eps_S;

S = y(1); I = y(2); R = y(3); mean_eps = y(4);

% incident infections
incident_infections = bet*I*mean_eps*S;

% recovery
recovery = gam*I;

dSdt = -incident_infections;

dIdt = incident_infections - recovery;

dRdt = recovery;

dmean_epsdt = - bet*I*init_variance_eps;

dydt = [dSdt; dIdt; dRdt; dmean_epsdt];