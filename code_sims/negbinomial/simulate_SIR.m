function dydt = simulate_SIR(t,y,params)

%% for classic SIR model
% no variation!

% parameters to local variables
bet= params.bet; 
gam = params.gam;

S = y(1); I = y(2); R = y(3);

% incident infections
incident_infections = bet*S*I;

% recovery
recovery = gam*I;

dSdt = -incident_infections;

dIdt = incident_infections - recovery;

dRdt = recovery;

dydt = [dSdt; dIdt; dRdt];