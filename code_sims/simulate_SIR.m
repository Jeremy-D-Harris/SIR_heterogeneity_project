function dydt = simulate_SIR(t,y,params)

% parameters to local variables
bet=params.beta; 
gam = params.gamma;

S = y(1); I = y(2); R = y(3);

dSdt = -bet*S*I;

dIdt = bet*S*I - gam*I;

dRdt = gam*I;

dydt = [dSdt; dIdt; dRdt];