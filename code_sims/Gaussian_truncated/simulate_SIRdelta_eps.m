function dydt = simulate_SIRdelta_eps(t,y,params)

% parameters to local variables
bet=params.beta; 
gam = params.gamma;
D = params.D;
E = params.E;

m = length(E);
n = length(D);


dist_S = y(1:m*n); dist_I = y(m*n+1:2*m*n); dist_R = y(2*m*n+1:3*m*n);

%reshaping
dist_S = reshape(dist_S,m,n);
dist_I = reshape(dist_I, m,n);
dist_R = reshape(dist_R,m,n);


%differential equations
dSdt = -bet*sum(sum(D.*dist_I,2)).*E.*dist_S;

dIdt = bet*sum(sum(D.*dist_I,2)).*E.*dist_S - gam.*dist_I; 

dRdt = gam.*dist_I;

%reshaping
dSdt = reshape(dSdt,[m*n,1]);
dIdt = reshape(dIdt, [m*n,1]);
dRdt = reshape(dRdt, [m*n,1]);


dydt = [dSdt; dIdt; dRdt];

