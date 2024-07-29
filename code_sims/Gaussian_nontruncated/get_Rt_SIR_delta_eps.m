function [Rt_calc] = get_Rt_SIR_delta_eps(params,y)

%% Calculation
% calculate effective reproduction number, R_t

% parameters to local variables
bet = params.beta; 
gam = params.gamma; 

for count=1:length(params.t_span)

% t = params.t_span(count);
S = y(count,1);
epsilon_S = y(count,4);
delta_I = y(count,5);

% transmissions
T = bet*S*delta_I*epsilon_S;

% transitions
Sigma = -gam;

NGM = -T*(1/Sigma); % next generation matrix: -T*Sigma^-1

eigen_values = eig(NGM);

Rt_calc(count) = max(eigen_values);

end