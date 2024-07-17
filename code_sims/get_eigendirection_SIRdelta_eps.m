function eigen_vector = get_eigendirection_SIRdelta_eps(params)


%% parameters
% parameters to local variables
bet=params.beta; 
gam = params.gamma;

epsilon_S = params.mu_eps_S; 
delta_I = params.mu_delta_I; 


%% SIR Version
% note - using this version
% disease free state
S = 1; I = 0; R = 0; %at t = 0


dSdt = [-bet*epsilon_S*delta_I*I,-bet*epsilon_S*delta_I*S,0];

dIdt = [bet*epsilon_S*delta_I*I,bet*epsilon_S*delta_I*S - gam,0];

dRdt = [0,gam,0];


A = [dSdt;dIdt;dRdt];

[eigen_directions, eigen_values] = eig(A); % get eigenvalues/eigenvectors
% eigen_directions
% eigen_values
[val, ind] = max(diag(eigen_values));
eigen_vector = eigen_directions(:,ind);% corresponds to max eigenvalue
end