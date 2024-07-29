function [eigen_vector] = get_eigendirection_mean_eps_delta(params)


%% parameters
% parameters to local variables
bet=params.bet; 
gam = params.gam;
N = params.N;

mean_eps_S = params.mean_eps_S; 
mean_delta_I = params.mean_delta_I; 


%% SIR Version
% note - using this version
% disease free state
S = N; I = 0; R = 0; %at t = 0


dSdt = [-bet*mean_eps_S*mean_delta_I*I,-bet*mean_eps_S*mean_delta_I*S,0];

dIdt = [bet*mean_eps_S*mean_delta_I*I,bet*mean_eps_S*mean_delta_I*S - gam,0];

dRdt = [0,gam,0];


A = [dSdt;dIdt;dRdt];

[eigen_directions, eigen_values] = eig(A); % get eigenvalues/eigenvectors
% eigen_directions
% eigen_values
[val, ind] = max(abs(diag(eigen_values)));
eigen_vector = eigen_directions(:,ind);% corresponds to max eigenvalue
end