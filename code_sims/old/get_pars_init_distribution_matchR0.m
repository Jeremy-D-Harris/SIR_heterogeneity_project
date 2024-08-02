function f = get_pars_init_distribution_matchR0(x,params)

% parameters
gam = params.gam;
eps = params.eps;
del = params.del;
dx = params.dx;
% init_joint_I = params.init_joint_I;
% init_joint_S = params.init_joint_S;

% parameters set
% set_corr_coeff = params.set_corr_coeff;
set_variance_eps = params.set_variance_eps;
set_variance_delta = params.set_variance_delta;

% intended values
intended_mean_eps = params.intended_mean_eps;
intended_mean_delta = params.intended_mean_delta;
intended_variance_eps = params.intended_variance_eps;
intended_variance_delta = params.intended_variance_delta;
intended_corr_coeff = params.intended_corr_coeff;
intended_R0 = params.intended_R0;

% x values
% mean_eps_S = x(1);
mean_eps_S = 1;
% mean_delta_S = x(2);
mean_delta_S = 1;
corr_coeff = x(1);
bet = x(2);

[X1, X2] = meshgrid(eps,del);
X = [X1(:) X2(:)];

mean_S = [mean_eps_S mean_delta_S];

offdiagonal_Sigma = corr_coeff*sqrt(set_variance_eps)*sqrt(set_variance_delta);

Sigma_S = [set_variance_eps, offdiagonal_Sigma; offdiagonal_Sigma, set_variance_delta]; %change off-diagonals for correlation

init_joint_S = mvnpdf(X, mean_S, Sigma_S);
init_joint_S = reshape(init_joint_S,length(del),length(eps))/sum(sum(init_joint_S))/dx/dx;

% marginals
init_marginal_eps_S = dx*sum(init_joint_S);
init_marginal_delta_I = dx*sum(init_joint_S,2)';

% Calculated Means: want to be = 1
mean_eps_S = dx*sum(eps.*init_marginal_eps_S);
mean_delta_S = dx*sum(del.*init_marginal_delta_I);

% Calculated variance in S
calc_variance_eps_S = dx*sum((eps- mean_eps_S*ones(size(eps))).^2.*init_marginal_eps_S);
variance_delta_S = dx*sum((del- mean_delta_S*ones(size(del))).^2.*init_marginal_delta_I);

% Calculated covariance in S
covariance_S = dx*dx*(del- mean_delta_S*ones(size(del)))*init_joint_S*(eps- mean_eps_S*ones(size(eps)))';

% corr_coeff: want vary from:
% -0.6, -0.3, 0.3, 0.6
calc_corr_coeff = covariance_S/sqrt(calc_variance_eps_S)/sqrt(variance_delta_S);

%% calculate R0

% transmissions
% from analysis: equations 16
calc_R0 = (bet/gam)*(mean_eps_S*mean_delta_S + calc_corr_coeff);

%% fuction to minimize
% f  = (mean_eps_S - intended_mean_eps)^2 + (mean_delta_S - intended_mean_delta)^2 + (calc_corr_coeff - intended_corr_coeff)^2;
% f  = (mean_eps_S - intended_mean_eps)^2 + (mean_delta_S - intended_mean_delta)^2 + (calc_corr_coeff - intended_corr_coeff)^2+(calc_variance_eps_S - intended_variance_eps)^2+(variance_delta_S - intended_variance_delta)^2;
% f  = (mean_eps_S - intended_mean_eps)^2 + (mean_delta_S - intended_mean_delta)^2 + (calc_corr_coeff - intended_corr_coeff)^2*(calc_variance_eps_S - intended_variance_eps)^2*(variance_delta_S - intended_variance_delta)^2;

% f  =  (calc_R0 - intended_R0)^2+(mean_eps_S - intended_mean_eps)^2 + (calc_corr_coeff - intended_corr_coeff)^2*(calc_variance_eps_S - intended_variance_eps)^2*(variance_delta_S - intended_variance_delta)^2;
f  =  (calc_R0 - intended_R0)^2 + (calc_corr_coeff - intended_corr_coeff)^2*(calc_variance_eps_S - intended_variance_eps)^2*(variance_delta_S - intended_variance_delta)^2;

