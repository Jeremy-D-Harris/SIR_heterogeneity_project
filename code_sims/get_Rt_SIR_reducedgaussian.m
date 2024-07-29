function Rt_calc = get_Rt_SIR_reducedgaussian(params,y)

%% Calculation
% calculate effective reproduction number, R_t
% classic S-I-R model, Rt = (beta/gamma)*S(t)

% parameters to local variables
bet = params.bet;
gam = params.gam;
N = params.N;
init_mean_eps = params.mean_eps_S;
init_variance_eps = params.variance_eps_S;

for count=1:length(params.t_span)

    % t = params.t_span(count);
    S = y(count,1);
    % mean_eps_S = y(count,4);
    mean_eps_S = 2*init_variance_eps*(log(S)+init_mean_eps^2/2/init_variance_eps-log(N));

    % transmissions: eps_bar(t) = S(t)^1/k for gamma 
    T = bet*S*mean_eps_S;

    % transitions
    Sigma = -gam;

    NGM = -T*(1/Sigma); % next generation matrix: -T*Sigma^-1

    eigen_values = eig(NGM);

    Rt_calc(count) = max(eigen_values);

end