function Rt_calc = get_Rt_SIR_reducedgamma(params,y)

%% Calculation
% calculate effective reproduction number, R_t
% classic S-I-R model, Rt = (beta/gamma)*S(t)

% parameters to local variables
bet = params.bet;
gam = params.gam;
shape_eps = params.shape_eps;

for count=1:length(params.t_span)

    % t = params.t_span(count);
    S = y(count,1);

    % transmissions: eps_bar(t) = S(t)^1/k for gamma 
    T = bet*S^(1+1/shape_eps)/N;

    % transitions
    Sigma = -gam;

    NGM = -T*(1/Sigma); % next generation matrix: -T*Sigma^-1

    eigen_values = eig(NGM);

    Rt_calc(count) = max(eigen_values);

end