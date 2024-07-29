function Rt_calc = get_Rt_SIR_eps(params,y)


%% calculate effective reproduction number, R_t
% for susc. variation alone

% parameters to local variables
bet = params.bet;
gam = params.gam;


for count=1:length(params.t_span)

    % t = params.t_span(count);
    this_S_traj = y(count,1);
    this_mean_epsilon_S = y(count,4);
    % this_mean_delta_I = y(count,5);


    % transmissions
    T = bet*this_mean_epsilon_S*this_S_traj;

    % transitions
    Sigma = -gam;

    % next generation matrix: -T*Sigma^-1
    NGM = -T*(1/Sigma);

    eigen_values = eig(NGM);

    Rt_calc(count) = max(eigen_values);


end