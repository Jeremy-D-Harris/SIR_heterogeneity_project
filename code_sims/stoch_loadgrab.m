%utility function to ascertain parameters to initialize stochastic
%simulations
function [bet,gam,susceptibility,transmissibility,initjoint] = stoch_loadgrab(filename)
    load(filename,"params","results");
    bet = params.bet;
    gam = params.gam;
    susceptibility = params.eps;
    transmissibility = params.del;
    initjoint = squeeze(results.S_traj_eps_delta_array(1,:,:)); %epsilon columns, delta on rowss
end