function [FOS_traj] = getFinalOutbreakSize(I_traj, R_traj)
    FOS_traj = I_traj + R_traj;
end
