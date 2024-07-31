% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation

%%

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'NewFigure2_strength_forward_020124';

%create color gradiets
c2 = [133,192,249]/255; % light blue
c1 = [15,32,128]/255; % dark blue

depth = 9;
[grad1,im]=colorGradient(c2,c1,depth);
% [grad2,im]=colorGradient(c2,c3,depth);

colors_rgb = grad1;

% define ending time:  
this_t_end_plt = 300;


%%
% need to load all three cases first:
file_location = '../data/forward/';

for count = 1:5
    
    if count==1
        % (1) positive
        infile_independent = 'dataPosCforwardHighVar.mat';
        load(strcat(file_location,infile_independent));
        
    elseif count==2
        
        % (2) positive
        infile_positivecorrelations = 'dataPosC0.4forwardHighVar.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
    elseif count==3
        % (2) independent
        infile_positivecorrelations = 'dataIndforwardHighVar.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
    elseif count==4
        % (2) negative
        infile_positivecorrelations = 'dataNegC-0.4forwardHighVar.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
    else
        % (2) negative
        infile_positivecorrelations = 'dataNegCforwardHighVar.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
    end
    
    this_I_traj = data.I_traj;
    I_traj(:,count) = this_I_traj;
    
    this_S_traj = data.S_traj;
    S_traj(:,count) = this_S_traj;
    
    % finer mesh for finding intercepts
    n_time_pts = 401;
    t_span_finer = linspace(params.t_span(1),params.t_span(end),n_time_pts);
    this_I_traj_interp = interp1(params.t_span,this_I_traj,t_span_finer,'spline');
    I_traj_interp(:,count) = this_I_traj_interp;
    
    % align trajectories: find initial time point
    ind = find(this_I_traj_interp>10^-7);
    this_ind_shift_finer = ind(1);
    ind_t_shift(count) = this_ind_shift_finer;
    this_t_shift_finer = t_span_finer(this_ind_shift_finer);
    t_shift_finer(count) = this_t_shift_finer;
    
    this_t_shift_plt = params.t_span-this_t_shift_finer*ones(size(params.t_span));
    t_shift_plt(count,:) = this_t_shift_plt;
    
    % finer mesh for t plot
    this_t_shift_plt_finer = t_span_finer-this_t_shift_finer*ones(size(t_span_finer));
    t_shift_plt_finer(count,:) = this_t_shift_plt_finer;
    
    % now interpolate trajectories
    this_S_traj_interp = interp1(this_t_shift_plt,this_S_traj,this_t_shift_plt_finer,'spline');
    S_traj_interp(:,count) = this_S_traj_interp;
    
    % mean susceptibility over time
    this_mu_epsilon_S_traj = data.mu_epsilon_S_traj;
    mu_epsilon_S_traj(:,count) = this_mu_epsilon_S_traj;
    this_mu_epsilon_S_traj_interp = interp1(this_t_shift_plt,this_mu_epsilon_S_traj,this_t_shift_plt_finer,'spline');
    mu_epsilon_S_traj_interp(:,count) = this_mu_epsilon_S_traj_interp;
    
    % mean potential transmissibility over time
    this_mu_delta_S_traj = data.mu_delta_S_traj;
    mu_delta_S_traj(:,count) = this_mu_delta_S_traj;
    this_mu_delta_S_traj_interp = interp1(this_t_shift_plt,this_mu_delta_S_traj,this_t_shift_plt_finer,'spline');
    mu_delta_S_traj_interp(:,count) = this_mu_delta_S_traj_interp;
    
    % delta trajectory and corresponding intercept
    this_mu_delta_I_traj = data.mu_delta_I_traj;
    mu_delta_I_traj(:,count) = this_mu_delta_I_traj;
    this_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_mu_delta_I_traj,this_t_shift_plt_finer,'spline');
    mu_delta_I_traj_interp(:,count) = this_mu_delta_I_traj_interp;
    
    
    this_total_incidence = transpose(params.beta*(this_mu_delta_I_traj.*this_I_traj.*this_mu_epsilon_S_traj.*this_S_traj));
    total_incidence(count,:) = this_total_incidence;
    this_total_incidence_interp = interp1(this_t_shift_plt,this_total_incidence,this_t_shift_plt_finer,'spline');
    total_incidence_interp(count,:) = this_total_incidence_interp;
    
    
    
    marg_S_e_traj = data.marg_S_e_traj;
    for kk=1:length(params.t_span)
        this_S_marg(1,:) = sum(data.S_traj_discrete(kk,:,:),3);
        marg_S_e_traj(kk,:) = (1/data.S_traj(kk))*this_S_marg;
    end
    
    for counter=1:length(marg_S_e_traj)
        
        this_marg_S_e_traj_interp = interp1(this_t_shift_plt,marg_S_e_traj(:,counter),this_t_shift_plt_finer,'spline');
        marg_S_e_traj_interp(:,counter) = this_marg_S_e_traj_interp;
        
    end
    
    
    marg_I_d_traj = data.marg_I_d_traj;
    for kk=1:length(params.t_span)
        this_I_marg(1,:) = sum(data.I_traj_discrete(kk,:,:),2);
        marg_I_d_traj(kk,:) = (1/data.I_traj(kk))*this_I_marg;
    end
    
    for counter=1:length(marg_I_d_traj)
        
        this_marg_I_d_traj_interp = interp1(this_t_shift_plt,marg_I_d_traj(:,counter),this_t_shift_plt_finer,'spline');
        marg_I_d_traj_interp(:,counter) = this_marg_I_d_traj_interp;
        
    end
    
    
    %calculate correlation coefficients
    % joint distribution at t0% ~t1
    joint_S_eps_delta(:,:) = data.S_traj_discrete(1,:,:);
    
    
    marg_S_d_traj = data.marg_S_d_traj;
    for kk=1:length(params.t_span)
        this_S_marg(1,:) = sum(data.S_traj_discrete(kk,:,:),2);
        marg_S_d_traj(kk,:) = (1/data.S_traj(kk))*this_S_marg;
    end
    
    %calculate correlation coefficients
    this_cov = sum(sum(((params.E(:,1)- mu_epsilon_S_traj(1)*ones(size(params.E(:,1))))*(params.D(1,:)- mu_delta_S_traj(1)*ones(size(params.D(1,:))))).*joint_S_eps_delta),2);
    this_var_e = sum((params.E(:,1)'- mu_epsilon_S_traj(1)*ones(size(params.E(:,1)'))).^2.*marg_S_e_traj(1,:));
    this_sd_e = sqrt(this_var_e);
    init_stand_dev_eps(1,count) = this_sd_e;
    
    this_var_d = sum((params.D(1,:)- mu_delta_S_traj(1)*ones(size(params.D(1,:)))).^2.*marg_S_d_traj(1,:),2);
    this_sd_d = sqrt(this_var_d);
    init_stand_dev_delta(1,count) = this_sd_d;
    
    this_corrcoef = this_cov/(this_sd_e*this_sd_d);
    corrcoef(:,count) = this_corrcoef;
    
    R0_calc(1,count) = params.beta*this_mu_epsilon_S_traj(1)*this_mu_delta_I_traj(6)/params.gamma;
    
    R0_analytic(1,count) = params.beta*(this_mu_epsilon_S_traj(1)*this_mu_delta_S_traj(1)+this_corrcoef*this_sd_e*this_sd_d)/params.gamma;
    
end

% corrcoef = flip(corrcoef);
% corrcoef_vary = -0.65:0.05:0.65;
corrcoef_vary = -1:0.05:1;

mean_stand_dev_eps = mean(init_stand_dev_eps);
mean_stand_dev_delta = mean(init_stand_dev_delta);

R0List_analytic = (params.beta/params.gamma)*(mu_epsilon_S_traj(1)*mu_delta_S_traj(1)+corrcoef_vary*mean_stand_dev_eps*mean_stand_dev_delta);
% R0List_analytic = R0List_analytic;



%% Plotting
f1 = figure(1); set(f1, 'Position', [100 500 900 350]);

%%
for count = 1:5
    
    
    % %% cumulative infections
    % subplot(1,2,3);
    % this_p(count) =semilogy(t_shift_plt_finer(count,:),total_incidence_interp(count,:),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5); hold on;
    % %     this_p(count).Color(4)=0.8;
    % axis([0 300 10^-7 0.3*10^-1]);
    % xlabel('Time (days)');
    % ylabel({'Incident Infections, $\eta(t)$'},'Interpreter','Latex');
    % f1=gca;
    % f1.LineWidth = 1;
    % f1.FontSize = 16;
    % f1.FontWeight = 'normal';
    % f1.FontName = 'Times';
    
    
    % if count==5
    % 
    % 
    %     txt = {'A'};
    %     text(0.025,1.035,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    % 
    % 
    %     % xticks([0 100 200 300]);
    %     % yticks([10^-7 10^-6 10^-5 10^-4 10^-3 10^-2]);
    %     % set(f1,'yticklabel',[{'10^{-7}'},{'10^{-6}'},{'10^{-5}'},{'10^{-4}'},{'10^{-3}'},{'10^{-2}'},{''}]);
    % 
    % 
    % end
    
    %%
    subplot(1,2,2);
    
    if count==1
        this_q(1) = plot(corrcoef_vary, R0List_analytic,'k--','LineWidth', 2.5); hold on;
        %         plot(corrcoef, this_R0_analytic,'k.','MarkerSize', 20); hold on;
    end
    plot(corrcoef(:,count), R0_calc(count),'.', 'Color', colors_rgb(2*count-1,:), 'MarkerSize', 30);hold on;
    
    axis([-1 1 1.2 2.8]);
    xlabel({'Correlation Coefficient, $\rho$'},'Interpreter','Latex');
    ylabel({'Basic Reproduction Number, $\mathcal R_0$'},'Interpreter','Latex');
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if count==5
        txt = {'B'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',18,'FontWeight','bold');
        
        
        % legend_char_q = strcat('Eq. (14), $\sigma_\varepsilon(0) = \sigma_\delta(0) = $',num2str(mean_stand_dev_delta,'%2.2f'));
        legend_char_q = 'Analytic $\mathcal R_0$';
        legend(this_q,{legend_char_q}, 'Interpreter','Latex','FontSize',16, 'Location','NorthWest');
        legend boxoff;
        
    end
    
    
    %% now plot linear scale
    subplot(1,2,1);
    
    this_h(count) =plot(t_shift_plt_finer(count,:),total_incidence_interp(count,:),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5); hold on;
    %     this_h(count).Color(4)=0.8;
    %     semilogy(t_shift_plt_finer(count,ind_t_int(count))*ones(length(t_shift_plt_finer)),...
    %         linspace(10^-4,total_incidence_interp(count,ind_t_int(count)),length(t_shift_plt_finer)),'--','Color',colors_rgb(count,:),'LineWidth',2); hold on;
    axis([0 this_t_end_plt 0 0.025]);
    xlabel('Time (days)'); ylabel({'Incident Infections, $\eta(t)$'},'Interpreter','Latex');
    % ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if count==5
        
        txt = {'A'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',18,'FontWeight','bold');
        box on
        
        legend_char1 = strcat('$\rho = $',num2str(round(corrcoef(:,1),2)));
        legend_char2 = strcat('$\rho = $',num2str(round(corrcoef(:,2),2)));
        legend_char3 = strcat('$\rho = $',num2str(round(corrcoef(:,3),2)));
        legend_char4 = strcat('$\rho = $',num2str(round(corrcoef(:,4),2)));
        legend_char5 = strcat('$\rho = $',num2str(round(corrcoef(:,5),2)));
        %legend_char4 = 'Classic SIR';
        
        legend(this_h,{legend_char1,legend_char2,legend_char3,legend_char4,legend_char5}, 'Interpreter','Latex','FontSize',16, 'Location','Northeast');
        
    end
    %
    
    
end


if save_fig_ans==1
    
    figures_location = './figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));
    
end