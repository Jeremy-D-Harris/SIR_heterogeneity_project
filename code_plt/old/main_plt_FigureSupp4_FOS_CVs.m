% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation

%%

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'FigureS4_fosandcvs_080922';

% light blue, medium blue, dark blue, gray, black
colors_rgb = [133,192,249; 74,112,188.5; 15,32,128; 166 166 166 ; 0,0,0]/255;
% define ending time:  
this_t_end_plt = 300;


%%
% need to load all three cases first:
file_location = '../data/gaussianTruncatedVar1/';

for count = 1:3
    
        
    if count==1
        % (1) positive
        infile_positivecorrelations = 'dataPosCinverseHighVar.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
        
    elseif count==2
        % (2) independent
        infile_independent = 'dataIndinverseHighVar.mat';
        load(strcat(file_location,infile_independent));
        infile_svonly = 'dataIVGaussianSusVarOnlyInverseHighVar';
        variable = 'data';
        data_susvaronly= load(strcat(file_location,infile_svonly),variable);
        FOSed_traj_susvaronly = data_susvaronly.data.FOSed_traj;
        mu_epsilon_S_traj_susvaronly = data_susvaronly.data.mu_epsilon_S_traj;
        marg_S_e_traj_susvaronly = data_susvaronly.data.marg_S_e_traj;
        for kk=1:length(params.t_span)
        this_S_marg_susvaronly(1,:) = sum(data_susvaronly.data.S_traj_discrete(kk,:,:),3);
        marg_S_e_traj_susvaronly(kk,:) = (1/data_susvaronly.data.S_traj(kk))*this_S_marg_susvaronly;
    end

        
        
    else
        % (3) negative
        infile_negativecorrelations = 'dataNegCinverseHighVar.mat';
        load(strcat(file_location,infile_negativecorrelations));
        
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
    ind = find(this_I_traj_interp>10^-5);
    this_ind_shift_finer = ind(1);
    ind_t_shift(count) = this_ind_shift_finer;
    
    % align trajectories manually too
    this_t_shift_finer = t_span_finer(this_ind_shift_finer)+0;
    %     if count==1
    %         this_t_shift_finer = t_span_finer(this_ind_shift_finer);
    %     elseif count==2
    %         this_t_shift_finer = t_span_finer(this_ind_shift_finer);
    %     else
    %         this_t_shift_finer = t_span_finer(this_ind_shift_finer);
    %     end
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
    %     this_mu_epsilon_S_traj = this_mu_epsilon_S_traj/this_mu_epsilon_S_traj(this_ind_shift_finer);
    %     scaled_mu_epsilon_S_traj(:,count) = this_mu_epsilon_S_traj;
    
    this_mu_epsilon_S_traj_interp = interp1(this_t_shift_plt,this_mu_epsilon_S_traj,this_t_shift_plt_finer,'spline');
    mu_epsilon_S_traj_interp(:,count) = this_mu_epsilon_S_traj_interp;
    
    % mean potential transmissibility over time
    this_mu_delta_S_traj = data.mu_delta_S_traj;
    mu_delta_S_traj(:,count) = this_mu_delta_S_traj;
    %     this_mu_epsilon_S_traj = this_mu_epsilon_S_traj/this_mu_epsilon_S_traj(this_ind_shift_finer);
    %     scaled_mu_epsilon_S_traj(:,count) = this_mu_epsilon_S_traj;
    
    this_mu_delta_S_traj_interp = interp1(this_t_shift_plt,this_mu_delta_S_traj,this_t_shift_plt_finer,'spline');
    mu_delta_S_traj_interp(:,count) = this_mu_delta_S_traj_interp;
    
    %     %% find mean susceptibility level for each case
    %     ind = find(this_mu_epsilon_S_traj_interp<epsilon_level);
    %     this_ind_t_int = ind(1);
    %     ind_t_int(count) = this_ind_t_int;
    %
    %     this_t_int = this_t_shift_plt_finer(this_ind_t_int);
    %     t_int(count) = this_t_shift_plt_finer(this_ind_t_int);
    %
    %     this_mu_epsilon = mu_epsilon_S_traj_interp(this_ind_t_int);
    %     mu_epsilon(count) = this_mu_epsilon;
    
    % delta trajectory and corresponding intercept
    this_mu_delta_I_traj = data.mu_delta_I_traj;
    mu_delta_I_traj(:,count) = this_mu_delta_I_traj;
    
    this_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_mu_delta_I_traj,this_t_shift_plt_finer,'spline');
    mu_delta_I_traj_interp(:,count) = this_mu_delta_I_traj_interp;
    
    %     this_mu_delta = this_mu_delta_I_traj_interp(this_ind_t_int);
    %     mu_delta(count) = this_mu_delta;
    
    this_beta_mu_delta_I_traj = params.beta*this_mu_delta_I_traj;
    beta_mu_delta_I_traj(:,count) = this_beta_mu_delta_I_traj;
    
    this_beta_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_beta_mu_delta_I_traj,this_t_shift_plt_finer,'spline');
    beta_mu_delta_I_traj_interp(:,count) = this_beta_mu_delta_I_traj_interp;
    
    
    %     this_beta_mu_delta = (params.beta*this_mu_epsilon_S_traj(this_ind_shift_finer))*this_mu_delta;
    %     beta_mu_delta(count) = this_beta_mu_delta;
    
    %     transmission_rate(count)=params.beta*this_mu_epsilon_S_traj(this_ind_t_shift);
    
    this_total_incidence = transpose(params.beta*(this_mu_delta_I_traj.*this_I_traj.*this_mu_epsilon_S_traj.*this_S_traj));
    total_incidence(count,:) = this_total_incidence;
    this_total_incidence_interp = interp1(this_t_shift_plt,this_total_incidence,this_t_shift_plt_finer,'spline');
    total_incidence_interp(count,:) = this_total_incidence_interp;
    
    this_cumulative_infections = data.FOSed_traj';
    cumulative_infections(count,:) = this_cumulative_infections;
    
    final_outbreak_size(count) = this_cumulative_infections(end);
    
    if count==3
        final_outbreak_size(4)=data.FOS_SIR_traj(end);
        final_outbreak_size(5)=FOSed_traj_susvaronly(end);
    end
    
    
    %homogeneous
    this_I_traj_homogeneous = data.I_traj_v;
    this_S_traj_homogeneous = data.S_traj_v;
    this_total_incidence_homogeneous = params_v.beta*(this_I_traj_homogeneous.*this_S_traj_homogeneous);
    
    
    marg_S_e_traj = data.marg_S_e_traj;
    for kk=1:length(params.t_span)
        this_S_marg(1,:) = sum(data.S_traj_discrete(kk,:,:),3);
        marg_S_e_traj(kk,:) = (1/data.S_traj(kk))*this_S_marg;
    end
    
    for counter=1:length(marg_S_e_traj)
        
        this_marg_S_e_traj_interp = interp1(this_t_shift_plt,marg_S_e_traj(:,counter),this_t_shift_plt_finer,'spline');
        marg_S_e_traj_interp(:,counter) = this_marg_S_e_traj_interp;
        
    end
    
    %this_marginal_susceptibility = marg_S_e_traj_interp(this_ind_t_int,:)';
    %marginal_susceptibility(:,count) = this_marginal_susceptibility;
    
    
    marg_I_d_traj = data.marg_I_d_traj;
    for kk=1:length(params.t_span)
        this_I_marg(1,:) = sum(data.I_traj_discrete(kk,:,:),2);
        marg_I_d_traj(kk,:) = (1/data.I_traj(kk))*this_I_marg;
    end
    
    for counter=1:length(marg_I_d_traj)
        
        this_marg_I_d_traj_interp = interp1(this_t_shift_plt,marg_I_d_traj(:,counter),this_t_shift_plt_finer,'spline');
        marg_I_d_traj_interp(:,counter) = this_marg_I_d_traj_interp;
        
    end
    
    
    
    
    %this_marginal_transmissibility = marg_I_d_traj_interp(this_ind_t_int,:)';
    %marginal_transmissibility(:,count) = this_marginal_transmissibility;
    
    % for plotting
    %     ind = find(params_iv.e > epsilon_level);
    %     this_ind_eps_int = ind(1);
    %     ind_eps_int(count) =  this_ind_eps_int;
    %
    %     ind = find(params_iv.e > mu_delta(count));
    %     this_ind_delta_int = ind(1);
    %     ind_delta_int(count) = this_ind_delta_int;
    %
    %     % interpolate trajectories too:
    %     this_S_traj_interp = interp1(this_t_shift_plt,this_S_traj,t_shift_plt_finer);
    %     this_I_traj_interp = interp1(this_t_shift_plt,this_I_traj,t_shift_plt_finer);
    %
    %calculate correlation coefficients
    % joint distribution at t0% ~t1
    joint_S_eps_delta(:,:) = data.S_traj_discrete(1,:,:);
    
    
    marg_S_d_traj = data.marg_S_d_traj;
    for kk=1:length(params.t_span)
        this_S_marg(1,:) = sum(data.S_traj_discrete(kk,:,:),2);
        marg_S_d_traj(kk,:) = (1/data.S_traj(kk))*this_S_marg;
    end
    
    %calculate correlation coefficients
    cov = sum(sum(((params.E(:,1)- mu_epsilon_S_traj(1)*ones(size(params.E(:,1))))*(params.D(1,:)- mu_delta_S_traj(1)*ones(size(params.D(1,:))))).*joint_S_eps_delta),2);
    var_e = sum((params.E(:,1)'- mu_epsilon_S_traj(1)*ones(size(params.E(:,1)'))).^2.*marg_S_e_traj(1,:));
    sd_e = sqrt(var_e);
    var_d = sum((params.D(1,:)- mu_delta_S_traj(1)*ones(size(params.D(1,:)))).^2.*marg_S_d_traj(1,:),2);
    sd_d = sqrt(var_d);
    corrcoef(:,count) = cov/(sd_e*sd_d);
    
    
    %% calculate variances
    if count == 1
        Svar = zeros(length(params.t_span),3);
        Sdvar = zeros(length(params.t_span),3);
        Ivar = zeros(length(params.t_span),3);
    end
    
    
    for tt = 1:length(params.t_span)
        Svar(tt,count) = sum((params.E(:,1)'- this_mu_epsilon_S_traj(tt)*ones(size(params.E(:,1)'))).^2.*marg_S_e_traj(tt,:));
    end

    if count == 3
        Svar_susvaronly = zeros(length(params.t_span),1);
        for tt = 1:length(params.t_span)
            Svar_susvaronly(tt) = sum((params.E(:,1)'- mu_epsilon_S_traj_susvaronly(tt)*ones(size(params.E(:,1)'))).^2.*marg_S_e_traj_susvaronly(tt,:));
        end
    end
    

    
    for tt = 1:length(params.t_span)
        Sdvar(tt,count) = sum((params.D(1,:)- this_mu_delta_S_traj(tt)*ones(size(params.D(1,:)))).^2.*marg_S_d_traj(tt,:));
    end

    
    for tt = 1:length(params.t_span)
        Ivar(tt,count) = sum((params.D(1,:)- this_mu_delta_I_traj(tt)*ones(size(params.D(1,:)))).^2.*marg_I_d_traj(tt,:));
    end

    
end





%% Plotting
f1 = figure(1); set(f1, 'Position', [100 500 1200 350]);

% f1 = figure(1); set(f1, 'Position', [100 500 1200 350]);

%%
for count = 1:3
    
    
    %% FOS
    subplot(1,3,1);
    this_h(count+2) = plot(t_shift_plt(count,:), cumulative_infections(count,:),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
    axis([0 this_t_end_plt 0 1]);
    xlabel('Time (days)'); %ylabel('Cumulative Infections');
    ylabel({'Cumulative Infections, $\int_0^t I(s) \, ds$ '},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if count==3
        
        
        txt = {'A'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
        
        this_h(1) = plot(t_shift_plt(2,:), data.FOS_SIR_traj,'Color',colors_rgb(4,:),'LineWidth',2.5); hold on;
        this_h(2) =plot(t_shift_plt(2,:),  FOSed_traj_susvaronly,'--','Color',colors_rgb(5,:),'LineWidth',2.5); hold on;
        yticks(0:0.2:1);
        set(f1,'yticklabel',[{'0'},{'0.2'},{'0.4'},{'0.6'},{'0.8'},{'1.0'}]);
        %         set(f1,'xticklabel',{[]},'yticklabel',[{'0'},{''},{'0.2'},{''},{'0.4'},{''},{'0.6'},{''},{'0.8'},{''},{'1.0'}]);
        
        legend_char1 = 'Classic SIR';
        legend_char2 = 'Variation in Susceptibility';
        legend_char5 = 'Negative Correlation, $\rho < 0$';
        legend_char4 = 'No Correlation, $\rho = 0$';
        legend_char3 = 'Positive Correlation, $\rho > 0$';
        
        
        legend(this_h,{legend_char1,legend_char2,legend_char3,legend_char4, legend_char5}, 'Interpreter','Latex','Location',[0.226 0.155 0.1 0.2],'FontSize',10);
        
        
    end
    
    %% CV Susceptibility
    subplot(1,3,2);
    
    this_p(count) = plot(t_shift_plt(count,:),Svar(:,count)./(mu_epsilon_S_traj(:,count).^2), 'Color', colors_rgb(count,:), 'LineWidth',2);hold on;
    
    axis([0 this_t_end_plt 0 1]);
    xlabel({'Time (days)'},'Interpreter','Latex'); 
    ylabel({'Coefficient of Variation (Squared)'},'Interpreter','Latex');
    title('Susceptibility');
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if count==3
        this_p(4) = plot(t_shift_plt(2,:),Svar_susvaronly./(mu_epsilon_S_traj_susvaronly.^2),'--','Color',colors_rgb(5,:),'LineWidth',2);
        txt = {'B'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
        box on
        
        
%         legend_char1 = 'Positive Correlation, $\rho > 0$';
%         legend_char2 = 'No Correlation, $\rho = 0$';
%         legend_char3 = 'Negative Correlation, $\rho < 0$';
%         legend_char4 = 'Variation in Susceptibility';
%         
%         legend(this_p,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','NorthWest','FontSize',10);
        
        
    
        
    end
    
    
    %% CV Transmissibility
    subplot(1,3,3);
    
    
    this_q(count) = plot(t_shift_plt(count,:),Ivar(:,count)./(mu_delta_I_traj(:,count).^2), 'Color', colors_rgb(count,:), 'LineWidth',2);hold on;
    
    axis([0 this_t_end_plt 0 1]);
    xlabel({'Time (days)'},'Interpreter','Latex'); 
    ylabel({'Coefficient of Variation (Squared)'},'Interpreter','Latex');
    title('Transmissibility');
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if count==3
        txt = {'C'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
        box on
        
        
%         legend_char1 = 'Positive Correlation, $\rho > 0$';
%         legend_char2 = 'No Correlation, $\rho = 0$';
%         legend_char3 = 'Negative Correlation, $\rho < 0$';
%         
%         legend(this_q,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthWest','FontSize',10);
        
        
    
        
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