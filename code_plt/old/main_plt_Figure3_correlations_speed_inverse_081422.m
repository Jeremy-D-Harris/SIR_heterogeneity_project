% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation

%%

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure3_correlations_speed_inverse_081422';


% light blue, medium blue, dark blue
colors_rgb = [133,192,249; 74,112,188.5; 15,32,128]/255;

% define ending time:
this_t_end_plt = 200;

var_eps = linspace(0.999,0.49,201);

ind_level = 40; %18
epsilon_level = var_eps(ind_level);


%%
% need to load all three cases first:
file_location = '../data/gaussianTruncatedVar1/';

for count = 1:3
    
    if count==1
        % (1) positive
        infile_independent = 'dataPosCinverseHighVar.mat';
        % infile_independent = 'dataIVGaussianSusVarOnlyInverseHighVar.mat';
        %         infile_independent = 'dataPos0.8Figure3.mat';
        load(strcat(file_location,infile_independent));
        
        
    elseif count==2
        % (2) independent
        infile_positivecorrelations = 'dataIndinverseHighVar.mat';
        %         infile_positivecorrelations = 'dataIndFigure3.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
        
    else
        % (3) negative
        infile_negativecorrelations = 'dataNegCinverseHighVar.mat';
        %         infile_negativecorrelations = 'dataNeg-0.8Figure3.mat';
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
    ind = find(this_I_traj_interp>10^-7);
    this_ind_shift_finer = ind(1);
    ind_t_shift(count) = this_ind_shift_finer;
    
    % align trajectories manually too
    this_t_shift_finer = t_span_finer(this_ind_shift_finer)+22;
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
    
    % find mean susceptibility level for each case
    ind = find(this_mu_epsilon_S_traj_interp<epsilon_level);
    this_ind_t_int = ind(1);
    ind_t_int(count) = this_ind_t_int;
    
    this_t_int = this_t_shift_plt_finer(this_ind_t_int);
    t_int(count) = this_t_shift_plt_finer(this_ind_t_int);
    
    this_mu_epsilon = mu_epsilon_S_traj_interp(this_ind_t_int);
    mu_epsilon(count) = this_mu_epsilon;
    
    % delta trajectory and corresponding intercept
    this_mu_delta_I_traj = data.mu_delta_I_traj;
    mu_delta_I_traj(:,count) = this_mu_delta_I_traj;
    
    this_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_mu_delta_I_traj,this_t_shift_plt_finer, 'spline');
    mu_delta_I_traj_interp(:,count) = this_mu_delta_I_traj_interp;
    
    this_mu_delta = this_mu_delta_I_traj_interp(this_ind_t_int);
    mu_delta(count) = this_mu_delta;
    
    this_beta_mu_delta_I_traj = params.beta*this_mu_delta_I_traj;
    beta_mu_delta_I_traj(:,count) = this_beta_mu_delta_I_traj;
    
    this_beta_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_beta_mu_delta_I_traj,this_t_shift_plt_finer,'spline');
    beta_mu_delta_I_traj_interp(:,count) = this_beta_mu_delta_I_traj_interp;
    
    
    %     this_beta_mu_delta = (params.beta*this_mu_epsilon_S_traj(this_ind_shift_finer))*this_mu_delta;
    %     beta_mu_delta(count) = this_beta_mu_delta;
    
    %     transmission_rate(count)=params.beta*this_mu_epsilon_S_traj(this_ind_t_shift);
    
    this_total_incidence = transpose(params.beta*(this_mu_delta_I_traj.*this_I_traj.*this_mu_epsilon_S_traj.*this_S_traj));
    total_incidence(count,:) = this_total_incidence;
    this_total_incidence_interp = interp1(this_t_shift_plt,this_total_incidence,this_t_shift_plt_finer, 'spline');
    total_incidence_interp(count,:) = this_total_incidence_interp;
    %     mu_delta_I_traj = data.mu_delta_I_traj;
    %     beta_mu_delta_I_traj = transmission_rate(count)*mu_delta_I_traj;
    %     beta_mu_delta(count) = beta_mu_delta_I_traj(ind_t_int);
    
    marg_S_e_traj = data.marg_S_e_traj;
    for kk=1:length(params.t_span)
        this_S_marg(1,:) = sum(data.S_traj_discrete(kk,:,:),3);
        marg_S_e_traj(kk,:) = (1/data.S_traj(kk))*this_S_marg;
    end
    
    for counter=1:length(marg_S_e_traj)
        
        this_marg_S_e_traj_interp = interp1(this_t_shift_plt,marg_S_e_traj(:,counter),this_t_shift_plt_finer);
        marg_S_e_traj_interp(:,counter) = this_marg_S_e_traj_interp;
        
    end
    
    this_marginal_susceptibility = marg_S_e_traj_interp(this_ind_t_int,:)';
    marginal_susceptibility(:,count) = this_marginal_susceptibility;
    
    
    marg_I_d_traj = data.marg_I_d_traj;
    for kk=1:length(params.t_span)
        this_I_marg(1,:) = sum(data.I_traj_discrete(kk,:,:),2);
        marg_I_d_traj(kk,:) = (1/data.I_traj(kk))*this_I_marg;
    end
    
    for counter=1:length(marg_I_d_traj)
        
        this_marg_I_d_traj_interp = interp1(this_t_shift_plt,marg_I_d_traj(:,counter),this_t_shift_plt_finer);
        marg_I_d_traj_interp(:,counter) = this_marg_I_d_traj_interp;
        
    end
    
    this_marginal_transmissibility = marg_I_d_traj_interp(this_ind_t_int,:)';
    marginal_transmissibility(:,count) = this_marginal_transmissibility;
    
    % for plotting
    ind = find(params_iv.e > epsilon_level);
    this_ind_eps_int = ind(1);
    ind_eps_int(count) =  this_ind_eps_int;
    
    ind = find(params_iv.e > mu_delta(count));
    this_ind_delta_int = ind(1);
    ind_delta_int(count) = this_ind_delta_int;
    
    % interpolate trajectories too:
    this_S_traj_interp = interp1(this_t_shift_plt,this_S_traj,t_shift_plt_finer,'spline');
    
end


%% Plotting
f1 = figure(1); set(f1, 'Position', [100 500 1200 700]);

%% panel A
subplot(2,3,1);
for count = 1:3
    
    this_p(4-count) = plot(t_shift_plt_finer(4-count,:), total_incidence_interp(4-count,:),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    
end

for count = 1:3
    
    plot(t_shift_plt_finer(4-count,ind_t_int(4-count))*ones(1,10),linspace(10^-7,total_incidence_interp(4-count,ind_t_int(4-count)),10),'--','Color',colors_rgb(4-count,:),'LineWidth',1); hold on;
    
end

axis([0 this_t_end_plt 0 0.02]);
xlabel('Time (days)'); ylabel({'Incident infections, $\eta(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'A'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


legend_char1 = 'Positive Correlation, $\rho > 0$';
legend_char2 = 'No Correlation, $\rho = 0$';
legend_char3 = 'Negative Correlation, $\rho < 0$';

legend(this_p,{legend_char1,legend_char2, legend_char3},'Interpreter','Latex', 'Position', [0.225 0.855 0.095 0.06],'FontSize',9);


%% panel B:
% now plot mean susceptibility
subplot(2,3,2);
for count = 1:3
    
    this_q(4-count)=plot(t_shift_plt_finer(4-count,:), mu_epsilon_S_traj_interp(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    
end

for count = 1:3
    
    plot(t_shift_plt_finer(4-count,ind_t_int(4-count))*ones(1,10),linspace(0.5,mu_epsilon_S_traj_interp(ind_t_int(4-count),4-count),10),'--','Color',colors_rgb(4-count,:),'LineWidth',1); hold on;
    
end

plot(t_shift_plt_finer(count,:),epsilon_level*ones(size(t_shift_plt_finer)),'k','LineWidth',1);
axis([0 this_t_end_plt 0.5 1]);
xlabel('Time (days)'); ylabel({'Mean Susceptiblity, $\bar{\varepsilon}(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'B'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

title('Susceptibility','FontWeight','normal');

%% panel C:
% effective transmission rate
subplot(2,3,3);

for count = 1:3
    
    this_r(4-count)=plot(t_shift_plt_finer(4-count,:),beta_mu_delta_I_traj_interp(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    %     this_r(count).Color(4) = 0.8;
    
end

for count = 1:3
    
    plot(t_shift_plt_finer(4-count,ind_t_int(4-count))*ones(1,10),...
        linspace(0.1,beta_mu_delta_I_traj_interp(ind_t_int(4-count),4-count),10),'--','Color',colors_rgb(4-count,:),'LineWidth',1); hold on;
    plot(t_shift_plt_finer(4-count,:),beta_mu_delta_I_traj_interp(ind_t_int(4-count),4-count)*ones(length(t_shift_plt_finer)),'Color',colors_rgb(4-count,:),'LineWidth',1); hold on;
    
end

axis([0 this_t_end_plt 0.1 0.3]);
xlabel('Time (days)'); ylabel({'Effective Transmission Rate, $\beta\,\bar{\delta}_I(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'C'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

title('Transmissibility','FontWeight','normal');

%% panel D:
% Susceptible population
subplot(2,3,4);
for count = 1:3
    
    this_s(4-count)=plot(t_shift_plt_finer(4-count,:), S_traj_interp(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
end

for count = 1:3
    
    plot(t_shift_plt_finer(4-count,ind_t_int(4-count))*ones(1,10),linspace(0,S_traj_interp(ind_t_int(4-count),4-count),10),'--','Color',colors_rgb(4-count,:),'LineWidth',1); hold on;
    
end

axis([0 this_t_end_plt 0.2 1]);
xlabel('Time (days)'); ylabel({'Susceptible Population, $S(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'D'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


%% panel E:
% now plot marginals for susceptibility
subplot(2,3,5);

for count = 1:3
    
    this_t(count)=plot(params_iv.e,marginal_susceptibility(:,count),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
    
end

plot(epsilon_level,0,'.','Color',colors_rgb(1,:),'LineWidth',2,'MarkerSize',20);
plot(epsilon_level,0,'o','Color',colors_rgb(2,:),'LineWidth',2,'MarkerSize',8);
plot(epsilon_level,0,'o','Color',colors_rgb(3,:),'LineWidth',2,'MarkerSize',12);


axis([0 3 1e-10 1.2e-2]);
yticks([]);
xlabel('Susceptibility, $\varepsilon$','interpreter','latex');
ylabel({'Population Density, $g_S(t_1,\varepsilon)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

title('Susceptibility','FontWeight','normal');

txt = {'E'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

%% panel F:
% now plot marginals for transmissibiilty
subplot(2,3,6);
for count = 1:3
    
    this_u(4-count)=plot(params_iv.e,marginal_transmissibility(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    plot(mu_delta(4-count),0,'.','Color',colors_rgb(4-count,:),'LineWidth',2,'MarkerSize',20);
    
end

axis([0 3 1e-10 1.2e-2]);
yticks([]);
xlabel('Transmissibility, $\delta$','interpreter','latex');
ylabel({'Population Density, $h_I(t_1,\delta)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

title('Transmissibility','FontWeight','normal');

txt = {'F'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


if save_fig_ans==1
    
    figures_location = './figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));
    
end