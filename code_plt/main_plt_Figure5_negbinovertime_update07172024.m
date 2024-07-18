% % plot marginal distributions over time

%%

clear all; close all; clc;

% colors_rgb = [133,192,249; 74,112,188.5; 15,32,128; 166 166 166]/255;
% colors_rgb = [133,192,249; 74,112,188.5; 15,32,128]/255;
colors_rgb = [74,112,188.5; 133,192,249]/255;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure5_negbinovertime_07172024';

% define ending time:
this_t_end_plt = 200;

% spacing/scaling subplots
frac_spacing = 0.74;
frac_scaling = 0.21;

% var_eps = fliplr([0.66 0.74 0.82 0.895]); %linspace(0.999,0.49,201);
var_eps = fliplr([0.3 0.5 0.70 0.95]);

txt_eps1 = ['$\bar{\varepsilon}(t_0) = ', num2str(var_eps(1),'%1.2f'), '$ '];
txt_eps2 = ['$\bar{\varepsilon}(t_1) = ', num2str(var_eps(2),'%1.2f'), '$ '];
txt_eps3 = ['$\bar{\varepsilon}(t_2) = ', num2str(var_eps(3),'%1.2f'), '$ '];
txt_eps4 = ['$\bar{\varepsilon}(t_3) = ', num2str(var_eps(4),'%1.2f'), '$ '];


%%
% need to load all three cases first:
% which_results
file_location = '../data/negbinomial/';

for count = 1:3
    
    if count==1
        
        infile_independent = 'negbinomial_overdispersed_uncorr_update07172024.mat';
        
        % infile_independent = '.mat';
        load(strcat(file_location,infile_independent));
        
        
    elseif count==2
        
        infile_positivecorrelations = 'negbinomial_overdispersed_poscorr_update07172024.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
        
    else
        % get classic
        infile_negativecorrelations = 'negbinomial_Poissonlike.mat';
        load(strcat(file_location,infile_negativecorrelations));
        
    end
    
    this_I_traj = results.I_traj';
    I_traj(count,:) = this_I_traj;
    
    this_S_traj = results.S_traj';
    S_traj(count,:) = this_S_traj;
    
    % finer mesh for finding intercepts
    n_time_pts = 401;
    t_span_finer = linspace(params.t_span(1),params.t_span(end),n_time_pts);
    this_I_traj_interp = interp1(params.t_span,this_I_traj,t_span_finer,'spline');
    I_traj_interp(count,:) = this_I_traj_interp;
    
    % align trajectories: find initial time point
    ind = find(this_I_traj_interp>10^-8);
    this_ind_shift_finer = ind(1);
    ind_t_shift(count) = this_ind_shift_finer;
    
    % align trajectories manually too
    this_t_shift_finer = t_span_finer(this_ind_shift_finer);
    %     if count==1
    %         this_t_shift_finer = t_span_finer(this_ind_shift_finer)+38;
    %     elseif count==2
    %         this_t_shift_finer = t_span_finer(this_ind_shift_finer)+38.25;
    %     else
    %         this_t_shift_finer = t_span_finer(this_ind_shift_finer)+37.75;
    %     end
    t_shift_finer(count) = this_t_shift_finer;
    
    this_t_shift_plt = params.t_span-this_t_shift_finer*ones(size(params.t_span));
    t_shift_plt(count,:) = this_t_shift_plt;
    
    % finer mesh for t plot
    this_t_shift_plt_finer = t_span_finer-this_t_shift_finer*ones(size(t_span_finer));
    t_shift_plt_finer(count,:) = this_t_shift_plt_finer;
    
    % now interpolate trajectories
%     this_S_traj_interp = interp1(this_t_shift_plt,this_S_traj,this_t_shift_plt_finer,'spline');
%     S_traj_interp(count,:) = this_S_traj_interp;
    
    % mean susceptibility over time
    this_mu_epsilon_S_traj = results.mu_epsilon_S_traj;
    mu_epsilon_S_traj(count,:) = this_mu_epsilon_S_traj;
    %     this_mu_epsilon_S_traj = this_mu_epsilon_S_traj/this_mu_epsilon_S_traj(this_ind_shift_finer);
    %     scaled_mu_epsilon_S_traj(:,count) = this_mu_epsilon_S_traj;
    
%     this_mu_epsilon_S_traj_interp = interp1(this_t_shift_plt,this_mu_epsilon_S_traj,this_t_shift_plt_finer,'spline');
%     mu_epsilon_S_traj_interp(count,:) = this_mu_epsilon_S_traj_interp;
    
    
    % delta trajectory and corresponding intercept
    this_mu_delta_I_traj = results.mu_delta_I_traj;
    mu_delta_I_traj(count,:) = this_mu_delta_I_traj;
    
    this_var_delta_I_traj = results.var_delta_I_traj;
    var_delta_I_traj(count,:) = this_var_delta_I_traj;
    
    this_disp_delta_I_traj = this_var_delta_I_traj.^2./this_mu_delta_I_traj;
    disp_delta_I_traj(count,:) = this_disp_delta_I_traj;
    
%     this_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_mu_delta_I_traj,this_t_shift_plt_finer,'spline');
%     mu_delta_I_traj_interp(count,:) = this_mu_delta_I_traj_interp;
    %     collect_mu_delta_I_traj_interp(count_eps,:,count) = this_mu_delta_I_traj_interp;
    
    
    this_beta_mu_delta_I_traj = params.beta*this_mu_delta_I_traj;
    beta_mu_delta_I_traj(count,:) = this_beta_mu_delta_I_traj;
    
%     this_beta_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_beta_mu_delta_I_traj,this_t_shift_plt_finer,'spline');
%     beta_mu_delta_I_traj_interp(count,:) = this_beta_mu_delta_I_traj_interp;
    
    
    this_total_incidence = transpose(params.beta*(this_mu_delta_I_traj.*this_I_traj.*this_mu_epsilon_S_traj.*this_S_traj));
    total_incidence(count,:) = this_total_incidence;
%     this_total_incidence_interp = interp1(this_t_shift_plt,this_total_incidence,this_t_shift_plt_finer,'spline');
%     total_incidence_interp(count,:) = this_total_incidence_interp;
    this_FOS = results.FOS_traj';
    FOS_collect(count,:) = this_FOS;
    
    if count==3
       this_FOS_v = results.FOS_traj'; 
    end
    
    %homogeneous
%     this_I_traj_homogeneous = results.I_traj_v;
%     this_S_traj_homogeneous = results.S_traj_v;
%     total_incidence_homogeneous = params_v.beta*(this_I_traj_homogeneous.*this_S_traj_homogeneous);
%     total_incidence_homogeneous_interp = interp1(this_t_shift_plt,total_incidence_homogeneous,this_t_shift_plt_finer,'spline');
    
    this_marg_eps_S_traj = results.marg_eps_S_traj;
    collect_marginal_susceptibility(count,:,:) = this_marg_eps_S_traj;
        
%     for kk=1:length(params.t_span)
%         this_S_marg(1,:) = sum(results.S_traj_discrete(kk,:,:),3);
%         marg_eps_S_traj(kk,:) = (1/results.S_traj(kk))*this_S_marg;
%     end
    
%     for counter=1:length(marg_eps_S_traj)
%         
%         this_marg_eps_S_traj_interp = interp1(this_t_shift_plt,marg_eps_S_traj(:,counter),this_t_shift_plt_finer,'spline');
%         marg_eps_S_traj_interp(counter,:) = this_marg_eps_S_traj_interp;
%         % time goes down!!
%     end
    
    %     this_marginal_susceptibility = marg_eps_S_traj_interp(this_ind_t_int,:)';
    
    %correlation,time,epsilon
%     marginal_susceptibility_interp(count,:,:) = marg_eps_S_traj_interp;
    
    
    this_marg_delta_I_traj = results.marg_delta_I_traj;
    collect_marginal_transmissibility(count,:,:) = this_marg_delta_I_traj;
%     for kk=1:length(params.t_span)
%         this_I_marg(1,:) = sum(results.I_traj_discrete(kk,:,:),2);
%         marg_delta_I_traj(kk,:) = (1/results.I_traj(kk))*this_I_marg;
%     end
    
%     for counter=1:length(marg_delta_I_traj)
%         
%         this_marg_delta_I_traj_interp = interp1(this_t_shift_plt,marg_delta_I_traj(:,counter),this_t_shift_plt_finer,'spline');
%         marg_delta_I_traj_interp(:,counter) = this_marg_delta_I_traj_interp;
%         % time goes down!!
%     end
    
    %     this_marginal_transmissibility = marg_delta_I_traj_interp(this_ind_t_int,:)';
    %correlation,time,delta
%     marginal_transmissibility_interp(count,:,:) = marg_delta_I_traj_interp;
    
    
    % interpolate trajectories too:
%     this_S_traj_interp = interp1(this_t_shift_plt,this_S_traj,t_shift_plt_finer,'spline');
    %     this_I_traj_interp = interp1(this_t_shift_plt,this_I_traj,t_shift_plt_finer);
    
end


%% Plotting
f1 = figure(1);
set(f1, 'Position', [100 100 1000 1000]);
annotation('rectangle',[0.37 0.11 0.55 0.21], 'Color','black','LineWidth',2);
annotation('rectangle',[0.37 0.32 0.55 0.21], 'Color','black','LineWidth',2);
annotation('rectangle',[0.37 0.53 0.55 0.21], 'Color','black','LineWidth',2);
annotation('rectangle',[0.37 0.74 0.55 0.25], 'Color','black','LineWidth',2);

sp1=subplot(4,3,1);

old_pos = get(sp1, 'Position');
set(sp1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
old_pos = get(sp1, 'Position');


%%
f1=subplot(4,3,4);

% this_p(1) =plot(t_shift_plt_finer(1,:), total_incidence_homogeneous_interp,'Color',colors_rgb(4,:),'LineWidth',2.5); hold on;

box('off');
% this_FOS_v
for count = 1:2
    
    
%     this_p(5-count) =plot(t_shift_plt_finer(4-count,:), total_incidence_interp(4-count,:),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
%     this_p(5-count).Color(4)=0.85;
%     plot(t_shift_plt(count,:), total_incidence(count,:),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
%     this_p(5-count).Color(4)=0.85;
    this_p(count+1) = plot(t_shift_plt(count,:), FOS_collect(count,:),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
    
end

this_p(1) = plot(t_shift_plt(1,:)+22, this_FOS_v,'k--','LineWidth',3); hold on;
% axis([0 this_t_end_plt 0 0.025]);
xlim([0 this_t_end_plt]);
xlabel('Time (days)'); ylabel({'Cumulative';'Infections'});
% ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

xticks([0 50 100 150 200]);
% yticks([0 0.005 0.01 0.015 0.02 0.025]);
% set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'0.005'},{'0.01'},{'0.015'},{'0.02'},{'0.025'}]);
set(f1,'xticklabel',{[]});
box('off');

txt = {'0'};
text(-0.06,0.05,txt,'Units','normalized',...
    'FontSize',14,'FontWeight','normal','FontName', 'Times');

txt = {'A'};
text(0.025,1.045,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');



%% now plot mean susceptibility
f1=subplot(4,3,7);
for count=1:2
    this_q=plot(t_shift_plt(count,:), mu_epsilon_S_traj(count,:),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
%     this_q=plot(t_shift_plt_finer(count,:), mu_epsilon_S_traj_interp(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
%     this_q.Color(4)=0.85;
end

for kk=1:length(var_eps)
    plot(t_shift_plt_finer(count,:),var_eps(kk)*ones(size(t_shift_plt_finer)),'k','LineWidth',1); hold on;
    if kk==1
        text(0.6,0.91,txt_eps1,'Interpreter','Latex','Units','normalized','FontSize',11);
    elseif kk==2
        text(0.6,0.68,txt_eps2,'Interpreter','Latex','Units','normalized','FontSize',11);
    elseif kk==3
        text(0.6,0.50,txt_eps3,'Interpreter','Latex','Units','normalized','FontSize',11);
    else
        text(0.6,0.32,txt_eps4,'Interpreter','Latex','Units','normalized','FontSize',11);
    end
end


axis([0 this_t_end_plt 0 1.1]);
xlabel('Time (days)'); ylabel({'Mean' ; 'Susceptiblity, $\bar{\varepsilon}(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');

yticks([0.1 0.2 0.4 0.6 0.8 1 1.1]);
set(f1,'xticklabel',{[]},'yticklabel',[{''},{'0.2'},{'0.4'},{'0.6'},{'0.8'},{'1.0'},{' '}]);
box('off');


%% plot effective transmission rate
f1=subplot(4,3,10);
for count=1:2
%     this_q=plot(t_shift_plt_finer(4-count,:),beta_mu_delta_I_traj_interp(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    this_q=plot(t_shift_plt(count,:),beta_mu_delta_I_traj(count,:),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
%     this_q=plot(t_shift_plt(count,:),1./disp_delta_I_traj(count,:),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
    this_q.Color(4)=0.85;
end
% xlim([0 this_t_end_plt]);
axis([0 this_t_end_plt 0 0.8]);
xlabel('Time (days)'); ylabel({'Effective'; 'Transmission'; 'Rate, $\beta\,\bar{\delta}_I(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');


% yticks([0.2 0.25 0.3 0.35 0.4]);
% set(f1,'yticklabel',[{'0.20'},{'0.25'},{'0.3'},{'0.35'},{''}]);
box('off');

% txt = {'0.40'};
% text(-0.135,0.935,txt,'Units','normalized',...
%     'FontSize',14,'FontWeight','normal','FontName', 'Times');


%% handle marginals: plot at different epsilon levels
% find mean susceptibility level for each case
for count_eps = 1:length(var_eps)
    
    this_epsilon_level = var_eps(count_eps);
    
    for count=1:2
        
        
        ind = find(mu_epsilon_S_traj(count,:)<this_epsilon_level);
        this_ind_t_int = ind(1);
        
        this_t_int = t_shift_plt(count,this_ind_t_int);
        %         t_int(count_eps,count) = this_t_int;
        
        this_mu_epsilon = mu_epsilon_S_traj(count,this_ind_t_int);
        %         mu_epsilon(count) = this_mu_epsilon;
        collect_mu_eps(count_eps,count) = this_mu_epsilon;
        
        
        
        
        this_mu_delta = mu_delta_I_traj(count,this_ind_t_int);
        %         mu_delta(count) = this_mu_delta;
        collect_mu_delta(count_eps,count) = this_mu_delta;
        
        this_var_delta = var_delta_I_traj(count,this_ind_t_int);
        %         mu_delta(count) = this_mu_delta;
        collect_var_delta(count_eps,count) = this_var_delta;
        
        this_CV = sqrt(this_var_delta)/this_mu_delta;
        collect_CV_delta(count_eps,count) = this_CV;
        collect_disp_delta(count_eps,count) = 1/this_CV^2;
        
        this_p_nb = 1 - this_mu_delta/this_var_delta;
        collect_p_delta(count_eps,count) = this_p_nb;
        
        this_marginal_susceptibility = reshape(collect_marginal_susceptibility(count,this_ind_t_int,:),1,length(params.eps));
        collect_susceptibility_distribution(count_eps,count,:) = this_marginal_susceptibility;
        % results.marg_eps_S_traj(this_ind_t_int,:)*S_traj(count,this_ind_t_int);
        
%         collect_marginal_transmissibility
        
        this_marginal_transmissibility = reshape(collect_marginal_transmissibility(count,this_ind_t_int,:),1,length(params.del));
        collect_transmissibility_distribution(count_eps,count,:) = this_marginal_transmissibility;
%         results.marg_delta_I_traj(this_ind_t_int,:)*I_traj(count,this_ind_t_int);
        
%         collect_marginal_susceptibility(count_eps,count,:) = this_marginal_susceptibility;
%         collect_marginal_transmissibility(count_eps,count,:) = this_marginal_transmissibility;
        
        % for plotting
%         ind = find(params.eps > this_epsilon_level);
%         this_ind_eps_int = ind(1);
        %         ind_eps_int(count) =  this_ind_eps_int;
        
%         ind = find(params.del > this_mu_delta);
%         this_ind_delta_int = ind(1);
        %         ind_delta_int(count) = this_ind_delta_int;
        %         collect_ind_delta_int(count_eps,count) = this_ind_delta_int;
        
        
        
        
    end
    
end


%% now plot marginals for susceptibility
for count_eps = 1:length(var_eps)
    
    this_epsilon_level = var_eps(count_eps);
    
    for count=1:2
        
        plt_this_susceptibility_distribution = reshape(collect_susceptibility_distribution(count_eps,count,:),1,length(params.eps));
        
        subplot(4,3,(2+3*(count_eps-1)));

        
        if count==1
%             plot(params.eps, collect_mu_eps(count_eps,count)*exp(-collect_mu_eps(count_eps,count)*params.eps),'r.-','MarkerSize',20);
            this_q=plot(params.eps,plt_this_susceptibility_distribution,'.-','Color',colors_rgb(count,:),'LineWidth',2.5,'MarkerSize',20); hold on;
%             this_q.Color(4)=0.75;
            plot(this_epsilon_level,0,'.','Color',colors_rgb(count,:),'LineWidth',2,'MarkerSize',20);
        elseif count==2
            this_q=plot(params.eps,plt_this_susceptibility_distribution,'.-','Color',colors_rgb(count,:),'LineWidth',1.5,'MarkerSize',20); hold on;
            this_q.Color(4)=0.75;
            % this_poiss = plot(params.eps, poisspdf(0:(length(params.eps)-1),collect_mu_eps(count_eps,count)),'k.-','LineWidth',0.5,'MarkerSize',12);
            plot(this_epsilon_level,0,'o','Color',colors_rgb(count,:),'LineWidth',2,'MarkerSize',8);
            
%             legend(this_poiss,)

        else
            
            plot(this_epsilon_level,0,'o','Color',colors_rgb(count,:),'LineWidth',2,'MarkerSize',12);
        end
        
        axis([0 5 0 1]);
%         xlim([0 5]);
        yticks([]);
        ylabel({['$g_S(t_',num2str(count_eps-1),',\varepsilon)$']},'interpreter','latex');
        f1=gca;
        f1.LineWidth = 1;
        f1.FontSize = 16;
        f1.FontWeight = 'normal';
        f1.FontName = 'Times';
        
        
        
        
        if count==2
            
            if count_eps == 1
                
                old_pos = get(f1, 'Position');
                set(f1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
                old_pos = get(f1, 'Position');
                
            else
                
                set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
                old_pos = get(f1, 'Position');
                
            end
            
            set(f1,'xticklabel',{[]},'yticklabel',{[]});
            box('off');
            
%             legend(this_poiss,'Poisson','Position',[1.2 0.82 0.12 0.1],'FontSize',14);
%             legend boxoff
            
        end
        
        
        if count_eps == 1
            
            txt = {'B'};
            text(0.025,1.045,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
            title('Susceptibility','FontSize',16,'FontWeight','normal');
            
            text(0.6,0.85,txt_eps1,'Interpreter','Latex','Units','normalized','FontSize',14);
            
            
        elseif count_eps == 2
            
            text(0.6,0.85,txt_eps2,'Interpreter','Latex','Units','normalized','FontSize',14);
            
            
        elseif count_eps == 3
            
            text(0.6,0.85,txt_eps3,'Interpreter','Latex','Units','normalized','FontSize',14);
            
            
        else
            
            text(0.6,0.85,txt_eps4,'Interpreter','Latex','Units','normalized','FontSize',14);
            
            xticks([0 1 2 3]);
            set(f1,'xticklabel',[{'0'},{'1'},{'2'},{'3'},{'4'},{'5'}]);
            
            xlabel('Susceptibility, $\varepsilon$','interpreter','latex');
        end
        
        
    end
    
end



%% now plot marginals for transmissibility
for count_eps = 1:length(var_eps)
    
    this_epsilon_level = var_eps(count_eps);
    
    subplot(4,3,(3+3*(count_eps-1)));
    for count=1:2
        
%         this_collect_mu_delta = collect_mu_delta(count_eps,4-count);
        this_collect_mu_delta = collect_mu_delta(count_eps,count);
        
        %         collect_marginal_transmissibility(count_eps,count,:)
%         plt_this_transmissibility_distribution(1,:) = collect_marginal_transmissibility(count_eps,4-count,:);
        plt_this_transmissibility_distribution = reshape(collect_transmissibility_distribution(count_eps,count,:),1,length(params.eps));
%         plot(params.eps, collect_mu_eps(count_eps,count)*exp(-collect_mu_eps(count_eps,count)*params.eps),'r.-','MarkerSize',20);
%         plot(params.del, poisspdf(0:(length(params.del)-1),collect_mu_delta(count_eps,count)),'k.-','MarkerSize',20); hold on;
%         this_q=plot(params.eps,plt_this_transmissibility_distribution,'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
        this_q=plot(params.del,plt_this_transmissibility_distribution,'.-','Color',colors_rgb(count,:),'LineWidth',2.5,'MarkerSize',20); hold on;
        this_q.Color(4)=0.8;
%         plot(params.del-1,gampdf(params.del,collect_disp_delta(count_eps,count), 1/collect_mu_delta(count_eps,count)),'k.-','LineWidth',0.5,'MarkerSize',12);
%         plot(params.del,nbinpdf(params.del,collect_disp_delta(count_eps,count),collect_p_delta(count_eps,count)),'k.-','LineWidth',0.5,'MarkerSize',12);
        
        plot(this_collect_mu_delta, 0,'.','Color',colors_rgb(count,:),'LineWidth',2,'MarkerSize',30);

        axis([0 5 0 0.1]);
%         xlim([0 5]);
        yticks([]);
        %         xlabel('Susceptibility, $\varepsilon$','interpreter','latex');
        ylabel({['$h_I(t_',num2str(count_eps-1),',\delta)$']},'interpreter','latex');
        f2=gca;
        f2.LineWidth = 1;
        f2.FontSize = 16;
        f2.FontWeight = 'normal';
        f2.FontName = 'Times';
        
        
        if count==2
            
            if count_eps == 1
                
                old_pos = get(f2, 'Position');
                set(f2,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
                old_pos = get(f2, 'Position');
                
            else
                
                set(f2,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
                old_pos = get(f2, 'Position');
                
            end
            
            set(f2,'xticklabel',{[]},'yticklabel',{[]});
            box('off');
            
        end
        
        
        
    end
    
    if count_eps == 1
        
        txt = {'C'};
        text(0.025,1.045,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        title('Transmissibility','FontSize',16,'FontWeight','normal');
        
%         text(0.6,0.85,txt_eps1,'Interpreter','Latex','Units','normalized','FontSize',14);
        
    elseif count_eps == 2
        
%         text(0.6,0.85,txt_eps2,'Interpreter','Latex','Units','normalized','FontSize',14);
        
        
    elseif count_eps == 3
        
%         text(0.6,0.85,txt_eps3,'Interpreter','Latex','Units','normalized','FontSize',14);
        
        
    else
        
%         text(0.6,0.85,txt_eps4,'Interpreter','Latex','Units','normalized','FontSize',14);
        
        xticks([0 1 2 3 4 5]);
        set(f2,'xticklabel',[{'0'},{'1'},{'2'},{'3'},{'4'},{'5'}]);
        
        xlabel('Transmissibility, $\delta$','interpreter','latex');
        
        
    end
    
end


%%
sp1 = subplot(4,3,1);
delete(sp1);
legend_char1 = 'Classic SIR';
legend_char2 = 'No Correlation, $\rho = 0$';
legend_char3 = 'Positive Correlation, $\rho > 0$';

% legend_char4 = 'Negative Correlation, $\rho < 0$';


legend(this_p,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Position',[0.17 0.82 0.12 0.08],'FontSize',14);

%% save figure ?
if save_fig_ans==1
    
    figures_location = './figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));
    
end