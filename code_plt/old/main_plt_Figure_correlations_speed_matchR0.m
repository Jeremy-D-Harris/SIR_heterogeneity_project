% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation

%%

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure_correlations_speed_matchR0_080124';


% light blue, medium blue, dark blue
colors_rgb = [133,192,249; 74,112,188.5; 15,32,128]/255;

% define ending time:
this_t_end_plt = 200;

var_eps = linspace(0.999,0.49,201);

ind_level = 18; %18
epsilon_level = var_eps(ind_level);


%%
% need to load all three cases first:
file_location = '../data/';

for count = 1:3
    
    if count==1
        % (1) positive correlation
        infile_independent = 'GaussianPositiveCorrelation_0pt6_matchR0.mat';
        % infile_independent = 'dataIVGaussianSusVarOnlyInverseHighVar.mat';
        %         infile_independent = 'dataPos0.8Figure3.mat';
        load(strcat(file_location,infile_independent));
        
        
    elseif count==2
        % (2) no correlation
        infile_positivecorrelations = 'GaussianNoCorrelation.mat';
        %         infile_positivecorrelations = 'dataIndFigure3.mat';
        load(strcat(file_location,infile_positivecorrelations));
        
        
    else
        % (3) negative correlation
        infile_negativecorrelations = 'GaussianNegativeCorrelation_0pt6_matchR0.mat';
        %         infile_negativecorrelations = 'dataNeg-0.8Figure3.mat';
        load(strcat(file_location,infile_negativecorrelations));
        
    end
    
    
    
    
end


% for plotting
    ind = find(params.eps > epsilon_level);
    this_ind_eps_int = ind(1);
    ind_eps_int(count) =  this_ind_eps_int;
    
    % ind = find(params.del > mu_delta(count));
    % this_ind_delta_int = ind(1);
    ind_delta_int(count) = this_ind_eps_int;


%% Plotting
f1 = figure(1); set(f1, 'Position', [100 500 1200 700]);

%% panel A
subplot(2,3,1);
for count = 1:4
    
    this_p(count) = plot(params.t_span, total_incidence(:,count),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
    
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