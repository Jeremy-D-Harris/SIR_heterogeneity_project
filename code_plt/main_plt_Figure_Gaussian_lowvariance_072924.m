% function void = main_plt_Figure_Gaussian_lowvariance_072924(void)

% plot results SIR model with transmissibility & susceptibility variation


%% set up
% plot results from Gamma Distribution

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure2_Gaussian_lowvariance_072924';

% medium blue, black, grey, violet, green
colors_rgb = [74,112,188.5; 0 0 0; 166 166 166; 169,90,161;0, 158, 115]/255;

% define ending time:
this_t_end_plt = 200;

epsilon_levels = [1.01, 0.9985];


% load results from file
file_location = '../data/';


% results from simulations
infile_results = 'GaussianNoCorrelation_lowvariance.mat';
load(strcat(file_location,infile_results));

fprintf('Plotting Results of Gaussian Distribution... \n');

%% Plotting
% f1 = figure(1); set(f1, 'Position', [100 500 1000 700]);
f1 = figure(1); set(f1, 'Position', [100 500 1200 700]);

ind_time_pt = 40;

% panel A: incident infections
subplot(2,3,1);
this_p(1) = plot(params.t_span, results_classic.total_incidence','Color',colors_rgb(3,:),'LineWidth',2.5); hold on;
this_p(2) = plot(params.t_span, results.total_incidence','Color',colors_rgb(1,:),'LineWidth',2.5); hold on;
this_p(3) = plot(params.t_span, results_var_susc.total_incidence','--','Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
this_p(4) = plot(params.t_span, results_reduced.total_incidence,':','Color',colors_rgb(5,:),'MarkerSize',6,'LineWidth',2);
for kk = 1:2
    if kk==1

        % plot(0,this_total_incidence_interp(this_ind_t_int(kk)),'o','Color',colors_rgb(2,:),'MarkerSize',8, 'LineWidth',2); hold on;
        plot(0,results.total_incidence(1),'o','Color',colors_rgb(2,:),'MarkerSize',8, 'LineWidth',2); hold on;
        text(0.01,0.074,'$t_0$','Interpreter','Latex','Units','normalized','FontSize',16)
    else


        plot(params.t_span(ind_time_pt),results.total_incidence(ind_time_pt),'o','Color',colors_rgb(4,:),'MarkerSize',8, 'LineWidth',2); hold on;
        text(0.21,0.1,'$t_1$','Interpreter','Latex','Units','normalized','FontSize',16)
    end

end


axis([0 this_t_end_plt 0 0.025]);
xlabel('Time (days)'); ylabel({'Incident infections, $\eta(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'A'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

legend_char1 = 'SIR';
legend_char2 = 'Variation ($\varepsilon$,\,$\delta$)';
legend_char3 = 'Variation ($\varepsilon$)';
legend_char4 = 'Reduced Model';

legend(this_p,{legend_char1,legend_char2, legend_char3, legend_char4}, 'Position', [0.245 0.84 0.09 0.07],'FontSize',10,'Interpreter','Latex');



% panel B: CV^2 Susceptibility
subplot(2,3,2);
this_q(2) = plot(params.t_span, results.CV2_eps_S_traj,'Color',colors_rgb(1,:),'LineWidth',2.5);hold on;
this_q(1) = plot(params.t_span, results_var_susc.CV2_eps_S_traj,'--','Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
this_q(3) = plot(params.t_span, params.variance_eps_S./results.mean_eps_S_traj.^2,':','Color',colors_rgb(5,:),'LineWidth',2);hold on;

% text(0.7,0.6,'$\frac{1}{k_{\varepsilon}} = \frac{1}{3}$','Interpreter','Latex','Units','normalized','FontSize',16,'Color',colors_rgb(5,:))
% this_p(1).Color(4) = 0.8;

axis([0 this_t_end_plt 0 0.5]);
xlabel('Time (days)');
ylabel({'Coefficient of Variation (Squared)'});
title('Susceptibility','FontWeight','normal');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'B'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

% legend_char1 = 'SIR';
legend_char1 = 'Variation Susceptibility ($\varepsilon$,\,$\delta$)';
legend_char2 = 'Variation Susceptibility ($\varepsilon$)';
legend_char3 = 'Reduced Model';

legend(this_q,{legend_char1,legend_char2, legend_char3},'Position',[0.502 0.855 0.085 0.06],'FontSize',10,'Interpreter','latex');

%% panel C: CV transmissibility
subplot(2,3,3);
this_r(2) = plot(params.t_span, results.CV2_delta_I_traj,'Color',colors_rgb(1,:),'LineWidth',2.5);hold on;
this_r(1) = plot(params.t_span, results.CV2_delta_S_traj,'--','Color',colors_rgb(3,:),'LineWidth',2.5);hold on;
% this_r(1).Color(4) = 0.85;
% this_r(3) = plot(params.t_span, 1/params.shape_delta*ones(size(params.t_span)),':','Color',colors_rgb(5,:),'LineWidth',2);hold on;
this_r(3) = plot(params.t_span, params.variance_delta_S/params.mean_delta_S^2*ones(size(params.t_span)),':','Color',colors_rgb(5,:),'LineWidth',2);hold on;

% text(0.7,0.125,'$\frac{1}{k_\delta} = \frac{1}{10}$','Interpreter','Latex','Units','normalized','FontSize',16,'Color',colors_rgb(5,:))

axis([0 this_t_end_plt 0 0.1]);
xlabel('Time (days)');
ylabel({'Coefficient of Variation (Squared)'});
title('Transmissibility','FontWeight','normal');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'C'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


legend_char1 = 'Variation Potential Transmissibility $(\varepsilon,\,\delta)$ ';
legend_char2 = 'Variation Effective Transmissibility $(\varepsilon,\,\delta)$, ';
legend_char3 = 'Reduced Model';

legend(this_r,{legend_char1,legend_char2,legend_char3}, 'interpreter','latex','Position',[0.755 0.855 0.095 0.06],'FontSize',10);


% panel D: Initial joint distribution
subplot(2,3,4);
imagesc(params.eps,params.del,results.init_joint_S);
set(gca,'YDir','normal');
colormap(parula);
f1=gca;
xticks([0 1 2 3]);
yticks([0 1 2 3]);
axis([0 3 0 3]);

xlabel('Susceptibility, $\varepsilon$','interpreter','latex');
ylabel('Potential Transmissibility, $\delta$','interpreter','latex');
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

title('Initial Joint Distribution, $f_S(0,\varepsilon,\delta)$','interpreter','latex','FontSize',12)

txt = {'D'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


% panel E: marginals for susceptibility
subplot(2,3,5);

this_s(1) = plot(params.eps,results.marginal_eps_S_traj(1,:),'Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
this_s(2) = plot(params.eps,results.marginal_eps_S_traj(ind_time_pt,:),'--','Color',colors_rgb(4,:),'LineWidth',2.5); hold on;

plot(results.mean_eps_S_traj(1),0,'o','Color',colors_rgb(2,:),'LineWidth',2,'MarkerSize',14);
plot(results.mean_eps_S_traj(ind_time_pt),0,'o','Color',colors_rgb(4,:),'LineWidth',2,'MarkerSize',9);


text(0.33,0.07,'$\overline{\varepsilon}$','Interpreter','Latex','Units','normalized','FontSize',16)

axis([0 3 0 1.5]);
yticks([]);
xlabel('Susceptibility, $\varepsilon$','interpreter','latex');
ylabel({'Population Density'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

title('Susceptibility','FontWeight','normal');

txt = {'E'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

legend_char1 = 'Marginal Susceptibility at $t_0$, $g_S(t_0,\varepsilon)$ ';
legend_char2 = 'Marginal Susceptibility at $t_1$, $g_S(t_1,\varepsilon)$ ';

legend(this_s,{legend_char1,legend_char2}, 'Interpreter','Latex','Position',[0.488 0.393 0.09 0.05],'FontSize',10);

% panel F: marginals for transmissibility
subplot(2,3,6);
this_t(1) = plot(params.eps,results.marginal_delta_S_traj(1,:),'Color',colors_rgb(2,:),'LineWidth', 2.5); hold on;
this_t(2) = plot(params.eps,results.marginal_delta_I_traj(ind_time_pt,:),'--','Color',colors_rgb(4,:),'LineWidth',2.5); hold on;

plot(results.mean_delta_S_traj(1),0,'o','Color',colors_rgb(2,:),'LineWidth',2,'MarkerSize',12);
plot(results.mean_delta_I_traj(ind_time_pt),0,'o','Color',colors_rgb(4,:),'LineWidth',2,'MarkerSize',9);

text(0.33,0.07,'$\overline{\delta}$','Interpreter','Latex','Units','normalized','FontSize',16)

axis([0 3 0 2.4]);

yticks([]);
xlabel('Transmissibility, $\delta$','interpreter','latex');
ylabel({'Population Density'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

title('Transmissibility','FontWeight','normal');

txt = {'F'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

legend_char1 = 'Marginal Potential Transmissibility at $t_0$, $h_S(t_0,\delta)$';
legend_char2 = 'Marginal Effective Transmissibility at $t_1$, $h_I(t_1,\delta)$';

legend(this_t,{legend_char1,legend_char2}, 'Interpreter','Latex','Position',[0.76 0.393 0.06 0.05],'FontSize',10);


%% save figure?
if save_fig_ans==1

    figures_location = './../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n');
    fprintf(strcat(figures_location,'\n\n'));

end