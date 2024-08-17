% function void = main_plt_Figure_Gaussian_072924(void)

% plot results SIR model with transmissibility & susceptibility variation


%% set up
% plot results from Gamma Distribution

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure5_Gaussian_indexsensitivity_stochastic_072924';

% medium blue, black, grey, violet, green
colors_rgb = [74,112,188.5; 0 0 0; 166 166 166; 169,90,161;0, 158, 115]/255;

% define ending time:
this_t_end_plt = 200;

% load results from file
file_location = '../data/';


% results from simulations
% infile_results = 'GaussianNoCorrelation_N10000_delta2.mat';
infile_results = 'GaussianNoCorrelation_N10000.mat';
load(strcat(file_location,infile_results));

total_incidence(:,1) = results_classic.total_incidence/params.N;
total_incidence(:,2) = results.total_incidence/params.N;

eps = params.eps;
del = params.del;

ind_temp = find(eps > 1);
ind_eps = ind_temp(1);

ind_temp = [];
set_mean_delta_I = 2;
ind_temp = find(del > set_mean_delta_I);
ind_delta_I_2 = ind_temp(1);

set_mean_delta_I = 0.5;
% ind_eps = find(eps > 1);
ind_temp = [];
ind_temp = find(del > set_mean_delta_I);
ind_delta_I_0pt5 = ind_temp(1);


fprintf('Plotting Results of Gaussian Distribution... \n');

%% Plotting
% f1 = figure(1); set(f1, 'Position', [100 500 1000 700]);
f1 = figure(1); set(f1, 'Position', [100 500 1200 700]);

% ind_time_pt = 50;



%% panel A: Initial joint distribution
subplot(2,2,1);
imagesc(params.eps,params.del,results.init_joint_S); hold on;
plot(params.eps(ind_eps),params.del(ind_delta_I_2), '.', 'MarkerSize',40, 'Color',colors_rgb(2,:));
plot(params.eps(ind_eps),params.del(ind_delta_I_0pt5), '.', 'MarkerSize',40, 'Color',colors_rgb(4,:)); hold on;


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

txt = {'A'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


%% panel B: incident infections

infile_results = 'GaussianNoCorrelation_N10000_delta2.mat';
load(strcat(file_location,infile_results));

total_incidence(:,3) = results.total_incidence/params.N;

infile_results = 'GaussianNoCorrelation_N10000_delta0pt5.mat';
load(strcat(file_location,infile_results));

total_incidence(:,4) = results.total_incidence/params.N;


subplot(2,2,2);
% this_p(1) = plot(params.t_span, results_classic.total_incidence/params.N','Color',colors_rgb(3,:),'LineWidth',2.5); hold on;
this_p(1) = plot(params.t_span, total_incidence(:,1),'Color',colors_rgb(3,:),'LineWidth',2.5); hold on;


this_p(2) = plot(params.t_span, total_incidence(:,2),'Color',colors_rgb(1,:),'LineWidth',2.5); hold on;

this_p(4) = plot(params.t_span, total_incidence(:,4),'Color',colors_rgb(4,:),'LineWidth',2.5); hold on;

this_p(3) = plot(params.t_span, total_incidence(:,3),'Color',colors_rgb(2,:),'LineWidth',2.5); hold on;

% this_p(3) = plot(params.t_span, results_var_susc.total_incidence/params.N','--','Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
% this_p(4) = plot(params.t_span, results_reduced.total_incidence/params.N,':','Color',colors_rgb(5,:),'MarkerSize',6,'LineWidth',2);
% for kk = 1:2
%     if kk==1
% 
%         % plot(0,this_total_incidence_interp(this_ind_t_int(kk)),'o','Color',colors_rgb(2,:),'MarkerSize',8, 'LineWidth',2); hold on;
%         plot(0,results.total_incidence(1),'o','Color',colors_rgb(2,:),'MarkerSize',8, 'LineWidth',2); hold on;
%         text(0.01,0.074,'$t_0$','Interpreter','Latex','Units','normalized','FontSize',16)
%     else
% 
% 
%         % plot(params.t_span(ind_time_pt),results.total_incidence(ind_time_pt),'o','Color',colors_rgb(4,:),'MarkerSize',8, 'LineWidth',2); hold on;
%         % text(0.21,0.1,'$t_1$','Interpreter','Latex','Units','normalized','FontSize',16)
%     end
% 
% end


axis([0 this_t_end_plt 0 0.025]);
xlabel('Time (days)'); ylabel({'Incident infections, $\eta(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'B'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

legend_char1 = 'SIR Model';
legend_char2 = 'Variation ($\varepsilon$,\,$\delta$)';
legend_char3 = 'Variation ($\varepsilon$)';
legend_char4 = 'Reduced Model';

% legend(this_p,{legend_char1,legend_char2, legend_char3, legend_char4}, 'Position', [0.245 0.84 0.09 0.07],'FontSize',10,'Interpreter','Latex');






%% save figure?
if save_fig_ans==1

    figures_location = './../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n');
    fprintf(strcat(figures_location,'\n\n'));

end