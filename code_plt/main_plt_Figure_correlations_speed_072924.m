% function void = main_plt_Figure_correlations_speed(void)
% want to plot correlations vs speed & strength

%%

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure4_correlation_speedstrength_072924';

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
% load results from file
file_location = '../data/';


% results from simulations



%% (1) positive
this_infile = 'GaussianPositiveCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));
% infile_results = 'Gaussian_lowvariance_update072924.mat';
% load(strcat(file_location,infile_results));

R0_collect(1) = results.Rt_traj(1);
total_incidence_collect(:,1) = results.total_incidence;



%% (2) positive
this_infile = 'GaussianPositiveCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

R0_collect(2) = results.Rt_traj(1);
total_incidence_collect(:,2) = results.total_incidence;



% (3) independent
this_infile = 'Gaussian_update072924.mat';
load(strcat(file_location,this_infile));

R0_collect(3) = results.Rt_traj(1);
total_incidence_collect(:,3) = results.total_incidence;



% (4) negative
this_infile = 'GaussianNegativeCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

R0_collect(4) = results.Rt_traj(1);
total_incidence_collect(:,4) = results.total_incidence;



% (5) negative
this_infile = 'GaussianNegativeCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));

R0_collect(5) = results.Rt_traj(1);
total_incidence_collect(:,5) = results.total_incidence;




%% analytic R0
corrcoef_vary = -1:0.05:1;

variance_eps_S = 0.49;
variance_delta_S = 0.34;

stan_dev_eps_S = sqrt(variance_eps_S);
stan_dev_delta_S = sqrt(variance_delta_S);

R0List_analytic = (params.bet/params.gam)*(1+corrcoef_vary*stan_dev_eps_S*stan_dev_delta_S);
% R0List_analytic = R0List_analytic;



%% Plotting
f1 = figure(1); set(f1, 'Position', [100 500 900 350]);


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



if save_fig_ans==1

    figures_location = './figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));

end