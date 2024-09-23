% want to plot correlations vs speed & strength

%%

clear all; close all; clc;

save_fig_ans = 1;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure4_correlation_speedstrength_072924';

%create color gradiets
c2 = [133,192,249]/255; % light blue
c1 = [15,32,128]/255; % dark blue

depth = 9;
grad1=colorGradient(c2,c1,depth);

colors_rgb = grad1;


%% load results from file
file_location = '../data/';


%% (1) positive
this_infile = 'GaussianPositiveCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));

corr_coef(1) = params.corr_coeff;
R0_collect(1) = results.Rt_traj(1);
total_incidence_collect(:,1) = results.total_incidence;


%% (2) positive
this_infile = 'GaussianPositiveCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

corr_coef(2) = params.corr_coeff;
R0_collect(2) = results.Rt_traj(1);
total_incidence_collect(:,2) = results.total_incidence;


%% (3) independent
this_infile = 'GaussianNoCorrelation.mat';
load(strcat(file_location,this_infile));

% should update?
corr_coef(3) = 0; %params.corr_coeff;
R0_collect(3) = results.Rt_traj(1);
total_incidence_collect(:,3) = results.total_incidence;


%% (4) negative
this_infile = 'GaussianNegativeCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

corr_coef(4) = params.corr_coeff;
R0_collect(4) = results.Rt_traj(1);
total_incidence_collect(:,4) = results.total_incidence;


%% (5) negative
this_infile = 'GaussianNegativeCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));

corr_coef(5) = params.corr_coeff;
R0_collect(5) = results.Rt_traj(1);
total_incidence_collect(:,5) = results.total_incidence;


% analytic R0
corrcoef_vary = -1:0.05:1;

variance_eps_S = 0.49;
variance_delta_S = 0.34;

R0_analytic = (params.bet/params.gam)*(1+corrcoef_vary*sqrt(variance_eps_S)*sqrt(variance_delta_S));




%% Plotting
X = get(0,'ScreenPixelsPerInch'); %determine screen pixels per inch (96 on windows, 72 on mac os)
factor = X/72;
f1 = figure(1); set(f1, 'Position', [100 500 factor*900 factor*350]);


%% plot incidence on linear scale
subplot(1,2,1);

for count = 1:5
    this_h(count) =plot(params.t_span,total_incidence_collect(:,count),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5); hold on;
end

axis([0 300 0 0.025]);
xlabel('Time (days)'); ylabel({'Incident Infections, $\eta(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'A'};
text(0.025,1.035,txt,'Units','normalized','FontSize',18,'FontWeight','bold');
box on

legend_char1 = strcat('$\rho = \,$',num2str(corr_coef(1),'%0.1f'));
legend_char2 = strcat('$\rho = \,$',num2str(corr_coef(2),'%0.1f'));
legend_char3 = strcat('$\rho = \,$',num2str(corr_coef(3),'%0.1f'));
legend_char4 = strcat('$\rho = \,$',num2str(corr_coef(4),'%0.1f'));
legend_char5 = strcat('$\rho = \,$',num2str(corr_coef(5),'%0.1f'));

legend(this_h,{legend_char1,legend_char2,legend_char3,legend_char4,legend_char5}, 'Interpreter','Latex','FontSize',16, 'Location','Northeast');


%% Panel B
subplot(1,2,2);

for count = 1:5

    if count == 1
        plot(corrcoef_vary, R0_analytic,'k--','LineWidth', 2.5); hold on;
        plot(corr_coef(count), R0_collect(count),'.', 'Color', colors_rgb(2*count-1,:), 'MarkerSize', 30);hold on;
    else
        plot(corr_coef(count), R0_collect(count),'.', 'Color', colors_rgb(2*count-1,:), 'MarkerSize', 30);hold on;
    end

end

axis([-0.65 0.65 1.5 2.5]);
xticks([-0.6 -0.3 0 0.3 0.6]);
xlabel({'Correlation Coefficient, $\rho$'},'Interpreter','Latex');
ylabel({'Basic Reproduction Number, $\mathcal R_0$'},'Interpreter','Latex');

f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'B'};
text(0.025,1.035,txt,'Units','normalized','FontSize',18,'FontWeight','bold');



if save_fig_ans==1

    figures_location = '../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));

end