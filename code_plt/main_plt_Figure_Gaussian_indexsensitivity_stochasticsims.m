% function void = main_plt_Figure_Gaussian_072924(void)

% plot results SIR model with transmissibility & susceptibility variation


%% set up
% plot results from Gamma Distribution

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure5_Gaussian_indexsensitivity_stochastic_080124';

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
f1 = figure(1); set(f1, 'Position', [100 500 1000 850]);

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
p(1) = plot(params.t_span, total_incidence(:,1),'Color',colors_rgb(3,:),'LineWidth',2.5); hold on;

p(2) = plot(params.t_span, total_incidence(:,2),'Color',colors_rgb(1,:),'LineWidth',2.5); hold on;

p(3) = plot(params.t_span, total_incidence(:,3),'Color',colors_rgb(2,:),'LineWidth',2.5); hold on;

p(4) = plot(params.t_span, total_incidence(:,4),'Color',colors_rgb(4,:),'LineWidth',2.5); hold on;



axis([0 this_t_end_plt 0 0.025]);
xlabel('Time (days)'); ylabel({'Incident infections, $\eta(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'B'};
text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

legend_char1 = 'SIR';
legend_char2 = 'Eigendistribution';
legend_char3 = 'Index Case: $\varepsilon = 1$, $\delta = 2$';
legend_char4 = 'Index Case: $\varepsilon = 1$, $\delta = 0.5$';

legend(p, {legend_char1,legend_char2,legend_char3,legend_char4},'Location','northeast','Interpreter','latex');

%% stochastic stuff

%create color gradiets
c2 = [133,192,249]/255; % light blue
c1 = [15,32,128]/255; % dark blue

depth = 9;
[grad1,im]=colorGradient(c2,c1,depth);
% [grad2,im]=colorGradient(c2,c3,depth);

colors_rgb = grad1;

% %% load results from file
% file_location = '../data/';

load(append(file_location,"stochastic_output_31Jul2024125047.mat"))

%% Plotting
% f1 = figure(1); set(f1, 'Position', [100 500 900 350]);

Outbreak_thresh = 50; %50 cases denotes an outbreak

%% plot outbreak proportion, (more than threshold infections) 
subplot(2,2,3);

Traj = length(stochSIR_sir.Final_Size);
labels = ["SIR" "\rho=-0.6" "\rho=0" "\rho=0.6"];
outbreakproportion = [sum(stochSIR_sir.Final_Size>Outbreak_thresh) sum(stochSIR_negcorr.Final_Size>Outbreak_thresh) sum(stochSIR_uncorr.Final_Size>Outbreak_thresh) sum(stochSIR_poscorr.Final_Size>Outbreak_thresh)]/Traj;
b = bar(labels,outbreakproportion);%,'Interpreter','latex')
ylabel("Outbreak probability (>50 cases)")
b.FaceColor = 'flat';
b.CData(1,:) = [166 166 166]/255;
b.CData(2,:) =  colors_rgb(9,:); 
b.CData(3,:) =  colors_rgb(5,:);
b.CData(4,:) =  colors_rgb(1,:);
txt = {'C'};
text(0.025,1.035,txt,'Units','normalized','FontSize',18,'FontWeight','bold');

%calc 95% confidence intervals
pest = outbreakproportion;
C = 0.95;
zscore = C+(1-C)/2;
margin = zscore*( sqrt( (pest.*(1-pest))/Traj ) );

upper = pest+margin;
lower = pest-margin;

%% plot outbreak size, (when more than threshold infections)

colorsbox = [[166 166 166]/255; colors_rgb(9,:); colors_rgb(5,:); colors_rgb(1,:)];
outbreaklist = [stochSIR_sir.Final_Size; stochSIR_negcorr.Final_Size; stochSIR_uncorr.Final_Size; stochSIR_poscorr.Final_Size];
outbreaklist(outbreaklist<=50) = NaN;
outbreaklistproportion = outbreaklist/10000;
ax = subplot(2,2,4);
x = 1:4;
offset = 0.2;
hold(ax);
for ii=1:4
    b2 = boxchart(x(ii)*ones(size(outbreaklistproportion(ii,:)))-offset, outbreaklistproportion(ii,:),'BoxFaceAlpha',1, 'BoxFaceColor', colorsbox(ii,:),'Notch','off');hold on;
    b2.BoxMedianLineColor = [1 1 1];
    b2.JitterOutliers = 'on';
    b2.MarkerStyle = '.';
    b2.MarkerColor = colorsbox(ii,:);
end

%.need to mark on where the deterministic models fall.
load(append(file_location,"GaussianNoCorrelation.mat"),"results_classic")
sirRfinal = results_classic.R_traj(end);
load(append(file_location,"GaussianNoCorrelation.mat"),"results")
UncorrRfinal = results.R_traj(end);
load(append(file_location,"GaussianPositiveCorrelation_0pt6.mat"),"results")
PoscorrRfinal = results.R_traj(end);
load(append(file_location,"GaussianNegativeCorrelation_0pt6.mat"),"results")
NegcorrRfinal = results.R_traj(end);

hold on
plot(1+offset,0.8,"*",'Color',colorsbox(1,:))
plot(2+offset,NegcorrRfinal,"*",'Color',colorsbox(2,:))
plot(3+offset,UncorrRfinal,"*",'Color',colorsbox(3,:))
plot(4+offset,PoscorrRfinal,"*",'Color',colorsbox(4,:))

box on;
ylabel("Final outbreak proportion")
xticks([1 2 3 4])
xticklabels(labels)

txt = {'D'};
text(0.025,1.035,txt,'Units','normalized','FontSize',18,'FontWeight','bold');


%display statistics
disp("median outbreak proportions (SIR,\rho=-0.6,\rho=0,\rho=0.6)")
median(outbreaklistproportion,2,"omitmissing")
disp("quantile outbreak proportions (SIR,\rho=-0.6,\rho=0,\rho=0.6) at [.025 .25 .50 .75 .975]")
quantile(outbreaklistproportion',[.025 .25 .50 .75 .975]) %quantiles for box and whiskers
disp("iqr outbreak proportions (SIR,\rho=-0.6,\rho=0,\rho=0.6)")
iqr(outbreaklistproportion')  %interquartile range







%% save figure?
if save_fig_ans==1

    figures_location = './../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n');
    fprintf(strcat(figures_location,'\n\n'));

end