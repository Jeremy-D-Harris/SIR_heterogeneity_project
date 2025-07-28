% want to plot correlations vs speed & strength

%%

clear all; close all; clc;

%%
save_fig_ans = 1;
% save figure:
% 0 = no, 1 = yes

%%

figure_name = 'Figure4_correlation_speedstrength_061825';

%create color gradiets
c2 = [133,192,249]/255; % light blue
c1 = [15,32,128]/255; % dark blue

depth = 9;
grad1=colorGradient(c2,c1,depth);

colors_rgb = grad1;
grey = [166 166 166]/255;

gamma = 0.1; %recovery rate (per day)

%% load results from file
file_location = '../data/';


%% (1) positive
this_infile = 'GaussianPositiveCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));

corr_coef(1) = params.corr_coeff;
R0_collect(1) = results.Rt_traj(1);
total_incidence_collect(:,1) = results.total_incidence;
total_incidence_collect_classic(:,1) = results_classic.total_incidence;


%% run corresponding classic SIR
% to match R0

% define sim params
params.bet = R0_collect(1)*params.gam;

% run model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%find init S,I,R from classic
init_conds_SIR_classic = [results_classic.S_traj(1);results_classic.I_traj(1);results_classic.R_traj(1)];

[t,y_traj_classic] = ode45(@(t,y)simulate_SIR(t,y,params), params.t_span, init_conds_SIR_classic, options);

S_traj_SIR_classic = y_traj_classic(:,1);
I_traj_SIR_classic = y_traj_classic(:,2);
% R_traj_SIR_classic = y_traj_classic(:,3);

% total incidence
total_incidence_classic = params.bet*I_traj_SIR_classic.*S_traj_SIR_classic;

%store
total_incidence_collect_r0(:,1) = total_incidence_classic;
%%


%% (2) positive
this_infile = 'GaussianPositiveCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

corr_coef(2) = params.corr_coeff;
R0_collect(2) = results.Rt_traj(1);
total_incidence_collect(:,2) = results.total_incidence;
total_incidence_collect_classic(:,2) = results_classic.total_incidence;

%% run corresponding classic SIR
% to match R0

% define sim params
params.bet = R0_collect(2)*params.gam;

% run model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%find init S,I,R from classic
init_conds_SIR_classic = [results_classic.S_traj(1);results_classic.I_traj(1);results_classic.R_traj(1)];

[t,y_traj_classic] = ode45(@(t,y)simulate_SIR(t,y,params), params.t_span, init_conds_SIR_classic, options);

S_traj_SIR_classic = y_traj_classic(:,1);
I_traj_SIR_classic = y_traj_classic(:,2);
% R_traj_SIR_classic = y_traj_classic(:,3);

% total incidence
total_incidence_classic = params.bet*I_traj_SIR_classic.*S_traj_SIR_classic;

%store
total_incidence_collect_r0(:,2) = total_incidence_classic;
%%


%% (3) independent
this_infile = 'GaussianNoCorrelation.mat';
load(strcat(file_location,this_infile));

% should update?
corr_coef(3) = params.corr_coeff;
R0_collect(3) = results.Rt_traj(1);
total_incidence_collect(:,3) = results.total_incidence;
total_incidence_collect_classic(:,3) = results_classic.total_incidence;

%% run corresponding classic SIR
% to match R0

% define sim params
params.bet = R0_collect(3)*params.gam;

% run model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%find init S,I,R from classic
init_conds_SIR_classic = [results_classic.S_traj(1);results_classic.I_traj(1);results_classic.R_traj(1)];

[t,y_traj_classic] = ode45(@(t,y)simulate_SIR(t,y,params), params.t_span, init_conds_SIR_classic, options);

S_traj_SIR_classic = y_traj_classic(:,1);
I_traj_SIR_classic = y_traj_classic(:,2);
% R_traj_SIR_classic = y_traj_classic(:,3);

% total incidence
total_incidence_classic = params.bet*I_traj_SIR_classic.*S_traj_SIR_classic;

%store
total_incidence_collect_r0(:,3) = total_incidence_classic;
%%


%% (4) negative
this_infile = 'GaussianNegativeCorrelation_0pt3.mat';
load(strcat(file_location,this_infile));

corr_coef(4) = params.corr_coeff;
R0_collect(4) = results.Rt_traj(1);
total_incidence_collect(:,4) = results.total_incidence;
total_incidence_collect_classic(:,4) = results_classic.total_incidence;


%% run corresponding classic SIR
% to match R0

% define sim params
params.bet = R0_collect(4)*params.gam;

% run model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%find init S,I,R from classic
init_conds_SIR_classic = [results_classic.S_traj(1);results_classic.I_traj(1);results_classic.R_traj(1)];

[t,y_traj_classic] = ode45(@(t,y)simulate_SIR(t,y,params), params.t_span, init_conds_SIR_classic, options);

S_traj_SIR_classic = y_traj_classic(:,1);
I_traj_SIR_classic = y_traj_classic(:,2);
% R_traj_SIR_classic = y_traj_classic(:,3);

% total incidence
total_incidence_classic = params.bet*I_traj_SIR_classic.*S_traj_SIR_classic;

%store
total_incidence_collect_r0(:,4) = total_incidence_classic;
%%

%% (5) negative
this_infile = 'GaussianNegativeCorrelation_0pt6.mat';
load(strcat(file_location,this_infile));

corr_coef(5) = params.corr_coeff;
R0_collect(5) = results.Rt_traj(1);
total_incidence_collect(:,5) = results.total_incidence;
total_incidence_collect_classic(:,5) = results_classic.total_incidence;

%% run corresponding classic SIR
% to match R0

% define sim params
params.bet = R0_collect(5)*params.gam;

% run model
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

%find init S,I,R from classic
init_conds_SIR_classic = [results_classic.S_traj(1);results_classic.I_traj(1);results_classic.R_traj(1)];

[t,y_traj_classic] = ode45(@(t,y)simulate_SIR(t,y,params), params.t_span, init_conds_SIR_classic, options);

S_traj_SIR_classic = y_traj_classic(:,1);
I_traj_SIR_classic = y_traj_classic(:,2);
% R_traj_SIR_classic = y_traj_classic(:,3);

% total incidence
total_incidence_classic = params.bet*I_traj_SIR_classic.*S_traj_SIR_classic;

%store
total_incidence_collect_r0(:,5) = total_incidence_classic;
%%


%% calculate analytic results for panel B%

% analytic R0
corrcoef_vary = -1:0.05:1;

variance_eps_S = 0.49;
variance_delta_S = 0.34;
params.bet = 0.2;
R0_analytic = (params.bet/params.gam)*(1+corrcoef_vary*sqrt(variance_eps_S)*sqrt(variance_delta_S)); %remember <epsilon(0)> = <delta_S(0)> = 1

%% Plotting
X = get(0,'ScreenPixelsPerInch'); %determine screen pixels per inch (96 on windows, 72 on mac os)
 factor = X/72;
 fig1 = figure(1); set(fig1, 'Position', [100 500 factor*900 factor*350]);
 fig1.OuterPosition(4) = 650;

%%%%%%%%%Fig4a%%%%%%%%%%%%
LTD = 1.5; %line thickness SIR dots

ax1 = subplot(5,2,1);
count = 1;
sirbg = plot(params.t_span,total_incidence_collect_classic(:,count),'Color',grey,'LineWidth',2.5);hold on;
sirr0 = plot(params.t_span,total_incidence_collect_r0(:,count),':','Color','black','LineWidth',LTD);
varepsdel = plot(params.t_span,total_incidence_collect(:,count),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5);

dummyg = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');

axis([0 300 0 0.03]);
xticklabels([]);
lgd=legend([varepsdel,sirbg,sirr0],{'Variation ($\varepsilon,\delta$):','SIR with matching $\beta, \gamma \ $ ','SIR with matching $\mathcal R_0$'},'Interpreter','Latex','Location','northeast','Box','off','FontSize',12);

ah1 = axes('position',get(gca,'position'),'visible','off');
lgd2 = legend(ah1, dummyg, '$\rho = 0.6$','Interpreter','Latex','Location','northeast','Box','off','FontSize',14);
lgd2.TextColor = colors_rgb(2*count-1,:);

f1=ah1;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

f1=ax1;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'A'};
text(0.025,1.2,txt,'Units','normalized','FontSize',18,'FontWeight','bold');
box on

ax2 = subplot(5,2,3);
count = 2;
plot(params.t_span,total_incidence_collect_classic(:,count),'Color',grey,'LineWidth',2.5);hold on;
plot(params.t_span,total_incidence_collect_r0(:,count),':','Color','black','LineWidth',LTD);
plot(params.t_span,total_incidence_collect(:,count),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5);

axis([0 300 0 0.03]);
xticklabels([])
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
lgd=legend(dummyh,{'$\rho = 0.3$'},'Interpreter','Latex','Location','northeast','Box','off','FontSize',14);
lgd.TextColor = colors_rgb(2*count-1,:);

f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

ax3 = subplot(5,2,5);
count = 3;
plot(params.t_span,total_incidence_collect_classic(:,count),'Color',grey,'LineWidth',2.5);hold on;
plot(params.t_span,total_incidence_collect_r0(:,count),':','Color','black','LineWidth',LTD); 
plot(params.t_span,total_incidence_collect(:,count),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5);

axis([0 300 0 0.03]);
xticklabels([]);
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
lgd=legend(dummyh,{'$\rho = 0.0$'},'Interpreter','Latex','Location','northeast','Box','off','FontSize',14);
lgd.TextColor = colors_rgb(2*count-1,:);
ylabel({'Incident Infections, $\eta(t)$'},'Interpreter','Latex');

f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

ax4 = subplot(5,2,7);
count = 4;
plot(params.t_span,total_incidence_collect_classic(:,count),'Color',grey,'LineWidth',2.5);hold on;
plot(params.t_span,total_incidence_collect_r0(:,count),':','Color','black','LineWidth',LTD);
plot(params.t_span,total_incidence_collect(:,count),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5);

axis([0 300 0 0.03]);
xticklabels([])
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
lgd=legend(dummyh,{'$\rho = -0.3$'},'Interpreter','Latex','Location','northeast','Box','off','FontSize',14);
lgd.TextColor = colors_rgb(2*count-1,:);

f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

ax5 = subplot(5,2,9);
count = 5;
plot(params.t_span,total_incidence_collect_classic(:,count),'Color',grey,'LineWidth',2.5); hold on;
plot(params.t_span,total_incidence_collect_r0(:,count),':','Color','black','LineWidth',LTD);
plot(params.t_span,total_incidence_collect(:,count),'Color',colors_rgb(2*count-1,:),'LineWidth',2.5);

axis([0 300 0 0.03]);
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
lgd=legend(dummyh,{'$\rho = -0.6$'},'Interpreter','Latex','Location','northeast','Box','off','FontSize',14);
lgd.TextColor = colors_rgb(2*count-1,:);
xlabel('Time (days)');

f1=gca;
f1.LineWidth = 1;
f1.FontSize = 16;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

%tile the above subplots

%x,y,width,height
top = ax1.Position(2)+ax1.Position(4);
bottom = ax5.Position(2);
totalheightpersubplot = (top-bottom)/5;

ax5.Position(4) = totalheightpersubplot;
ax4.Position(4) = totalheightpersubplot;
ax4.Position(2) = ax5.Position(2)+ax5.Position(4);
ax3.Position(4) = totalheightpersubplot;
ax3.Position(2) = ax4.Position(2)+ax4.Position(4);
ax2.Position(4) = totalheightpersubplot;
ax2.Position(2) = ax3.Position(2)+ax3.Position(4);
ax1.Position(4) = totalheightpersubplot;
ax1.Position(2) = ax2.Position(2)+ax2.Position(4);
linkaxes([ax1 ax2 ax3 ax4 ax5],'x')

%%%%%%%%%Fig4b%%%%%%%%%%%

subplot(5,2,[2,4,6,8,10])

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

%save the figure

if save_fig_ans==1

    figures_location = '../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));

end