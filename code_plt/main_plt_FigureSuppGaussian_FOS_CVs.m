% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation

%%

clear all; close all; clc;

save_fig_ans = 1;
% save figure:
% 0 = no, 1 = yes

figure_name = 'SuppFigure_Gaussian_FOS_CVs_091624';

% light blue, medium blue, dark blue, gray, black
colors_rgb = [133,192,249; 74,112,188.5; 15,32,128; 166 166 166 ; 0,0,0]/255;

% define ending time:  
this_t_end_plt = 200;


%%
% need to load all three cases first:
file_location = '../data/';

for count = 1:3

    if count==1
        % (1) positive correlations
        infile_independent = 'GaussianPositiveCorrelation_0pt6_matchR0.mat';
        load(strcat(file_location,infile_independent));


    elseif count==2
        % (2) no correlations
        infile_positivecorrelations = 'GaussianNoCorrelation.mat';
        load(strcat(file_location,infile_positivecorrelations));


    else
        % (3) negative correlations
        infile_negativecorrelations = 'GaussianNegativeCorrelation_0pt6_matchR0.mat';
        load(strcat(file_location,infile_negativecorrelations));

    end

    %% load parameters
    t_span = params.t_span;
    dt = t_span(2) - t_span(1); % should = 1 uniformly
    bet = params.bet;
    eps = params.eps;
    del = params.del;

    %% load trajectories
    this_incidence = results.total_incidence;
    incidence(:,count) = this_incidence;

    this_cum_infections = cumsum(this_incidence)*dt;
    cum_infections(:,count) = this_cum_infections;

    if count == 2
        cum_infections_classic = cumsum(results_classic.total_incidence)*dt;


        cum_infections_var_susc = cumsum(results_var_susc.total_incidence)*dt;
        CV2_eps_S_traj_var_susc = results_var_susc.CV2_eps_S_traj;


    end

    % CV2 susceptibility over time
    this_CV2_eps_S_traj= results.CV2_eps_S_traj;
    CV2_eps_S_traj(:,count) = this_CV2_eps_S_traj;

    % CV2 transmissibility over time
    this_CV2_delta_I_traj = results.CV2_delta_I_traj;
    CV2_delta_I_traj(:,count) = this_CV2_delta_I_traj;


end




%% Plotting
X = get(0,'ScreenPixelsPerInch'); %determine screen pixels per inch (96 on windows, 72 on mac os)
factor = X/72;
f1 = figure(1); set(f1, 'Position', [100 500 factor*1200 factor*350]);

%%
for count = 1:3


    %% FOS
    subplot(1,3,1);
    %this_h(count+2) = plot(t_span, cum_infections(:,count),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
    this_h(count+1) = plot(t_span, cum_infections(:,count),'Color',colors_rgb(count,:),'LineWidth',2.5); hold on;
    axis([0 this_t_end_plt 0 1]);
    xlabel('Time (days)'); %ylabel('Cumulative Infections');
    ylabel({'Cumulative Infections, $\int_0^t \eta(s) \, ds$ '},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';


    if count==3


        txt = {'A'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

        this_h(1) = plot(t_span, cum_infections_classic,'Color',colors_rgb(4,:),'LineWidth',2.5); hold on;
        %this_h(2) =plot(t_span, cum_infections_var_susc,'--','Color',colors_rgb(5,:),'LineWidth',2.5); hold on;
        yticks(0:0.2:1);
        set(f1,'yticklabel',[{'0'},{'0.2'},{'0.4'},{'0.6'},{'0.8'},{'1.0'}]);
        %         set(f1,'xticklabel',{[]},'yticklabel',[{'0'},{''},{'0.2'},{''},{'0.4'},{''},{'0.6'},{''},{'0.8'},{''},{'1.0'}]);

        legend_char1 = 'SIR';
        %legend_char2 = 'Variation in Susceptibility';
        legend_char5 = 'Negative Correlation, $\rho < 0$';
        legend_char4 = 'No Correlation, $\rho = 0$';
        legend_char3 = 'Positive Correlation, $\rho > 0$';


        %legend(this_h,{legend_char1,legend_char2,legend_char3,legend_char4, legend_char5}, 'Interpreter','Latex','Location',[0.147 0.725 0.1 0.2],'FontSize',10);
        %legend(this_h,{legend_char1,legend_char3,legend_char4, legend_char5}, 'Interpreter','Latex','Location',[0.147 0.725 0.1 0.2],'FontSize',10);
        lgd = legend(this_h,{legend_char1,legend_char3,legend_char4, legend_char5}, 'Interpreter','Latex','Location','NorthWest','FontSize',10);
        %lgd.Position(1) = 0.13;
        %lgd.Position(2) = 0.72;

    end

    %% CV Susceptibility
    subplot(1,3,2);

    this_p(count) = plot(t_span,CV2_eps_S_traj(:,count), 'Color', colors_rgb(count,:), 'LineWidth',2.5);hold on;

    axis([0 this_t_end_plt 0 1]);
    xlabel({'Time (days)'},'Interpreter','Latex');
    ylabel({'Coefficient of Variation (Squared)'},'Interpreter','Latex');
    title('Susceptibility');

    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';


    if count==2
        %this_p(4) = plot(t_span,CV2_eps_S_traj_var_susc,'--','Color',colors_rgb(5,:),'LineWidth',2.5);
        txt = {'B'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
        box on
    end


    %% CV Transmissibility
    subplot(1,3,3);


    this_q(count) = plot(t_span,CV2_delta_I_traj(:,count), 'Color', colors_rgb(count,:), 'LineWidth',2.5);hold on;

    axis([0 this_t_end_plt 0 1]);
    xlabel({'Time (days)'},'Interpreter','Latex');
    ylabel({'Coefficient of Variation (Squared)'},'Interpreter','Latex');
    title('Transmissibility');

    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';


    if count==3
        txt = {'C'};
        text(0.025,1.035,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
        box on
    end
    
end


%% save figure ?
if save_fig_ans==1

    figures_location = '../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));

end