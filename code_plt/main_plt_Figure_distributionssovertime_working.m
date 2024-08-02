% % plot marginal distributions over time

%%

clear all; close all; clc;

colors_rgb = [133,192,249; 74,112,188.5; 15,32,128; 166 166 166]/255;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure_correlations_distributionsovertime_matchR0';

% define ending time:
this_t_end_plt = 200;

% spacing/scaling subplots
frac_spacing = 0.74;
frac_scaling = 0.21;

% var_eps = fliplr([0.66 0.74 0.82 0.895]); %linspace(0.999,0.49,201);
epsilon_level = fliplr([0.66 0.80 0.90 0.999]);

txt_eps1 = ['$\bar{\varepsilon}(t_0) = ', num2str(epsilon_level(1),'%1.2f'), '$ '];
txt_eps2 = ['$\bar{\varepsilon}(t_1) = ', num2str(epsilon_level(2),'%1.2f'), '$ '];
txt_eps3 = ['$\bar{\varepsilon}(t_2) = ', num2str(epsilon_level(3),'%1.2f'), '$ '];
txt_eps4 = ['$\bar{\varepsilon}(t_3) = ', num2str(epsilon_level(4),'%1.2f'), '$ '];


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
    bet = params.bet;
    eps = params.eps;
    del = params.del;

    %% load trajectories
    this_I_traj = results.I_traj;
    I_traj(:,count) = this_I_traj;

    this_S_traj = results.S_traj;
    S_traj(:,count) = this_S_traj;

    % mean susceptibility over time
    this_mean_eps_S_traj = results.mean_eps_S_traj;
    mean_eps_S_traj(:,count) = this_mean_eps_S_traj;

    % find mean susceptibility level for each case
    ind = find(this_mean_eps_S_traj<epsilon_level);
    this_ind_t_int = ind(1);
    ind_t_int(count) = this_ind_t_int;

    this_t_int = t_span(this_ind_t_int);
    t_int(count) = this_t_int;

    this_mean_eps_int = this_mean_eps_S_traj(this_ind_t_int);
    mean_eps_int(count) = this_mean_eps_int;

    % delta trajectory and corresponding intercept
    this_mean_delta_I_traj = results.mean_delta_I_traj;
    mean_delta_I_traj(:,count) = this_mean_delta_I_traj;

    this_mean_delta_int = this_mean_delta_I_traj(this_ind_t_int);
    mean_delta_int(count) = this_mean_delta_int;

    this_effective_transmission_rate_traj = bet*this_mean_delta_I_traj;
    effective_transmission_rate_traj(:,count) = this_effective_transmission_rate_traj;

    this_total_incidence = bet*(this_mean_delta_I_traj.*this_I_traj.*this_mean_eps_S_traj.*this_S_traj);
    total_incidence(:,count) = this_total_incidence;

    if count == 2
        total_incidence_classic = results_classic.total_incidence;
    end

    this_marginal_eps_S_traj = results.marginal_eps_S_traj(this_ind_t_int,:);
    marginal_susceptibility(count,:) = this_marginal_eps_S_traj;

    this_marginal_delta_I_traj = results.marginal_delta_I_traj(this_ind_t_int,:);

    marginal_transmissibility(count,:) = this_marginal_delta_I_traj;



end



%% Plotting
f1 = figure(1);
set(f1, 'Position', [100 100 1000 1000]);
% annotation('rectangle',[0.37 0.11 0.55 0.21], 'Color','black','LineWidth',2);
% annotation('rectangle',[0.37 0.32 0.55 0.21], 'Color','black','LineWidth',2);
% annotation('rectangle',[0.37 0.53 0.55 0.21], 'Color','black','LineWidth',2);
% annotation('rectangle',[0.37 0.74 0.55 0.25], 'Color','black','LineWidth',2);

sp1=subplot(4,3,1);

old_pos = get(sp1, 'Position');
set(sp1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
old_pos = get(sp1, 'Position');


%%
f1=subplot(4,3,4);

this_p(1) =plot(t_span, total_incidence_classic,'Color',colors_rgb(4,:),'LineWidth',2.5); hold on;

box('off');

for count = 1:3


    this_p(5-count) =plot(t_span, total_incidence(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    this_p(5-count).Color(4)=0.85;


end
axis([0 this_t_end_plt 0 0.025]);
xlabel('Time (days)'); ylabel({'Incident';'Infections, $\eta(t)$'},'Interpreter','Latex');
% ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

xticks([0 50 100 150 200]);
yticks([0 0.005 0.01 0.015 0.02 0.025]);
set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'0.005'},{'0.01'},{'0.015'},{'0.02'},{'0.025'}]);
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
for count=1:3
    this_q=plot(t_span, mean_eps_S_traj(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    this_q.Color(4)=0.85;
end

for kk=1:length(epsilon_level)
    plot(t_span,epsilon_level(kk)*ones(size(t_span)),'k','LineWidth',1); hold on;
    if kk==1
        text(0.1,0.89,txt_eps1,'Interpreter','Latex','Units','normalized','FontSize',11);
    elseif kk==2
        text(0.1,0.71,txt_eps2,'Interpreter','Latex','Units','normalized','FontSize',11);
    elseif kk==3
        text(0.1,0.55,txt_eps3,'Interpreter','Latex','Units','normalized','FontSize',11);
    else
        text(0.1,0.31,txt_eps4,'Interpreter','Latex','Units','normalized','FontSize',11);
    end
end


axis([0 this_t_end_plt 0.5 1.1]);
xlabel('Time (days)'); ylabel({'Mean' ; 'Susceptiblity, $\bar{\varepsilon}(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');

yticks([0.5 0.6 0.8 1 1.1]);
set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'0.6'},{'0.8'},{'1.0'},{' '}]);
box('off');


%% plot effective transmission rate
f1=subplot(4,3,10);
for count=1:3
    this_q=plot(t_span,effective_transmission_rate_traj(:,4-count),'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
    this_q.Color(4)=0.85;
end
axis([0 this_t_end_plt 0.15 0.3]);
xlabel('Time (days)'); ylabel({'Effective'; 'Transmission'; 'Rate, $\beta\,\bar{\delta}_I(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
old_pos = get(f1, 'Position');


yticks([0.15 0.2 0.25 0.3]);
set(f1,'yticklabel',[{'0.15'},{'0.20'},{'0.25'},{''}]);
box('off');

txt = {'0.30'};
text(-0.135,0.935,txt,'Units','normalized',...
    'FontSize',14,'FontWeight','normal','FontName', 'Times');


%% handle marginals: plot at different epsilon levels
% find mean susceptibility level for each case
for count_eps = 1:length(epsilon_level)

    this_epsilon_level = epsilon_level(count_eps);

    for count=1:3


        ind = find(mean_eps_S_traj(:,count)<this_epsilon_level);
        this_ind_t_int = ind(1);

        this_t_int = t_span(this_ind_t_int);
        %         t_int(count_eps,count) = this_t_int;

        this_mean_eps_S_int = mean_eps_S_traj(this_ind_t_int,count);
        %         mu_epsilon(count) = this_mu_epsilon;

        this_mean_delta_I_int = mean_delta_I_traj(this_ind_t_int,count);
        %         mu_delta(count) = this_mu_delta;
        mean_delta_I_int(count_eps,count) = this_mean_delta_I_int;

        this_marginal_susceptibility(1,:) = marginal_susceptibility(count,:);
        this_marginal_transmissibility(1,:) = marginal_transmissibility(count,:);

        collect_marginal_susceptibility(count_eps,count,:) = this_marginal_susceptibility;
        collect_marginal_transmissibility(count_eps,count,:) = this_marginal_transmissibility;

        % for plotting
        ind = find(eps > this_epsilon_level);
        this_ind_eps_int = ind(1);
        ind_eps_int(count) =  this_ind_eps_int;

        ind = find(del > this_mean_delta_I_int);
        this_ind_delta_int = ind(1);
        ind_delta_int(count) = this_ind_delta_int;
        collect_ind_delta_int(count_eps,count) = this_ind_delta_int;




    end

end


%% now plot marginals for susceptibility
for count_eps = 1:length(epsilon_level)

    this_epsilon_level = epsilon_level(count_eps);

    for count=1:3

        plt_this_susceptibility_distribution(1,:) = collect_marginal_susceptibility(count_eps,4-count,:);

        subplot(4,3,(2+3*(count_eps-1)));
        this_q=plot(eps,plt_this_susceptibility_distribution,'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
        % this_q.Color(4)=0.85;

        if count==1
            plot(this_epsilon_level,0,'.','Color',colors_rgb(4-count,:),'LineWidth',2,'MarkerSize',20);
        elseif count==2
            plot(this_epsilon_level,0,'o','Color',colors_rgb(4-count,:),'LineWidth',2,'MarkerSize',8);
        else
            plot(this_epsilon_level,0,'o','Color',colors_rgb(4-count,:),'LineWidth',2,'MarkerSize',12);
        end

        axis([0 3 0 1]);
        yticks([]);
        ylabel({['$g_S(t_',num2str(count_eps-1),',\varepsilon)$']},'interpreter','latex');
        f1=gca;
        f1.LineWidth = 1;
        f1.FontSize = 16;
        f1.FontWeight = 'normal';
        f1.FontName = 'Times';




        if count==3

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
            set(f1,'xticklabel',[{'0'},{'1'},{'2'},{'3'}]);

            xlabel('Susceptibility, $\varepsilon$','interpreter','latex');
        end


    end

end



%% now plot marginals for transmissibility
for count_eps = 1:length(epsilon_level)

    this_epsilon_level = epsilon_level(count_eps);

    subplot(4,3,(3+3*(count_eps-1)));
    for count=1:3

        this_collect_mean_delta_int = mean_delta_I_int(count_eps,4-count);

        %         collect_marginal_transmissibility(count_eps,count,:)
        plt_this_transmissibility_distribution(1,:) = collect_marginal_transmissibility(count_eps,4-count,:);

        this_q=plot(eps,plt_this_transmissibility_distribution,'Color',colors_rgb(4-count,:),'LineWidth',2.5); hold on;
        this_q.Color(4)=0.8;

        plot(this_collect_mean_delta_int, 0,'.','Color',colors_rgb(4-count,:),'LineWidth',2,'MarkerSize',30);

        axis([0 3 0 1]);
        yticks([]);
        %         xlabel('Susceptibility, $\varepsilon$','interpreter','latex');
        ylabel({['$h_I(t_',num2str(count_eps-1),',\delta)$']},'interpreter','latex');
        f1=gca;
        f1.LineWidth = 1;
        f1.FontSize = 16;
        f1.FontWeight = 'normal';
        f1.FontName = 'Times';


        if count==3

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

        end



    end

    if count_eps == 1

        txt = {'C'};
        text(0.025,1.045,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        title('Transmissibility','FontSize',16,'FontWeight','normal');

        text(0.6,0.85,txt_eps1,'Interpreter','Latex','Units','normalized','FontSize',14);

    elseif count_eps == 2

        text(0.6,0.85,txt_eps2,'Interpreter','Latex','Units','normalized','FontSize',14);


    elseif count_eps == 3

        text(0.6,0.85,txt_eps3,'Interpreter','Latex','Units','normalized','FontSize',14);


    else

        text(0.6,0.85,txt_eps4,'Interpreter','Latex','Units','normalized','FontSize',14);

        xticks([0 1 2 3]);
        set(f1,'xticklabel',[{'0'},{'1'},{'2'},{'3'}]);

        xlabel('Transmissibility, $\delta$','interpreter','latex');


    end

end


%%
sp1 = subplot(4,3,1);
delete(sp1);
legend_char1 = 'Classic SIR';
legend_char2 = 'Positive Correlation, $\rho > 0$';
legend_char3 = 'No Correlation, $\rho = 0$';
legend_char4 = 'Negative Correlation, $\rho < 0$';


legend(this_p,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Position',[0.17 0.82 0.12 0.1],'FontSize',14);

%% save figure ?
if save_fig_ans==1

    figures_location = '../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));

end