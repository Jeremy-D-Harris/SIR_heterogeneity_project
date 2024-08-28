% % plot bivaraite distributions over time

%%

clear all; close all; clc;

colors_rgb = [133,192,249; 74,112,188.5; 15,32,128; 166 166 166]/255;

save_fig_ans = 1;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure_correlations_distributionsovertime_matchR0_NegBin_bivariate';


%% loading

% define ending time:
this_t_end_plt = 200;

% spacing/scaling subplots
frac_spacing = 0.74;
frac_scaling = 0.21;

% var_eps = fliplr([0.66 0.74 0.82 0.895]); %linspace(0.999,0.49,201);
epsilon_level = fliplr([0.51 0.65 0.80 0.999]);

txt_eps1 = ['$\bar{\varepsilon}(t_0) = ', num2str(epsilon_level(1),'%1.2f'), '$ '];
txt_eps2 = ['$\bar{\varepsilon}(t_1) = ', num2str(epsilon_level(2),'%1.2f'), '$ '];
txt_eps3 = ['$\bar{\varepsilon}(t_2) = ', num2str(epsilon_level(3),'%1.2f'), '$ '];
txt_eps4 = ['$\bar{\varepsilon}(t_3) = ', num2str(epsilon_level(4),'%1.2f'), '$ '];


%%
% need to load all three cases first:
file_location = '../data/';

for count = 1:2

    if count==1
        
        % (2) positive correlations
        infile_positivecorrelations = 'NegBin_PositiveCorrelation.mat';
        load(strcat(file_location,infile_positivecorrelations));


    elseif count==2
        
        
        % (1) no correlations 
        infile_independent = 'NegBin_NoCorrelation.mat';
        load(strcat(file_location,infile_independent));


    % else
    %     % (3) negative correlations
    %     infile_negativecorrelations = 'GaussianNegativeCorrelation_0pt6_matchR0.mat';
    %     load(strcat(file_location,infile_negativecorrelations));

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


    for iii = 1:length(epsilon_level)
        epstime = find(results.mean_eps_S_traj<epsilon_level(iii)); %indices with epsilon level less than value
        associatedtimings(count,iii) = params.t_span(epstime(1));
        arrayS(count,iii,:,:) = results.S_traj_eps_delta_array(epstime(1),:,:);
        arrayI(count,iii,:,:) = results.I_traj_eps_delta_array(epstime(1),:,:);
    end
        

end


%% plotting

f1 = figure(1);
%set(f1, 'Position', [100 100 1000 1000]);
set(f1, 'Position', [100 100 1500 1000]);

sp1=subplot(3,4,2);

%old_pos = get(sp1, 'Position');
%set(sp1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
%old_pos = get(sp1, 'Position');


%%
%f1=subplot(3,4,1);

this_p(1) =plot(t_span, total_incidence_classic,'Color',colors_rgb(4,:),'LineWidth',2.5); hold on;

box('off');

for count = 1:2


    this_p(count+1) =plot(t_span, total_incidence(:,3-count),'Color',colors_rgb(3-count,:),'LineWidth',2.5); hold on;
    this_p(count+1).Color(4)=0.85;


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

%set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
%old_pos = get(f1, 'Position');
%% plot mean susceptibility
f1=subplot(3,4,3);
for count=1:2
    this_q=plot(t_span, mean_eps_S_traj(:,3-count),'Color',colors_rgb(3-count,:),'LineWidth',2.5); hold on;
    this_q.Color(4)=0.85;
end

for kk=1:length(epsilon_level)
    plot(t_span,epsilon_level(kk)*ones(size(t_span)),'k','LineWidth',1); hold on;
    if kk==1
        text(0.1,0.94,txt_eps1,'Interpreter','Latex','Units','normalized','FontSize',11);
    elseif kk==2
        text(0.1,0.7,txt_eps2,'Interpreter','Latex','Units','normalized','FontSize',11);
    elseif kk==3
        text(0.1,0.53,txt_eps3,'Interpreter','Latex','Units','normalized','FontSize',11);
    else
        text(0.1,0.38,txt_eps4,'Interpreter','Latex','Units','normalized','FontSize',11);
    end
end


axis([0 this_t_end_plt 0.2 1.1]);
xlabel('Time (days)'); ylabel({'Mean' ; 'Susceptiblity, $\bar{\varepsilon}(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

txt = {'B'};
text(0.025,1.045,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
%set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
%old_pos = get(f1, 'Position');

yticks([0.2 0.4 0.6 0.8 1]);
set(f1,'xticklabel',{[]},'yticklabel',[{''},{'0.4'},{'0.6'},{'0.8'},{'1.0'}]);
box('off');

txt = {'0.2'};
text(-0.11,0.06,txt,'Units','normalized',...
    'FontSize',14,'FontWeight','normal','FontName', 'Times');


%% plot effective transmission rate
f1=subplot(3,4,4);
for count=1:2
    this_q=plot(t_span,effective_transmission_rate_traj(:,3-count),'Color',colors_rgb(3-count,:),'LineWidth',2.5); hold on;
    this_q.Color(4)=0.85;
end
axis([0 this_t_end_plt 0.15 0.25]);
xlabel('Time (days)'); ylabel({'Effective'; 'Transmission'; 'Rate, $\beta\,\bar{\delta}_I(t)$'},'interpreter','latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


txt = {'C'};
text(0.025,1.045,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
%set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
%old_pos = get(f1, 'Position');


yticks([0.15 0.20 0.25]);
set(f1,'yticklabel',[{'0.15'},{'0.20'},{''}]);
box('off');

txt = {'0.25'};
text(-0.135,0.935,txt,'Units','normalized',...
    'FontSize',14,'FontWeight','normal','FontName', 'Times');

%% plotting bivariate distributions!
hh = 1;
for count=1:2
    for count_eps = 1:length(epsilon_level) %find epsilon levels

        this_epsilon_level = epsilon_level(count_eps);
        subplot(3,4,4+hh)
       
        eps_plt = params.eps;
        del_plt = params.del;

        climsS = [0 0.5];
        climsI = [0 0.04];
        %clims = empty;
        imagesc(eps_plt,del_plt,squeeze(arrayS(count,count_eps,:,:)), 'AlphaData', .9);
        hold on
        xline(epsilon_level(count_eps),LineWidth=2,LineStyle='--')
        yline(effective_transmission_rate_traj(associatedtimings(count,count_eps)+1),LineWidth=2,LineStyle='--')
        title({strcat('$\bar{\varepsilon}(t)=$ ',num2str(this_epsilon_level))},'interpreter','latex');
        % pcolor(eps,del,joint_S);
        %axis xy;
        set(gca,'YDir','normal');
        % colorbar;
        xlim([0 3]); ylim([0 3]);
        xlabel('susceptibility $\varepsilon$','interpreter','latex');

        if hh==1
            ylabel({'$\rho>0$';'transmissibility $\delta$'},'interpreter','latex');
        elseif hh == 5
            ylabel({'$\rho=0$';'transmissibility $\delta$'},'interpreter','latex');
        else
            ylabel({'transmissibility $\delta$'},'interpreter','latex');
        end
 hh = hh+1;

    end
end

sp1 = subplot(3,4,1);
delete(sp1);
legend_char1 = 'SIR';
legend_char2 = 'No Correlation, $\rho = 0$';
legend_char3 = 'Positive Correlation, $\rho > 0$';
legend(this_p,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Position',[0.17 0.82 0.12 0.1],'FontSize',14);







%% save figure ?
if save_fig_ans==1

    figures_location = '../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');
    saveas(f1,strcat(figures_location,figure_name),'pdf');

    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));

end
