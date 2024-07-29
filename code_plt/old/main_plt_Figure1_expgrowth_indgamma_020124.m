% simulate SIR model with transmissibility & susceptibility variation
% discretized version with correlation

%%

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'NewFigure1_jointS_indgamma_expgrowth_020124';

% medium blue, black, grey, violet, green
colors_rgb = [74,112,188.5; 0 0 0; 166 166 166; 169,90,161;0, 158, 115]/255;

% define ending time:
this_t_end_plt = 200;

epsilon_levels = [1.01, 0.9985];


%%
% need to load all three cases first:
file_location = '../data/GammaInd/VaryingShapeParameterOfPotentialTrans/';


% independent discrete simulations
infile_discrete = 'dataIVGammaSCms1ss3mi1si10ICms1ss3mi1si2.mat';
load(strcat(file_location,infile_discrete));

% sus variation only
infile_svonly = 'dataIVGammaSusVarOnlySCss3ICss3.mat';
variable = 'data';
data_susvaronly= load(strcat(file_location,infile_svonly),variable);
mu_epsilon_S_traj_susvaronly = data_susvaronly.data.mu_epsilon_S_traj;
marg_S_e_traj_susvaronly = data_susvaronly.data.marg_S_e_traj;

% independent analytical
infile_analytical = 'gammaIndAnalS3M1.mat';
variable = 'data';
data_analytical = load(strcat(file_location,infile_analytical),variable);
data_analytical = data_analytical.data;
variable = 'params';
params_analytical = load(strcat(file_location,infile_analytical),variable);
params_analytical = params_analytical.params;



this_I_traj = data.I_traj;
this_I_traj_homogeneous = data.I_traj_v;
this_I_traj_susvaronly = data.I_traj_sv;
this_I_traj_analytical = data_analytical.I_traj;

this_S_traj = data.S_traj;
this_S_traj_homogeneous = data.S_traj_v;
this_S_traj_susvaronly = data.S_traj_sv;
this_S_traj_analytical = data_analytical.S_traj;

% finer mesh for finding intercepts
n_time_pts = 401;
t_span_finer = linspace(params.t_span(1),params.t_span(end),n_time_pts);
this_I_traj_interp = interp1(params.t_span,this_I_traj,t_span_finer,'spline');

% shift trajectory: find initial time point
ind = find(this_I_traj_interp>10^-5);
% ind = find(this_I_traj_interp>0);
this_ind_shift_finer = ind(1);
this_t_shift_finer = t_span_finer(this_ind_shift_finer)+0;%38.25;
this_t_shift_plt = params.t_span-this_t_shift_finer*ones(size(params.t_span));
this_t_shift_plt_finer = t_span_finer-this_t_shift_finer*ones(size(t_span_finer));

% now interpolate trajectories
% this_S_traj_interp = interp1(this_t_shift_plt,this_S_traj,this_t_shift_plt_finer);

% mean susceptibility over time
this_mu_epsilon_S_traj = data.mu_epsilon_S_traj;
this_mu_epsilon_S_traj_interp = interp1(this_t_shift_plt,this_mu_epsilon_S_traj,this_t_shift_plt_finer,'spline');
this_mu_epsilon_S_traj_analytical = data_analytical.mu_epsilon_S_traj;
this_mu_epsilon_S_traj_interp_analytical = interp1(this_t_shift_plt,this_mu_epsilon_S_traj_analytical,this_t_shift_plt_finer,'spline');


% find mean susceptibility level for each case
this_ind_t_int= zeros(length(epsilon_levels),1);
this_t_int = zeros(length(epsilon_levels),1);
for ee = 1:length(epsilon_levels)
    ind = find(this_mu_epsilon_S_traj_interp<epsilon_levels(ee));
    this_ind_t_int(ee) = ind(1);
    this_t_int(ee) = this_t_shift_plt_finer(this_ind_t_int(ee));
end

% ind = find(this_mu_epsilon_S_traj_interp<epsilon_level);
% this_ind_t_int = ind(1);
% this_t_int = this_t_shift_plt_finer(this_ind_t_int);

% mean intrinsic transmissibility over time
this_mu_delta_S_traj = data.mu_delta_S_traj;
this_mu_delta_S_traj_interp = interp1(this_t_shift_plt,this_mu_delta_S_traj,this_t_shift_plt_finer,'spline');


% this_mu_epsilon = mu_epsilon_S_traj_interp(this_ind_t_int);

% delta trajectory and corresponding intercept
this_mu_delta_I_traj = data.mu_delta_I_traj;
this_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_mu_delta_I_traj,this_t_shift_plt_finer,'spline');
this_mu_delta = this_mu_delta_I_traj_interp(this_ind_t_int);
this_mu_delta_I_traj_analytical = data_analytical.mu_delta_I_traj;
this_mu_delta_I_traj_interp_analytical = interp1(this_t_shift_plt,this_mu_delta_I_traj_analytical,this_t_shift_plt_finer,'spline');

% effective transmission rate
this_beta_mu_delta_I_traj = params.beta*this_mu_delta_I_traj;
this_beta_mu_delta_I_traj_interp = interp1(this_t_shift_plt,this_beta_mu_delta_I_traj,this_t_shift_plt_finer,'spline');

this_total_incidence = transpose(params.beta*(this_mu_delta_I_traj.*this_I_traj.*this_mu_epsilon_S_traj.*this_S_traj));
this_total_incidence_interp = interp1(this_t_shift_plt,this_total_incidence,this_t_shift_plt_finer,'spline');

this_total_incidence_analytical = transpose(params.beta*(this_mu_delta_I_traj_analytical.*this_I_traj_analytical.*this_mu_epsilon_S_traj_analytical.*this_S_traj_analytical));
this_total_incidence_analytical_interp = interp1(this_t_shift_plt,this_total_incidence_analytical,this_t_shift_plt_finer,'spline');

this_total_incidence_susvaronly = transpose(params.beta*(mu_epsilon_S_traj_susvaronly.*this_I_traj_susvaronly.*this_S_traj_susvaronly));
this_total_incidence_susvaronly_interp = interp1(this_t_shift_plt,this_total_incidence_susvaronly,this_t_shift_plt_finer,'spline');

this_total_incidence_homogeneous = params_v.beta*(this_I_traj_homogeneous.*this_S_traj_homogeneous);
this_total_incidence_homogeneous_interp = interp1(this_t_shift_plt,this_total_incidence_homogeneous,this_t_shift_plt_finer,'spline');


%% susceptibility marginal
marg_S_e_traj = data.marg_S_e_traj;
for kk=1:length(params.t_span)
    this_S_marg(1,:) = sum(data.S_traj_discrete(kk,:,:),3);
    marg_S_e_traj(kk,:) = (1/data.S_traj(kk))*this_S_marg;
end
for counter=1:length(marg_S_e_traj)
    
    this_marg_S_e_traj_interp = interp1(this_t_shift_plt,marg_S_e_traj(:,counter),this_t_shift_plt_finer);
    marg_S_e_traj_interp(:,counter) = this_marg_S_e_traj_interp;
    
end

this_marginal_susceptibility = marg_S_e_traj_interp(this_ind_t_int,:)';
this_marginal_susceptibility_init = marg_S_e_traj_interp(90,:)';

this_mean_susceptibility_init = sum(params.E(:,1).*this_marginal_susceptibility_init);
this_mean_susceptibility = sum(params.E(:,1).*this_marginal_susceptibility);

%% susceptibility marginal for analytical

this_mean_susceptibility_analytical = this_mu_epsilon_S_traj_interp_analytical(this_ind_t_int);
this_mean_susceptibility_init_analytical = this_mu_epsilon_S_traj_interp_analytical(90);

this_marginal_susceptibility_analytical = ones(length(params.E(:,1)),length(this_ind_t_int));
for tt = 1:length(this_ind_t_int)
    dist = gampdf(params.E(:,1),params_analytical.shape, this_mean_susceptibility_analytical(tt)./params_analytical.shape);
    dist = dist./sum(dist);
    this_marginal_susceptibility_analytical(:,tt) = dist;   
end
dist = gampdf(params.E(:,1),params_analytical.shape,this_mean_susceptibility_init_analytical/params_analytical.shape);
dist = dist./sum(dist);
this_marginal_susceptibility_init_analytical = dist;

%% susceptibility marginal for susceptibility variation only

for kk=1:length(params.t_span)
    this_S_marg_susvaronly(1,:) = sum(data_susvaronly.data.S_traj_discrete(kk,:,:),3);
    marg_S_e_traj_susvaronly(kk,:) = (1/data_susvaronly.data.S_traj(kk))*this_S_marg_susvaronly;
end


%% intrinsic transmissibility marginal
marg_S_d_traj = data.marg_S_d_traj;
for kk=1:length(params.t_span)
    this_S_d_marg(1,:) = sum(data.S_traj_discrete(kk,:,:),2);
    marg_S_d_traj(kk,:) = (1/data.S_traj(kk))*this_S_d_marg;
end
for counter=1:length(marg_S_d_traj)
    
    this_marg_S_d_traj_interp = interp1(this_t_shift_plt,marg_S_d_traj(:,counter),this_t_shift_plt_finer);
    marg_S_d_traj_interp(:,counter) = this_marg_S_d_traj_interp;
    
end

% evaluated during exponential growth
this_intrinsic_marginal_transmissibility = marg_S_d_traj_interp(this_ind_t_int,:)';
d_delta = params.E(3,1)-params.E(2,1);

this_intrinsic_mean_transmissibility = sum(params.E(:,1).*this_intrinsic_marginal_transmissibility);


%% transmissibility marginal
marg_I_d_traj = data.marg_I_d_traj;
for kk=1:length(params.t_span)
    this_I_marg(1,:) = sum(data.I_traj_discrete(kk,:,:),2);
    marg_I_d_traj(kk,:) = (1/data.I_traj(kk))*this_I_marg;
end
for counter=1:length(marg_I_d_traj)
    
    this_marg_I_d_traj_interp = interp1(this_t_shift_plt,marg_I_d_traj(:,counter),this_t_shift_plt_finer);
    marg_I_d_traj_interp(:,counter) = this_marg_I_d_traj_interp;
    
    
end

this_marginal_transmissibility = marg_I_d_traj_interp(this_ind_t_int,:)';
this_mean_transmissibility = sum(params.D(1,:)'.*this_marginal_transmissibility);

%% transmissibility marginal for analytical

this_mean_transmissibility_analytical = this_mu_delta_I_traj_interp_analytical(this_ind_t_int);
this_mean_transmissibility_init_analytical = this_mu_delta_I_traj_interp_analytical(90);

this_marginal_transmissibility_analytical = ones(length(params.D(1,:)'),length(this_ind_t_int));
for tt = 1:length(this_ind_t_int)
    dist = gampdf(params.D(1,:)',params_analytical.shape_d, this_mean_transmissibility_analytical(tt)./params_analytical.shape_d);
    dist = dist./sum(dist);
    this_marginal_transmissibility_analytical(:,tt) = dist;   
end
dist = gampdf(params.D(1,:)',params_analytical.shape_d,this_mean_transmissibility_init_analytical/params_analytical.shape_d);
dist = dist./sum(dist);
this_marginal_transmissibility_init_analytical = dist;

%% for plotting
this_ind_eps_int = zeros(length(epsilon_levels),1);
this_ind_delta_int = zeros(length(this_mu_delta),1);
for ee = 1:length(epsilon_levels)
    ind = find(params.E(:,1) > epsilon_levels(ee));
    this_ind_eps_int(ee) = ind(1);
    
    ind = find(params.E(:,1) > this_mu_delta(ee));
    this_ind_delta_int(ee) = ind(1);
end 
% ind = find(params.E(:,1) > epsilon_level);
% this_ind_eps_int = ind(1);

% ind = find(params.E(:,1) > this_mu_delta);
% this_ind_delta_int = ind(1);

% interpolate trajectories too:
this_S_traj_interp = interp1(this_t_shift_plt,this_S_traj,this_t_shift_plt_finer);


% joint distribution at t0% ~t1
joint_S_eps_delta(:,:) = data.S_traj_discrete(1,:,:); %23

%% calculate variances
Svar = zeros(length(params.t_span),1);
for tt = 1:length(params.t_span)
    Svar(tt) = sum((params.E(:,1)'- this_mu_epsilon_S_traj(tt)*ones(size(params.E(:,1)'))).^2.*marg_S_e_traj(tt,:));
end

Svar_susvaronly = zeros(length(params.t_span),1);
for tt = 1:length(params.t_span)
    Svar_susvaronly(tt) = sum((params.E(:,1)'- mu_epsilon_S_traj_susvaronly(tt)*ones(size(params.E(:,1)'))).^2.*marg_S_e_traj_susvaronly(tt,:));
end

Svar_analytical = zeros(length(params.t_span),1);
for tt = 1:length(params.t_span)
    dist = gampdf(params.E(:,1),params_analytical.shape,this_mu_epsilon_S_traj_analytical(tt)/params_analytical.shape);
    dist = dist./sum(dist);

    Svar_analytical(tt) = sum((params.E(:,1)'- this_mu_epsilon_S_traj_analytical(tt)*ones(size(params.E(:,1)'))).^2.*transpose(dist));
end

Sdvar = zeros(length(params.t_span),1);
for tt = 1:length(params.t_span)
    Sdvar(tt) = sum((params.D(1,:)- this_mu_delta_S_traj(tt)*ones(size(params.D(1,:)))).^2.*marg_S_d_traj(tt,:));
end


Ivar = zeros(length(params.t_span),1);
for tt = 1:length(params.t_span)
    Ivar(tt) = sum((params.D(1,:)- this_mu_delta_I_traj(tt)*ones(size(params.D(1,:)))).^2.*marg_I_d_traj(tt,:));
end

Ivar_analytical = zeros(length(params.t_span),1);
for tt = 1:length(params.t_span)
    dist = gampdf(params.E(:,1),params_analytical.shape_d,this_mu_delta_I_traj_analytical(tt)/params_analytical.shape_d);
    dist = dist./sum(dist);

    Ivar_analytical(tt) = sum((params.E(:,1)'- this_mu_delta_I_traj_analytical(tt)*ones(size(params.E(:,1)'))).^2.*transpose(dist));
end

%% Plotting
% f1 = figure(1); set(f1, 'Position', [100 500 1000 700]);
f1 = figure(1); set(f1, 'Position', [100 500 1200 700]);


%% panel A: incident infections
subplot(2,3,1);
this_p(1) = plot(this_t_shift_plt_finer, this_total_incidence_homogeneous_interp,'Color',colors_rgb(3,:),'LineWidth',2.5); hold on;
this_p(3) = plot(this_t_shift_plt_finer, this_total_incidence_interp,'Color',colors_rgb(1,:),'LineWidth',2.5); hold on;
this_p(2) = plot(this_t_shift_plt_finer, this_total_incidence_susvaronly_interp,'--','Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
this_p(4) = plot(this_t_shift_plt_finer, this_total_incidence_analytical_interp,':','Color',colors_rgb(5,:),'MarkerSize',6,'LineWidth',2);
for ee = 1:length(epsilon_levels)
    if ee==1
    
    plot(0,this_total_incidence_interp(this_ind_t_int(ee)),'o','Color',colors_rgb(2,:),'MarkerSize',8, 'LineWidth',2); hold on;
    text(0.01,0.074,'$t_0$','Interpreter','Latex','Units','normalized','FontSize',16)
    else
        
    plot(this_t_shift_plt_finer(this_ind_t_int(ee)),this_total_incidence_interp(this_ind_t_int(ee)),'o','Color',colors_rgb(4,:),'MarkerSize',8, 'LineWidth',2); hold on;
    text(0.23,0.074,'$t_1$','Interpreter','Latex','Units','normalized','FontSize',16)
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

legend_char1 = 'SIR Model';
legend_char2 = 'Variation in Susceptibility (Discrete Model)';
legend_char3 = 'Uncorrelated Variation (Discrete Model)';
legend_char4 = 'Uncorrelated Variation (Reduced Model)';

legend(this_p,{legend_char1,legend_char2, legend_char3, legend_char4}, 'Position', [0.2 0.846 0.09 0.07],'FontSize',10);



%% panel B: CV Susceptibility
subplot(2,3,2);
this_q(2) =plot(this_t_shift_plt, Svar./this_mu_epsilon_S_traj.^2,'Color',colors_rgb(1,:),'LineWidth',2.5);hold on;
this_q(1) = plot(this_t_shift_plt, Svar_susvaronly./mu_epsilon_S_traj_susvaronly.^2,'--','Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
this_q(3) =plot(this_t_shift_plt, Svar_analytical./this_mu_epsilon_S_traj_analytical.^2,':','Color',colors_rgb(5,:),'LineWidth',2);hold on;

text(0.7,0.6,'$\frac{1}{k} = \frac{1}{3}$','Interpreter','Latex','Units','normalized','FontSize',16,'Color',colors_rgb(5,:))
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

% legend_char1 = 'SIR model';
legend_char1 = 'Variation in Susceptibility (Discrete Model)';
legend_char2 = 'Uncorrelated Variation (Discrete Model)';
legend_char3 = 'Uncorrelated Variation (Reduced Model)';

legend(this_q,{legend_char1,legend_char2, legend_char3},'Position',[0.484 0.855 0.085 0.06],'FontSize',10);

%% panel C: CV transmissibility
subplot(2,3,3);
this_r(2) = plot(this_t_shift_plt, Ivar./this_mu_delta_I_traj.^2,'Color',colors_rgb(1,:),'LineWidth',2.5);hold on;
this_r(1) = plot(this_t_shift_plt, Sdvar./this_mu_delta_S_traj.^2,'--','Color',colors_rgb(3,:),'LineWidth',2.5);hold on;
% this_r(1).Color(4) = 0.85;
this_r(3) = plot(this_t_shift_plt, Ivar_analytical./this_mu_delta_I_traj_analytical.^2,':','Color',colors_rgb(5,:),'LineWidth',2);hold on;

text(0.7,0.125,'$\frac{1}{k} = \frac{1}{10}$','Interpreter','Latex','Units','normalized','FontSize',16,'Color',colors_rgb(5,:))

axis([0 this_t_end_plt 0 0.5]);
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

 
legend_char1 = 'Potential Transmissibility (Discrete Model)'; 
legend_char2 = 'Effective Transmissibility (Discrete Model)'; 
legend_char3 = 'Effective Transmissibility (Reduced Model)'; 

legend(this_r,{legend_char1,legend_char2,legend_char3}, 'interpreter','latex','Position',[0.755 0.855 0.095 0.06],'FontSize',10);



%% panel D: Initial joint distribution
subplot(2,3,4);
imagesc(params.E(:,1),params.E(:,1),transpose(joint_S_eps_delta));
set(gca,'YDir','normal');
colormap(parula);
f1=gca;
xticks([0 1 2 3]); 
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



%% panel E: marginals for susceptibility
subplot(2,3,5);
for ee = 1:length(epsilon_levels)

    if ee==1
        this_r(1) = plot(params.E(:,1),this_marginal_susceptibility(:,ee),'Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
        plot(epsilon_levels(ee),0,'o','Color',colors_rgb(2,:),'LineWidth',2,'MarkerSize',14);
        
    else
        this_r(2) = plot(params.E(:,1),this_marginal_susceptibility(:,ee),'--','Color',colors_rgb(4,:),'LineWidth',2.5); hold on;
        plot(epsilon_levels(ee),0,'o','Color',colors_rgb(4,:),'LineWidth',2,'MarkerSize',9);
        
        this_r(3) = plot(params.E(:,1),this_marginal_susceptibility_analytical(:,ee),':','Color',colors_rgb(5,:),'LineWidth',2);hold on;
        plot(epsilon_levels(ee),0,'.','Color',colors_rgb(5,:),'MarkerSize',26);
        
    end
    
    text(0.33,0.07,'$\overline{\varepsilon}$','Interpreter','Latex','Units','normalized','FontSize',16)
    
end

axis([0 3 0 0.8e-1]);
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

legend_char1 = 'Marginal Susceptibility, $g_S(0,\varepsilon)$ (Discrete Model)';
legend_char2 = 'Marginal Susceptibility at $t_1$, $g_S(t_1,\varepsilon)$ (Discrete Model)';
legend_char3 = 'Marginal Susceptibility at $t_1$ (Reduced Model)';

legend(this_r,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Position',[0.485 0.38 0.095 0.06],'FontSize',10);

%% panel F: marginals for transmissibility
subplot(2,3,6);
this_s(1) = plot(params.E(:,1),marg_S_d_traj(1,:),'Color',colors_rgb(3,:),'LineWidth', 2.5); hold on;
% this_s(2) = plot(params.E(:,1),this_marginal_transmissibility(:,1),'Color',colors_rgb(2,:),'LineWidth',2.5); hold on;
this_s(2) = plot(params.E(:,1),this_marginal_transmissibility(:,2),'--','Color',colors_rgb(4,:),'LineWidth',2.5); hold on;
this_s(3) = plot(params.E(:,1),this_marginal_transmissibility_analytical(:,2),':','Color',colors_rgb(5,:),'LineWidth',2);hold on;

% plot(this_mu_delta(1),0,'o','Color',colors_rgb(1,:),'LineWidth',2,'MarkerSize',14);
plot(this_mu_delta(2),0,'o','Color',colors_rgb(3,:),'LineWidth',2,'MarkerSize',12);
plot(this_mu_delta(2),0,'o','Color',colors_rgb(4,:),'LineWidth',2,'MarkerSize',9);
plot(this_mu_delta(2),0,'.','Color',colors_rgb(5,:),'MarkerSize',26);

text(0.33,0.07,'$\overline{\delta}$','Interpreter','Latex','Units','normalized','FontSize',16)

axis([0 3 0 .8e-1]);
% yticks([0.5e-6 1.5e-6]);
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

legend_char1 = 'Potential Marginal Transmissibility, $h_S(0,\delta)$ (Discrete Model)';
% legend_char2 = 'Effective Marginal Transmissibility, $h_I(0,\delta)$ (Discrete Model)';
legend_char2 = 'Effective Marginal Transmissibility, $h_I(t_1,\delta)$ (Discrete Model)';
legend_char3 = 'Effective Marginal Transmissibility at $t_1$ (Reduced Model)';

legend(this_s,{legend_char1,legend_char2, legend_char3}, 'Interpreter','Latex','Position',[0.78 0.375 0.095 0.06],'FontSize',10);

%%
if save_fig_ans==1
    
    figures_location = './figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(figures_location,'\n\n'));
    
end