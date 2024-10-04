function [] = plt_distributions(eps_plt,del_plt,joint_S,joint_I,marg_eps_S,marg_delta_S,marg_eps_I,marg_delta_I, my_rgb_colors)

f1=figure; set(f1, 'Position', [10   50   850   700]);

%% top panels
subplot(2,2,1);
imagesc(eps_plt,del_plt,joint_S);
% pcolor(eps,del,joint_S);
%axis xy;
set(gca,'YDir','normal');
colorbar;
xlim([0 3]); ylim([0 3]);

xlabel('susceptibility $\varepsilon$','interpreter','latex');
ylabel({'transmissibility $\delta$'},'interpreter','latex');
zlabel('Population Fraction')
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

title('Joint Distribution in S');

subplot(2,2,2);
p(1)=plot(eps_plt,marg_eps_S,'.-','Color',my_rgb_colors(3,:),'LineWidth',2.5,'MarkerSize',20); hold on;
p(2)=plot(del_plt, marg_delta_S,'.-','Color',my_rgb_colors(2,:),'LineWidth',2.5,'MarkerSize',20); hold on;

xlim([0 3]);
title('Marginal Distributions in S');

xlabel('susceptibility $\varepsilon$ or transmissibility $\delta$','interpreter','latex');
ylabel({'Population Fraction'});

f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

legend(p,{'Susceptibility','Transmissibility'},'Location','NorthEast');

%% bottom panels
subplot(2,2,3);
imagesc(eps_plt,del_plt,joint_I);
% pcolor(eps,del,joint_S);
%axis xy;
set(gca,'YDir','normal');
colorbar;
% xlim([0 4.5]); ylim([0 4.5]);
xlim([0 3]); ylim([0 3]);

xlabel('susceptibility $\varepsilon$','interpreter','latex');
ylabel({'transmissibility $\delta$'},'interpreter','latex');
f2=gca;
f2.LineWidth = 1;
f2.FontSize = 14;
f2.FontWeight = 'normal';
f2.FontName = 'Times';

title('Joint Distribution in I');

subplot(2,2,4);
q(1)=plot(eps_plt,marg_eps_I,'.-','Color',my_rgb_colors(3,:),'LineWidth',2.5,'MarkerSize',20); hold on;
q(2)=plot(del_plt, marg_delta_I,'.-','Color',my_rgb_colors(2,:),'LineWidth',2.5,'MarkerSize',20); hold on;

xlim([0 3]); 
title('Marginal Distributions in I');

xlabel('susceptibility $\varepsilon$ or transmissibility $\delta$','interpreter','latex');
ylabel({'Probability'});

f2=gca;
f2.LineWidth = 1;
f2.FontSize = 14;
f2.FontWeight = 'normal';
f2.FontName = 'Times';

legend(q,{'Susceptibility','Transmissibility'},'Location','NorthEast');