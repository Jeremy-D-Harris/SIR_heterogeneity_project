function [] = plotInitDists(eps_plt, del, init_joint_S, init_joint_I,init_marg_eps_S,init_marg_delta_S,init_marg_eps_I,init_marg_delta_I, my_rgb_colors)
   


f1=figure(1); set(f1, 'Position', [10   50   840   350]);
    subplot(1,2,1);
    
    imagesc(eps_plt,del,init_joint_S);
    % pcolor(eps,del,joint_S);
    %axis xy;
    set(gca,'YDir','normal');
    colorbar;
    xlim([0 4.5]); ylim([0 4.5]);
    
    xlabel('susceptibility $\varepsilon$','interpreter','latex');
    ylabel({'transmissibility $\delta$'},'interpreter','latex');
    zlabel('Population Fraction')
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Joint Distribution in S');
    
    % cbh = colorbar();
    % set color range
    % caxis([0,.2])
    % set ticks
    % set(cbh, 'YTick', [0.001, 0.01, 0.05, .1], ...
    %     'YTickLabel', {'p=0.001', 'p=0.01', 'p=0.05', 'p=.1'})
    % cbh.Label.String = 'Population Density';
    
    
    % f1=figure(1);
    subplot(1,2,2);
    q(1)=plot(eps_plt,init_marg_eps_S,'.-','Color',my_rgb_colors(3,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    q(2)=plot(del, init_marg_delta_S,'.-','Color',my_rgb_colors(2,:),'LineWidth',2.5,'MarkerSize',20); hold on;
    % cat_margs = [marg_eps_S;marg_delta_S];
    % q=bar(eps,cat_margs);
    axis([0 4.5 0 1]);
    
    axis([0 4 0 1]);
    title('Marginal Distributions');
    xlabel('transmissibility d');
    
    xlabel('susceptibility $\varepsilon$ or transmissibility $\delta$','interpreter','latex');
    ylabel({'Population Fraction'});
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    legend(q,{'Susceptibility','Potential Transmissibility'},'Location','NorthEast');