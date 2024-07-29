function [] = plotInitDists(marg_e_s, marg_d_s, marg_e_i, marg_d_i, e, d, joint_s, joint_i, E, D, default_colors)
    f1=figure(1); set(f1, 'Position', [10   50   400   620]);

    subplot(2,1,1)
    q(1)=plot(d,marg_d_s,'--','Color',default_colors(1,:),'LineWidth',2); hold on;
    q(2)=plot(d, marg_d_i,'--','Color',default_colors(2,:),'LineWidth',2); hold on;

    xlabel('transmissibility d'); ylabel({'Population'; 'Fraction'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';

    legend(q,{'Marginal of Pop_S','Marginal of Pop_I'},'Location','SouthEast');


    subplot(2,1,2)
    q(1)=plot(e, marg_e_s,'--','Color',default_colors(3,:),'LineWidth',2); hold on;
    q(2)=plot(e, marg_e_i,'--','Color',default_colors(2,:),'LineWidth',2); hold on;

    xlabel('susceptibility e'); ylabel({'Population'; 'Fraction'});
    f2=gca;
    f2.LineWidth = 1;
    f2.FontSize = 14;
    f2.FontWeight = 'normal';
    f2.FontName = 'Times';

    legend(q,{'Marginal of Pop_S','Marginal of Pop_I'},'Location','SouthEast');

    f3 = figure(3); set(f3, 'Position', [450   50   400   620]);
    subplot(2,1,1)
    imagesc(e,d,joint_s');
    set(gca,'YDir','normal');
    colorbar;

    xlabel('susceptibility e'); ylabel({'transmissibility d'}); zlabel('Population Fraction')
    f4=gca;
    f4.LineWidth = 1;
    f4.FontSize = 14;
    f4.FontWeight = 'normal';
    f4.FontName = 'Times';

    title('Joint Distribution of Pop_S');


    subplot(2,1,2)
    imagesc(e,d,joint_i');
    set(gca,'YDir','normal');
    colorbar;

    xlabel('susceptibility e'); ylabel({'transmissibility d'}); zlabel('Population Fraction')
    f5=gca;
    f5.LineWidth = 1;
    f5.FontSize = 14;
    f5.FontWeight = 'normal';
    f5.FontName = 'Times';

    title('Joint Distribution of Pop_I');
    
end
