function [] = plt_initdistributions(eps_plt,del_plt,joint_S,joint_I, corr_coeff, count)


% f1=figure; set(f1, 'Position', [10   50   850   700]);

%% top panels
subplot(2,5,count);
imagesc(eps_plt,del_plt,joint_S);
% pcolor(eps,del,joint_S);
%axis xy;
set(gca,'YDir','normal');
% colorbar;
xlim([0 3]); ylim([0 3]);

xlabel('susceptibility $\varepsilon$','interpreter','latex');
ylabel({'transmissibility $\delta$'},'interpreter','latex');
% zlabel('Population Fraction')
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

% title('Joint Distribution in S');
if count == 3
    title('$\rho =\,0 $','Interpreter','latex');
else
    title(['$\rho =\, $', num2str(corr_coeff,'%0.2f')],'Interpreter','latex');
end

if count == 1
    txt = {'A'};
    text(0.025,1.075,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

elseif count == 2
    txt = {'B'};
    text(0.025,1.075,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

elseif count == 3
    txt = {'C'};
    text(0.025,1.075,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

elseif count == 4
    txt = {'D'};
    text(0.025,1.075,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

else
    txt = {'E'};
    text(0.025,1.075,txt,'Units','normalized','FontSize',16,'FontWeight','bold');

end



%% bottom panels
subplot(2,5,(count+5));
imagesc(eps_plt,del_plt,joint_I);
% pcolor(eps,del,joint_S);
%axis xy;
set(gca,'YDir','normal');
% colorbar;
% xlim([0 4.5]); ylim([0 4.5]);
xlim([0 3]); ylim([0 3]);

xlabel('susceptibility $\varepsilon$','interpreter','latex');
ylabel({'transmissibility $\delta$'},'interpreter','latex');
f2=gca;
f2.LineWidth = 1;
f2.FontSize = 14;
f2.FontWeight = 'normal';
f2.FontName = 'Times';

% title('Joint in I');

