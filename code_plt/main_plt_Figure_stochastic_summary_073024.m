% % plot stochastic simulation summary statistics

%%

clear all; close all; clc;

save_fig_ans = 0;
% save figure:
% 0 = no, 1 = yes

figure_name = 'Figure_stochastic_summary_073124';
figure_name2 = 'Figure_stochastic_histograms_091624';

%create color gradiets
c2 = [133,192,249]/255; % light blue
c1 = [15,32,128]/255; % dark blue

depth = 9;
[grad1,im]=colorGradient(c2,c1,depth);
% [grad2,im]=colorGradient(c2,c3,depth);

colors_rgb = grad1;

%% load results from file
file_location = '../data/';

load(append(file_location,"stochastic_output_31Jul2024125047.mat"))

%% Plotting
f1 = figure(1); set(f1, 'Position', [100 500 900 350]);

Outbreak_thresh = 50; %50 cases denotes an outbreak

%% plot outbreak proportion, (more than threshold infections) 
subplot(1,2,1);

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
txt = {'A'};
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
ax = subplot(1,2,2);
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

txt = {'B'};
text(0.025,1.035,txt,'Units','normalized','FontSize',18,'FontWeight','bold');


%display statistics
disp("median outbreak proportions (SIR,\rho=-0.6,\rho=0,\rho=0.6)")
median(outbreaklistproportion,2,"omitmissing")
disp("quantile outbreak proportions (SIR,\rho=-0.6,\rho=0,\rho=0.6) at [.025 .25 .50 .75 .975]")
quantile(outbreaklistproportion',[.025 .25 .50 .75 .975]) %quantiles for box and whiskers
disp("iqr outbreak proportions (SIR,\rho=-0.6,\rho=0,\rho=0.6)")
iqr(outbreaklistproportion')  %interquartile range


%% histogram of outbreak final size and final time
f2 = figure(2); set(f2, 'Position', [100 500 600 900]);

sizeVec = [stochSIR_sir.Final_Size stochSIR_negcorr.Final_Size stochSIR_uncorr.Final_Size stochSIR_poscorr.Final_Size];
timeVec = [stochSIR_sir.Final_Time stochSIR_negcorr.Final_Time stochSIR_uncorr.Final_Time stochSIR_poscorr.Final_Time];
sizexrange = [min(sizeVec) max(sizeVec)*1.1];
timexrange = [min(timeVec) max(timeVec)*1.1];
nbins = 100;
sedges = 0:max(sizeVec)/nbins:max(sizeVec);
tedges =  0:max(timeVec)/nbins:max(timeVec);

label_xpos = 0.025;
label_ypos = 1.06;
label_fontSz = 15;

leftxlabel = "Final outbreak size";
rightxlabel = "Outbreak duration (days)";


t = tiledlayout(4,2);
nexttile
%subplot(4,2,1)
histogram(stochSIR_sir.Final_Size,sedges,'FaceColor',colorsbox(1,:),'FaceAlpha',1)
xlim(sizexrange)
xlabel(leftxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"SIR",'Location','northeast','Box','off','FontSize',14)

txt = {'A'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

nexttile
%subplot(4,2,2)
histogram(stochSIR_sir.Final_Time,tedges,'FaceColor',colorsbox(1,:),'FaceAlpha',1)
xlim(timexrange)
xlabel(rightxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"SIR",'Location','northeast','Box','off','FontSize',14)

txt = {'B'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

nexttile
%subplot(4,2,3)
histogram(stochSIR_negcorr.Final_Size,sedges,'FaceColor',colorsbox(2,:),'FaceAlpha',1)
xlim(sizexrange)
xlabel(leftxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"\rho=-0.6",'Location','northeast','Box','off','FontSize',14)

txt = {'C'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

nexttile
%subplot(4,2,4)
histogram(stochSIR_negcorr.Final_Time,tedges,'FaceColor',colorsbox(2,:),'FaceAlpha',1)
xlim(timexrange)
xlabel(rightxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"\rho=-0.6",'Location','northeast','Box','off','FontSize',14)

txt = {'D'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

nexttile
%subplot(4,2,5)
histogram(stochSIR_uncorr.Final_Size,sedges,'FaceColor',colorsbox(3,:),'FaceAlpha',1)
xlim(sizexrange)
xlabel(leftxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"\rho=0",'Location','northeast','Box','off','FontSize',14)

txt = {'E'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

nexttile
%subplot(4,2,6)
histogram(stochSIR_uncorr.Final_Time,tedges,'FaceColor',colorsbox(3,:),'FaceAlpha',1)
xlim(timexrange)
xlabel(rightxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"\rho=0",'Location','northeast','Box','off','FontSize',14)

txt = {'F'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

nexttile
%subplot(4,2,7)
histogram(stochSIR_poscorr.Final_Size,sedges,'FaceColor',colorsbox(4,:),'FaceAlpha',1)
xlim(sizexrange)
xlabel(leftxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"\rho=0.6",'Location','northeast','Box','off','FontSize',14)


txt = {'G'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

nexttile
%subplot(4,2,8)
histogram(stochSIR_poscorr.Final_Time,tedges,'FaceColor',colorsbox(4,:),'FaceAlpha',1)
xlim(timexrange)
xlabel(rightxlabel,'FontSize',14,'FontName','Times');
dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend(dummyh,"\rho=0.6",'Location','northeast','Box','off','FontSize',14)

txt = {'H'};
text(label_xpos,label_ypos,txt,'Units','normalized','FontSize',label_fontSz,'FontWeight','bold');

t.TileSpacing = 'compact';
t.Padding = 'compact';



if save_fig_ans==1

    figures_location = '../figures/';
    saveas(f1,strcat(figures_location,figure_name),'epsc');

    fprintf('Figure saved:\n'); 
    fprintf(strcat(figure_name,'\n\n'));

    fprintf('Location:\n'); 
    fprintf(strcat(figures_location,'\n\n'));

     figures_location = '../figures/';
    saveas(f2,strcat(figures_location,figure_name2),'epsc');

    fprintf('Figure saved:\n'); 
    fprintf(strcat(figure_name2,'\n\n'));

    fprintf('Location:\n'); 
    fprintf(strcat(figures_location,'\n\n'));

end
