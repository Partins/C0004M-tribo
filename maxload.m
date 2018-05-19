%% Assessing the load carrying capacity (LCC) for a 1D tilted pad bearing
% This script defines the dimensionless LCC for a 1D titled pad bearing
% where the ratio between the leading edge (hl) and the trailing edge (ht)
% film thickness is hl/ht = 1+k
%
% Author: Andreas Almqvist 2015
clear all; % Matlabs first commandment
%% INPUT
% Specifying the range of slope parameters k
k = linspace(0.01,4,10001);
%% POSTPROCESSING
% Evaluating the dimensionless LCC for the chosen range of k
LCC = log(1+k)./k.^2 - 2./(2+k)./k;
%% SOLVER
% Finding the maximum dimensionless LCC and the index for which this
% maximum occurs
[maxN,indmax] = max(LCC);
%% PREPROCESSING AND VISUALIZATION
% Displaying the maximum dimensionless LCC and the corresponding k
disp([LCC(indmax),k(indmax)])
% Visualizing the LCC as a function of k
figure(1); % Initializing the figure;
clf; % Clearing the figure from previous plots
set(gcf,'color',[1,1,1]); % Specifying the background color of the fig
hold on; % Enabling plotting more than one line
% Variation in N with k
plot(k,LCC,'k-','linewidth',1);
% Illustrating the location of the max
plot(k(indmax),LCC(indmax),...
'ro','markersize',5,'markerfacecolor','k');
plot(k(k<k(indmax)),LCC(indmax).*ones(1,indmax-1),'r--');
plot([k(indmax),k(indmax)],[0,LCC(indmax)],'r--');
set(gca,'xtick',[0,1,2,3,4]); % Specifying xticks
set(gca,'ytick',[0,0.01,0.02,0.03]);
set(gca,'fontsize',12,'box','on'); % Setting the fontsize
% Writing labels for the maximum local on the axis
text(-0.05,LCC(indmax),['$',sprintf('%1.4f',LCC(indmax)),'$'],...
'horizontalalignment','right','interpreter','latex',...
'fontsize',16,'color','r');
text(k(indmax),-0.0012,['$',sprintf('%1.3f',k(indmax)),'$'],...
'horizontalalignment','center','verticalalignment','top',...
'interpreter','latex','fontsize',16,'color','r');
xlabel('$k$','interpreter','latex','fontsize',16); % Creating x-axis label
ylabel('$\overline{N}$','interpreter','latex','fontsize',16);