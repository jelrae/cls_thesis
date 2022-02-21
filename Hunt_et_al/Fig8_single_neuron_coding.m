clear all;close all;
bd = get_basedir;
addpath(genpath(bd));

%% load in neuronal analyses
load(fullfile(bd,'neuronal_regression_results','info_gathering_cue1_both.mat'));

% n.b. these results can also be obtained by running
% out = InfoGathering_Neurons_Cue1_regression(0);
% fields_to_workspace(out);
    
%% define areas. indices 1 = ACC; 2 = LPFC 3 = OFC

load('all_units_info.mat','all_units');
for i = 1:length(all_units)
    area(i) = all_units{i}.UnitInfo.area;
    average_firing_rate(i) =  all_units{i}.UnitInfo.average_firing_rate;
end
area_index = [3 2 1];
color_index = {[0.2081 0.1663 0.5292] [0.1986 0.7214 0.6310] [0.7763 0.7831 0.0538]};
area_name = {'OFC' 'DLPFC' 'ACC'}; %to correspond to area_index

area([87 find(average_firing_rate<1)]) = nan; %bad cells

%% bar plots of number of significant neurons for different variables of interest

% relevant contrasts of parameter estimates from regression
% (see InfoGathering_Neurons_Cue1_regression for details):
% 26 - value
% 41 - Left MINUS Right value
% 51 - Mag MINUS Prob
% 42 - Top MINUS Bottom value

figure;subplot(2,1,1);
set(gcf,'Position',[292   815   700   966]);
for i = 1:length(area_index)
    toplot(1,i) = mean(significant_neuron(area==area_index(i),26))*100;
    toplot(2,i) = mean(significant_neuron(area==area_index(i),41))*100;
    toplot(3,i) = mean(significant_neuron(area==area_index(i),51))*100;
    toplot(4,i) = mean(significant_neuron(area==area_index(i),42))*100;
end

b = bar(toplot);

legend(area_name,'Box','off');
ylabel('Significant Neurons (%)');
set(gca,'XTickLabel',{'Cue value' 'Action value' 'Attribute value' 'Spatial value'});
line([0 5],[5 5],'LineStyle','--','Color','k');
box off;
clear l
tidyfig(gcf,10);

%% plots a couple of example single neurons - units 681 and 705
out = plot_neurons_average_peristim(681,[]);fields_to_workspace(out);
cc = [[0.1:0.2:0.9]' 0.5*ones(5,1) [0.9:-0.2:0.1]'];

ymin = min(min([Cue1Lmn Cue1Rmn]))/1.1;
ymax = max(max([Cue1Lmn Cue1Rmn]))*1.3;

subplot(4,2,5);hold on;
title('Cue 1 value, on left');tidyfig;xlim([-100 400]);
for i = 1:5;
    p(i) = plot(timeind,Cue1Lmn(i,:),'Color',cc(i,:),'LineWidth',2);
end
l = legend(p(5:-1:1),{'Best', ' ', ' ', ' ', 'Worst'},'Location','northwest');
set(l,'Box','off');
ylim([ymin ymax]);xlim([-100 400]);ylabel('Firing rate (Hz)');
tidyfig(gcf,10);

subplot(4,2,6);hold on;
title('Cue 1 value, on right');tidyfig;xlim([-100 400]);
for i = 1:5;
    plot(timeind,Cue1Rmn(i,:),'Color',cc(i,:),'LineWidth',2);
end
ylim([ymin ymax]);xlim([-100 400]);
tidyfig(gcf,10);

out = plot_neurons_average_peristim(705,[]);fields_to_workspace(out);

ymin = min(min([Cue1Pmn Cue1Mmn]))/1.1;
ymax = max(max([Cue1Pmn Cue1Mmn]))*1.3;

subplot(4,2,7);hold on;
title('Cue 1 value, probability');tidyfig;xlim([-100 400]);
for i = 1:5;
    p(i) = plot(timeind,Cue1Pmn(i,:),'Color',cc(i,:),'LineWidth',2);
end
l = legend(p(5:-1:1),{'Best', ' ', ' ', ' ', 'Worst'},'Location','northwest');
set(l,'Box','off');
ylim([ymin ymax]);xlim([-100 400]);ylabel('Firing rate (Hz)');
tidyfig(gcf,10);

subplot(4,2,8);hold on;
title('Cue 1 value, magnitude');tidyfig;xlim([-100 400]);
for i = 1:5;
    plot(timeind,Cue1Mmn(i,:),'Color',cc(i,:),'LineWidth',2);
end
ylim([ymin ymax]);xlim([-100 400]);
tidyfig(gcf,10);

clear l

%% plot coefficient of partial determination for cue value, left/right action value, 
%  attribute value, top/bottom spatial value

figure;
subplot(2,2,1);
for i = 1:length(area_index);
    l(i,:) = plotmse(squeeze(varCPD(area == area_index(i),1,:))*100,1,...
        timebins,'Color',color_index{i});
    xlim([-100 500]);
end
legend(l(:,1),area_name,'Location','NorthWest','Box','Off');
ylabel(sprintf('CPD, %s \n(mean +/- s.e.)','%'));
title('Cue value');
tidyfig(gcf,10);

subplot(2,2,2);
for i = 1:length(area_index);
    plotmse(squeeze(varCPD(area == area_index(i),2,:))*100,1,timebins,'Color',color_index{i});
    ylim([0 1.5]);
    xlim([-100 500]);
end
ylabel(sprintf('CPD, %s \n(mean +/- s.e.)','%'));
title('Action value');
tidyfig(gcf,10);

subplot(2,2,3);
for i = 1:length(area_index);
    plotmse(squeeze(varCPD(area == area_index(i),3,:))*100,1,timebins,'Color',color_index{i});
    ylim([0 1.5]);
    xlim([-100 500]);
end
ylabel(sprintf('CPD, %s \n(mean +/- s.e.)','%'));
xlabel('Time relative to Cue 1 onset (ms)');
title('Attribute value');
tidyfig(gcf,10);

subplot(2,2,4);
for i = 1:length(area_index);
    plotmse(squeeze(varCPD(area == area_index(i),4,:))*100,1,timebins,'Color',color_index{i});
    ylim([0 1.5]);
    xlim([-100 500]);
end
ylabel(sprintf('CPD, %s \n(mean +/- s.e.)','%'));
xlabel('Time relative to Cue 1 onset (ms)');
title('Spatial value');
tidyfig(gcf,10);

%% methods figure, showing correlation in design matrix

figure;
imagesc(mean(DM_crosscorr.^2,3));
xlabel('Explanatory variable #');
ylabel('Explanatory variable #');
caxis([0 0.3]);
set(gca,'YTick',1:13,'XTick',1:13)
c = colorbar;
ylabel(c,sprintf('Mean correlation between EVs, r^{%0.4g}',2))
tidyfig(gcf,11)

%% binomial test p-vals - tests whether region has more significant neurons than chance for each regressor
for i = 1:3
    bino_pvals(1,i) = (1-binocdf(sum(significant_neuron(area==area_index(i),26)),sum(area==area_index(i)),0.05))*2;
    bino_pvals(2,i) = (1-binocdf(sum(significant_neuron(area==area_index(i),41)),sum(area==area_index(i)),0.05))*2; %action value
    bino_pvals(3,i) = (1-binocdf(sum(significant_neuron(area==area_index(i),51)),sum(area==area_index(i)),0.05))*2; %attribute value
    bino_pvals(4,i) = (1-binocdf(sum(significant_neuron(area==area_index(i),42)),sum(area==area_index(i)),0.05))*2; %spatial value
end

% calculate fraction of significant neurons
for i = 1:3
    fraction_sig(1,i) = sum(significant_neuron(area==area_index(i),26))./sum(area==area_index(i));
    fraction_sig(2,i) = sum(significant_neuron(area==area_index(i),41))./sum(area==area_index(i)); %action value
    fraction_sig(3,i) = sum(significant_neuron(area==area_index(i),51))./sum(area==area_index(i)); %attribute value
    fraction_sig(4,i) = sum(significant_neuron(area==area_index(i),42))./sum(area==area_index(i)); %spatial value
end

%% pairwise chi2 tests of whether different regions have different fractions of significant neurons

[chistatACC_vs_OFC(1),pACC_vs_OFC(1)] = chi2([sum(significant_neuron(area==1,26)) sum(significant_neuron(area==3,26))],[sum(area==1) sum(area==3)]);
[chistatACC_vs_OFC(2),pACC_vs_OFC(2)] = chi2([sum(significant_neuron(area==1,41)) sum(significant_neuron(area==3,41))],[sum(area==1) sum(area==3)]);
[chistatACC_vs_OFC(3),pACC_vs_OFC(3)] = chi2([sum(significant_neuron(area==1,51)) sum(significant_neuron(area==3,51))],[sum(area==1) sum(area==3)]);
[chistatACC_vs_OFC(4),pACC_vs_OFC(4)] = chi2([sum(significant_neuron(area==1,42)) sum(significant_neuron(area==3,42))],[sum(area==1) sum(area==3)]);

[chistatACC_vs_LPFC(1),pACC_vs_LPFC(1)] = chi2([sum(significant_neuron(area==1,26)) sum(significant_neuron(area==2,26))],[sum(area==1) sum(area==2)]);
[chistatACC_vs_LPFC(2),pACC_vs_LPFC(2)] = chi2([sum(significant_neuron(area==1,41)) sum(significant_neuron(area==2,41))],[sum(area==1) sum(area==2)]);
[chistatACC_vs_LPFC(3),pACC_vs_LPFC(3)] = chi2([sum(significant_neuron(area==1,51)) sum(significant_neuron(area==2,51))],[sum(area==1) sum(area==2)]);
[chistatACC_vs_LPFC(4),pACC_vs_LPFC(4)] = chi2([sum(significant_neuron(area==1,42)) sum(significant_neuron(area==2,42))],[sum(area==1) sum(area==2)]);

[chistatOFC_vs_LPFC(1),pOFC_vs_LPFC(1)] = chi2([sum(significant_neuron(area==3,26)) sum(significant_neuron(area==2,26))],[sum(area==3) sum(area==2)]);
[chistatOFC_vs_LPFC(2),pOFC_vs_LPFC(2)] = chi2([sum(significant_neuron(area==3,41)) sum(significant_neuron(area==2,41))],[sum(area==3) sum(area==2)]);
[chistatOFC_vs_LPFC(3),pOFC_vs_LPFC(3)] = chi2([sum(significant_neuron(area==3,51)) sum(significant_neuron(area==2,51))],[sum(area==3) sum(area==2)]);
[chistatOFC_vs_LPFC(4),pOFC_vs_LPFC(4)] = chi2([sum(significant_neuron(area==3,42)) sum(significant_neuron(area==2,42))],[sum(area==3) sum(area==2)]);


%% 1-way ANOVA for whether neurons encode attended value at different latencies

tt = timebins(11:end);
[~,tmax] = max(squeeze(varCPD(:,1,11:end)),[],2);
tmax = tt(tmax);
tvec = [];
gvec = [];
for i = 1:3
    ind = (area == area_index(i)) & significant_neuron(:,26)';
    tvec = [tvec tmax(ind)];
    gvec = [gvec i*ones(1,sum(ind))];
end

anova1(tvec,gvec); %needs matlab stats toolbox
