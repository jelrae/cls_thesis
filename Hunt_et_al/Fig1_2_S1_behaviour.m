clear all; close all;
[bd] = get_basedir;

%% load in behavioural data
addpath(genpath(fullfile(bd)));

load('Miles InfoGathering All Data');MCB=zInfoGatheringInfo;
load('Frank InfoGathering All Data');FCB=zInfoGatheringInfo;

clear zInfoGatheringInfo

%% process data
outM = InfoGathering_Behaviour(MCB);
outF = InfoGathering_Behaviour(FCB);

%% make figure 1

figure(1); clf
set(gcf,'Position',[10   979  976 726]);

%probability of choosing left/right on attribute trials
subplot(2,2,1);hold on;

diffmat = repmat(1:5,5,1)'-repmat(1:5,5,1);

plot(diffmat(:),outM.AvChoice12Att(:),'r.','MarkerSize',16);
dd = -4:0.01:4; b = outM.attchoice_b; l(1) = plot(dd,1./(1+exp(-(b(1)+dd.*b(2)))),'r','LineWidth',2);
plot(diffmat(:),outF.AvChoice12Att(:),'b.','MarkerSize',16);
dd = -4:0.01:4; b = outF.attchoice_b; l(2) = plot(dd,1./(1+exp(-(b(1)+dd.*b(2)))),'b','LineWidth',2);
xlabel(sprintf('Left - right picture\n rank difference'));
ylabel('p(Choice=Left)');
ll = legend(l,{'Subject M','Subject F'},'Location','SouthEast');
set(ll,'Box','off');
t = title('Attribute trials');
tidyfig(gcf,22);box off; axis square
set(t,'FontSize',16)

%probability of choosing 1/2 on option trials
subplot(2,2,2);hold on;

summat = repmat(-2:2,5,1)'+repmat(-2:2,5,1);

plot(summat(:),outM.AvChoice12Opt(:),'r.','MarkerSize',16);
dd = -4:0.01:4; b = outM.optchoice_b; l(1) = plot(dd,1./(1+exp(-(b(1)+dd.*b(2)))),'r','LineWidth',2);
plot(summat(:),outF.AvChoice12Opt(:),'b.','MarkerSize',16);
dd = -4:0.01:4; b = outF.optchoice_b; l(2) = plot(dd,1./(1+exp(-(b(1)+dd.*b(2)))),'b','LineWidth',2);
xlabel(sprintf('First + second picture\n rank sum (relative to average)'));
ylabel('p(Choice=Option 1)');
%ll = legend(l,{'Subject M','Subject F'},'Location','SouthEast');
%set(ll,'Box','off');
t = title('Option trials');
tidyfig(gcf,22);box off; axis square
set(t,'FontSize',16)

%%

figure;

subplot(2,2,1);
imagesc(1-outF.AvSaccDir12Att);
 axis square
xlabel('Cue 2 Rank');ylabel('Cue 1 Rank'); set(gca,'YAxisLocation','right')
colormap parula;caxis([0 0.85]);
tidyfig(gcf,22)
c = colorbar('WestOutside'); 
%axes(c);tidyfig(gcf,22);ylabel('p(Third saccade to option 1)','FontSize',14);

subplot(2,2,2);
imagesc(1-outM.AvSaccDir12Att);
 axis square
xlabel('Cue 2 Rank');ylabel('Cue 1 Rank'); set(gca,'YAxisLocation','right')
colormap parula;caxis([0 1]);
tidyfig(gcf,22)
c = colorbar('WestOutside'); 
%axes(c);tidyfig(gcf,22);ylabel('p(Third saccade to option 1)','FontSize',14);

subplot(2,2,3)
imagesc(outF.nPicView12Opt);
axis square
xlabel('Cue 2 Rank');ylabel('Cue 1 Rank'); set(gca,'YAxisLocation','right')
colormap parula;caxis([2 2.8]);
tidyfig(gcf,22)
c = colorbar('WestOutside');

subplot(2,2,4)
imagesc(outM.nPicView12Opt);
axis square
xlabel('Cue 2 Rank');ylabel('Cue 1 Rank'); set(gca,'YAxisLocation','right')
colormap parula;caxis([2.3 3.2]);
tidyfig(gcf,22)
c = colorbar('WestOutside');
%axes(c);tidyfig(gcf,22);ylabel(sprintf('Average number\n of cues viewed'),'FontSize',14);

%% Extended Data Figure 1a

figure;

RCmn = [outM.seq_choice_b outF.seq_choice_b];
RCse = [outM.seq_choice_b_stats.se outF.seq_choice_b_stats.se];

subplot(2,1,1);
h = barwitherr(RCse,RCmn);
box off; tidyfig(gcf,12);
l = legend({'Subject M' 'Subject F'});
set(l,'Box','Off');
set(gca,'XTickLabel',{'Left bias' 'Final cue' 'Cue n-1' 'Cue n-2' 'Cue n-3'},'XTickLabelRotation',30);
ylabel('Regression coefficient (+/- s.e.m.)');

RCmn = [[outM.choice_b(1:4); 0; outM.choice_b(5:8)] [outF.choice_b(1:4); 0; outF.choice_b(5:8)]];
RCse = [[outM.choice_b_stats.se(1:4); 0; outM.choice_b_stats.se(5:8)] [outF.choice_b_stats.se(1:4); 0; outF.choice_b_stats.se(5:8)]];

subplot(2,1,2);
h = barwitherr(RCse,RCmn);
box off; tidyfig(gcf,12);
set(gca,'XTickLabel',{'Left bias' 'L-R Probabilty' 'L-R Magnitude' 'Choose 1 bias' ''},'XTickLabelRotation',30);
xlabel('Option trials                                      Attribute trials     ');
ylabel('Regression coefficient (+/- s.e.m.)')
set(gcf,'Position',[320         526         705        1187]);

%% Extended Data Figure 1b/1c

figure;subplot(1,2,1);
NPVmn = [mean(outM.accuracyNPV*100); mean(outF.accuracyNPV*100)]';
NPVse = [std(outM.accuracyNPV*100)./sqrt(length(outM.accuracyNPV)); ...
    std(outF.accuracyNPV*100)./sqrt(length(outF.accuracyNPV))]';
barwitherr(NPVse,NPVmn);ylim([50 100]);
box off; tidyfig(gcf,12);
set(gca,'XTickLabel',[2 3 4]);
xlabel(sprintf('Number of pictures viewed\n before choice'));
ylabel('% optimal choices');

subplot(1,2,2);
DIFFmn = [mean(outM.mvd); mean(outF.mvd)]';
DIFFse = [std(outM.mvd)./sqrt(length(outM.mvd)); ...
    std(outF.mvd)./sqrt(length(outF.mvd))]';
barwitherr(DIFFse,DIFFmn);ylim([0.07 0.31]);
box off; tidyfig(gcf,12);
set(gca,'XTickLabel',[2 3 4]);
xlabel(sprintf('Number of pictures viewed\n before choice'));
ylabel('Inverse difficulty, |Left EV - Right EV| (a.u.)');
set(gcf,'Position',[1104  566  532 503]);

