clear all; close all;
[bd] = get_basedir;

addpath(genpath(bd));

%% load in neuronal analyses

load(fullfile(bd,'neuronal_regression_results','info_gathering_cue123_alltrials_attentional_final.mat'));

% n.b. these results can also be obtained by running
% out = InfoGathering_Neurons_Cue123_alltrials_attentional_value_final(0);
% fields_to_workspace(out);

 %% load in info about units

load('all_units_info.mat');
nUnits = length(all_units);

for i = 1:nUnits
    area(i) = all_units{i}.UnitInfo.area;
    average_firing_rate(i) = all_units{i}.UnitInfo.average_firing_rate;
    isfrank(i) = all_units{i}.UnitInfo.isfrank;
    ismiles(i) = all_units{i}.UnitInfo.ismiles;
end

area([87 653 find(average_firing_rate<1)]) = nan; %bad cells - 653 contains nans in T-stats

%% load in tstats 

for i = 1:length(t_stats1); tt1(:,:,i) = t_stats1{i}; end
for i = 1:length(t_stats2); tt2(:,:,i) = t_stats2{i}; end
for i = 1:length(t_stats3); tt3(:,:,i) = t_stats3{i}; end

%% define area/timebin of interest
aoi = 3; %area of interest - 1 = ACC, 2 = LPFC, 3= OFC
time_of_interest = 250; %time of interest post-cue, in ms
thisbin = find(timebins==time_of_interest); % find bin corresponding to t=time_of_interest

%% plot CPD across time

%Value of Cue 1, 2, and 3, timelocked to Cue 1
set(gcf,'Position',[109 214 1237 591]);
subplot(2,23,4:8);hold on;xlim([-100 400]);ylim([0 2.5]);
ll(1,:) =plotmse(squeeze(cpd1(area==aoi,1,:))*100,1,timebins,'Color','r');
ll(2,:) =plotmse(squeeze(cpd1(area==aoi,2,:))*100,1,timebins,'Color','g');
ll(3,:) =plotmse(squeeze(cpd1(area==aoi,3,:))*100,1,timebins,'Color','b');
xlabel(sprintf('Time after cue 1 (ms)'));
ylabel(sprintf('CPD (%%, mean +/- s.e.m.)'));tidyfig(gcf,14);
lll = legend(ll(:,1),{'Cue 1 Value','Cue 2 Value','Cue 3 Value'});
set(lll,'Box','off','Location','NorthWest');

%Value of Cue 1, 2, and 3, timelocked to Cue 2
subplot(2,23,10:14);hold on;xlim([-100 400]);ylim([0 2.5]);
plotmse(squeeze(cpd2(area==aoi,1,:))*100,1,timebins,'Color','r');
plotmse(squeeze(cpd2(area==aoi,2,:))*100,1,timebins,'Color','g');
plotmse(squeeze(cpd2(area==aoi,3,:))*100,1,timebins,'Color','b');
xlabel(sprintf('Time after cue 2 (ms)'));
set(gca,'YTickLabel',[]);
tidyfig(gcf,14);

%Value of Cue 1, 2, and 3, timelocked to Cue 3
subplot(2,23,16:20);hold on;xlim([-100 400]);ylim([0 2.5]);
plotmse(squeeze(cpd3(area==aoi,1,:))*100,1,timebins,'Color','r');
plotmse(squeeze(cpd3(area==aoi,2,:))*100,1,timebins,'Color','g');
plotmse(squeeze(cpd3(area==aoi,3,:))*100,1,timebins,'Color','b');
xlabel(sprintf('Time after cue 3 (ms)'));
set(gca,'YTickLabel',[]);
tidyfig(gcf,14);

%%%%%%%%%%%%%%% 
%parts b - f
%%%%%%%%%%%%%%%

axisoptions = {'Xlim',[-11 11],'Ylim',[-11 11]};

%figure;set(gcf,'Position',[36 477 1371 600]);
subplot(2,5,6);
scatter_plus_fit(squeeze(tt1(thisbin,19,area==aoi)),... %value of cue 1 at cue 1
                 squeeze(tt2(thisbin,20,area==aoi)),... %value of cue 2 at cue 2
                 [],[],axisoptions)
xlabel('V(Cue 1) at Cue 1');
ylabel('V(Cue 2) at Cue 2');
title('All trials');
axis square
tidyfig(gcf,14); box off;

subplot(2,5,7);
scatter_plus_fit(squeeze(tt2(thisbin,6,area==aoi)),... %value of cue 2 at cue 2, attribute trials
                 squeeze(tt2(thisbin,5,area==aoi)),... %value of cue 1 at cue 2, attribute trials
                 [],[],axisoptions)
xlabel('V(Cue 2) at Cue 2');
ylabel('V(Cue 1) at Cue 2');
title('Attribute trials');
axis square
tidyfig(gcf,14); box off;


subplot(2,5,8);
scatter_plus_fit(squeeze(tt3(thisbin,3,area==aoi)),... %value of cue 3 at cue 3, option trials
                 squeeze(tt3(thisbin,1,area==aoi)),... %value of cue 1 at cue 3, option trials
                 [],[],axisoptions);
xlabel('V(Cue 3) at Cue 3');
ylabel('V(Cue 1) at Cue 3');
title('Option trials');
axis square
tidyfig(gcf,14); box off

subplot(2,5,9);
scatter_plus_fit(squeeze(tt3(thisbin,3,area==aoi)),... %value of cue 3 at cue 3, option trials
                 squeeze(tt3(thisbin,2,area==aoi)),... %value of cue 2 at cue 3, option trials
                 [],[],axisoptions);
xlabel('V(Cue 3) at Cue 3');
ylabel('V(Cue 2) at Cue 3');
title('Option trials');
axis square
tidyfig(gcf,14); box off


subplot(2,5,10);
scatter_plus_fit(squeeze(tt3(thisbin,1,area==aoi)),... %value of cue 1 at cue 3, option trials
                 squeeze(tt3(thisbin,2,area==aoi)),... %value of cue 2 at cue 3, option trials
                 [],[],axisoptions);
xlabel('V(Cue 1) at Cue 3');
ylabel('V(Cue 2) at Cue 3');
title('Option trials');
axis square
tidyfig(gcf,14); box off

%% cross-correlation plots

figure;
set(gcf,'Position',[441 48 852   757]);
cid2 = [5 1 1 2];
cid1 = [6 2 3 3];
baseplot = 6;
jadd = [0 3 0 3];
iadd = [-1 -1 2 2];
cform_thresh = 0.2; %cluster-forming threshold for permutation test
nPerm = 1000; %number of permutations to perform
for count = 1:4
    CUEID1 = cid1(count); %rows
    CUEID2 = cid2(count); %columns
    clim = 0.3;
    
    [cmats{3,1},~,nullclust{3,1}] = cross_corr_tstats(t_stats1,t_stats1,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{3,2},~,nullclust{3,2}] = cross_corr_tstats(t_stats1,t_stats2,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{3,3},~,nullclust{3,3}] = cross_corr_tstats(t_stats1,t_stats3,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{2,1},~,nullclust{2,1}] = cross_corr_tstats(t_stats2,t_stats1,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{2,2},~,nullclust{2,2}] = cross_corr_tstats(t_stats2,t_stats2,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{2,3},~,nullclust{2,3}] = cross_corr_tstats(t_stats2,t_stats3,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{1,1},~,nullclust{1,1}] = cross_corr_tstats(t_stats3,t_stats1,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{1,2},~,nullclust{1,2}] = cross_corr_tstats(t_stats3,t_stats2,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    [cmats{1,3},~,nullclust{1,3}] = cross_corr_tstats(t_stats3,t_stats3,CUEID1,CUEID2,area==aoi,1,cform_thresh,nPerm);
    
    for i = 1:3
        for j = 1:3
            h= subplot(baseplot,baseplot,(i+iadd(count))*baseplot+jadd(count)+j);
            p = get(h,'Position');
            p(1) = p(1) + (2-j)*0.031;
            p(2) = p(2) + (i-2)*0.03;
            set(h,'Position',p);
            imagesc(timebins,timebins,cmats{i,j});
            caxis([-clim clim]);
            
            [b,l] = bwboundaries(abs(cmats{i,j})>cform_thresh);
            for c = 1:length(b)
                if sum(l(:)==c)>prctile(nullclust{i,j},99.9) %check if this cluster if larger than 99th prctile of null distribution
                    hold on;plot(timebins(b{c}(:,2)),timebins(b{c}(:,1)),'k','LineWidth',2)
                end
            end
            
            set(gca,'YDir','normal');
            if j==1
                ylabel(sprintf('Time after \nCue %0.0f (ms)',4-i));
            else
                set(gca,'YTickLabel',[]);
            end
            if i==3
                xlabel(sprintf('Time after \nCue %0.0f (ms)',j));
            else
                set(gca,'XTickLabel',[]);
            end
            tidyfig(gcf,10);
            
            
            axis square;
            if count == 1 && i == 2 && j == 1
                text(time_of_interest,time_of_interest,'b','FontSize',14,'FontName','Arial','FontWeight','normal','Color',[0 0 0]);
            elseif count == 2 && i == 2 && j == 1
                text(time_of_interest,time_of_interest,'b','FontSize',14,'FontName','Arial','FontWeight','normal','Color',[0 0 0]);
            elseif count == 1 && i == 2 && j == 2
                text(time_of_interest,time_of_interest,'c','FontSize',14,'FontName','Arial','FontWeight','normal','Color',[1 1 1]);
            elseif count == 3 && i == 1 && j == 3
                text(time_of_interest,time_of_interest,'d','FontSize',14,'FontName','Arial','FontWeight','normal','Color',[1 1 1]);
            elseif count == 4 && i == 1 && j == 3
                text(time_of_interest,time_of_interest,'e','FontSize',14,'FontName','Arial','FontWeight','normal','Color',[1 1 1]);
            elseif count == 2 && i == 1 && j == 3
                text(time_of_interest,time_of_interest,'f','FontSize',14,'FontName','Arial','FontWeight','normal','Color',[0 0 0]);
            end
        end
    end
    drawnow;
end

%% compare correlations between regions using Fisher's r-to-Z transform

thisbin = find(timebins==250); 
area_label = {'ACC' 'DLPFC' 'OFC'};
for aoi = 1:3
    [r,~,rlo,rup] = corrcoef(squeeze(tt2(thisbin,6,area==aoi)),squeeze(tt2(thisbin,5,area==aoi)));
    ccA(aoi) = r(1,2); ccAsem(aoi) = sqrt((1-r(1,2).^2)./(length(squeeze(tt2(thisbin,6,area==aoi)))-2));
    [r,~,rlo,rup] = corrcoef(squeeze(tt3(thisbin,3,area==aoi)),squeeze(tt3(thisbin,1,area==aoi)));
    ccB(aoi) = r(1,2); ccBsem(aoi) = sqrt((1-r(1,2).^2)./(length(squeeze(tt3(thisbin,3,area==aoi)))-2));
    [r,~,rlo,rup] = corrcoef(squeeze(tt3(thisbin,3,area==aoi)),squeeze(tt3(thisbin,2,area==aoi)));
    ccC(aoi) = r(1,2); ccCsem(aoi) = sqrt((1-r(1,2).^2)./(length(squeeze(tt3(thisbin,3,area==aoi)))-2));
end

for aoiA = 1:3
    for aoiB = 1:3
        %VCue2 minus VCue1 (@Cue 2), attribute trials
        zA(aoiA,aoiB) = (atanh(ccA(aoiA))-atanh(ccA(aoiB)))...
            ./sqrt(1/(sum(area==aoiA)-3) + 1/(sum(area==aoiB)-3));
        pA(aoiA,aoiB) = (1-normcdf(abs(zA(aoiA,aoiB)),0,1))*2;
        
        %VCue3 minus VCue1 (@Cue 3), option trials
        zB(aoiA,aoiB) = (atanh(ccB(aoiA))-atanh(ccB(aoiB)))...
            ./sqrt(1/(sum(area==aoiA)-3) + 1/(sum(area==aoiB)-3));
        pB(aoiA,aoiB) = (1-normcdf(abs(zB(aoiA,aoiB)),0,1))*2;
        
        %VCue3 minus VCue2 (@Cue 3), option trials
        zC(aoiA,aoiB) = (atanh(ccC(aoiA))-atanh(ccC(aoiB)))...
            ./sqrt(1/(sum(area==aoiA)-3) + 1/(sum(area==aoiB)-3));
        pC(aoiA,aoiB) = (1-normcdf(abs(zC(aoiA,aoiB)),0,1))*2;
    end
end

cc_stack = [ccA; ccB; ccC];
cc_err = [ccAsem; ccBsem; ccCsem];

%% plot differences between regions for value comparison figure

figure;
set(gcf,'Position',[1506  636 329  583]);
this_order = [3 2 1]; 
l = barwitherr(cc_err(:,this_order,:),cc_stack(:,this_order));
%ylim([-0.4 0.2]);
tidyfig(gcf,12); box off
set(gca,'XTickLabel',{'Attribute trials, Vcue2 vs. Vcue1' ...
                      'Option trials, Vcue3 vs. Vcue1' ...
                      'Option trials, Vcue3 vs. Vcue2'},'XTickLabelRotation',90);
ylabel('Correlation of regression coefficients across neurons (r)');
l = legend(area_label(this_order));
set(l,'box','off');
