clear all; close all;
bd = get_basedir;
addpath(genpath(bd));

%% load in neuronal analyses
load(fullfile(bd,'neuronal_regression_results','info_gathering_cue123_alltrials_attentional_final.mat'));

% n.b. these results can also be obtained by running
% out = InfoGathering_Neurons_Cue1_regression(0);
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

thresh = 1.64; %norminv(0.05)

v2code = squeeze(tt2(:,6,:));
v1code = squeeze(tt2(:,5,:));
vDcode = squeeze(tt2(:,31,:));
sigsigA = abs(v2code)>thresh&abs(v1code)>thresh&abs(vDcode)>thresh&sign(v2code)~=sign(v1code);
for i = 1:3
    nSig(i,1) = sum(any(sigsigA(:,area==i)));
    meanSig(i,1) = mean(any(sigsigA(:,area==i)));
end


v3code = squeeze(tt3(:,3,:));
v1code = squeeze(tt3(:,1,:));
vDcode = squeeze(tt3(:,32,:));
sigsigB = abs(v3code)>thresh&abs(v1code)>thresh&abs(vDcode)>thresh&sign(v3code)~=sign(v1code);
for i = 1:3
    nSig(i,2) = sum(any(sigsigB(:,area==i)));
    meanSig(i,2) = mean(any(sigsigB(:,area==i)));
end


v3code = squeeze(tt3(:,3,:));
v2code = squeeze(tt3(:,2,:));
vDcode = squeeze(tt3(:,33,:));
sigsigC = abs(v3code)>thresh&abs(v2code)>thresh&abs(vDcode)>thresh&sign(v3code)~=sign(v2code);
for i = 1:3
    nSig(i,3) = sum(any(sigsigC(:,area==i)));
    meanSig(i,3) = mean(any(sigsigC(:,area==i)));
end

figure;

area_index = [3 2 1];
color_index = {[0.2081 0.1663 0.5292] [0.1986 0.7214 0.6310] [0.7763 0.7831 0.0538]};
area_name = {'OFC' 'DLPFC' 'ACC'}; %to correspond to area_index

for i = 1:4
    if i == 1;
        subplot(2,3,1:3);
        tmp = sigsigA|sigsigB|sigsigC;
    elseif i==2;
        subplot(3,3,7);
        tmp = sigsigA;
    elseif i==3
        subplot(3,3,8);
        tmp = sigsigB;
    elseif i==4
        subplot(3,3,9);
        tmp = sigsigC;
    end
    
    for r = 1:3
        p = plot(timebins,conv(mean(tmp(:,area==area_index(r))')*100',ones(1,5)./5,'same'),'Color',color_index{r},'LineWidth',2);
        hold on;
        if i ==1
            set(p,'LineWidth',3);
        end
    end
    xlim([-100 400]);
    
    if i == 1;
        title('Value attended minus unattended coding (all cues)');
        xlabel('Time post-stimulus onset');
        ylabel('% single neurons passing criteria');
        tidyfig(gcf,10);box off;
        l = legend(area_name,'Location','northwest');
        set(l,'box','off');
    elseif i==2;
        title(sprintf('Cue 2 minus Cue 1,\n attribute trials'));
        xlabel('Time post-cue 2 onset');
        ylabel('% single neurons passing criteria');
        tidyfig(gcf,12);box off
    elseif i==3
        title(sprintf('Cue 3 minus Cue 1,\n option trials'));
        xlabel('Time post-cue 3 onset');
        tidyfig(gcf,12);box off
    elseif i==4
        title(sprintf('Cue 3 minus Cue 2,\n option trials'));
        xlabel('Time post-cue 3 onset');
        tidyfig(gcf,12);box off
    end
end

set(gcf,'Position',[1000  586 597 752]);

%% test for belief confirmation at cue 2/3 onset

thresh = 1.96; %norminv(0.05/2)
bccode2o = squeeze(tt2(:,13,:)); %option trials, cue 2
bccode2a = squeeze(tt2(:,14,:)); %attribute trials, cue 2
bccode3o = squeeze(tt3(:,15,:)); %option trials, cue 3
bccode3a = squeeze(tt3(:,16,:)); %attribute trials, cue 3

figure;set(gcf,'Position',[440   147   738   651]);
for i = 1:4
    if i == 1;
        subplot(2,2,1);
        tmp = abs(bccode2o)>thresh;
    elseif i==2;
        subplot(2,2,2);
        tmp = abs(bccode2a)>thresh;
    elseif i==3
        subplot(2,2,3);
        tmp = abs(bccode3o)>thresh;
    elseif i==4
        subplot(2,2,4);
        tmp = abs(bccode3a)>thresh;
    end

    for r = 1:3
        p = plot(timebins,conv(mean(tmp(:,area==area_index(r))')*100',ones(1,5)./5,'same'),'Color',color_index{r},'LineWidth',2);
        hold on;
    end
    
    xlim([-100 400]);ylim([0 22.5]);
    if i == 1;
        box off;
        title(sprintf('Belief confirmation at cue 2\n Option trials'));
        xlabel('Time post-cue 2  onset');
        ylabel('% single neurons significant');
        tidyfig(gcf,14);
        l = legend(area_name,'Location','northwest');
        set(l,'box','off');
        line([-100 400],[5 5],'Color','k','LineStyle','--','LineWidth',0.5);
    elseif i==2;
        title(sprintf('Belief confirmation at cue 2\n Attribute trials'));
        xlabel('Time post-cue 2 onset');
        tidyfig(gcf,14);box off
        line([-100 400],[5 5],'Color','k','LineStyle','--','LineWidth',0.5);
    elseif i==3
        title(sprintf('Belief confirmation at cue 3\n Option trials'));
        ylabel('% single neurons significant');
        xlabel('Time post-cue 3 onset');
        tidyfig(gcf,14);box off
        line([-100 400],[5 5],'Color','k','LineStyle','--','LineWidth',0.5);
    elseif i==4
        title(sprintf('Belief confirmation at cue 3\n Attribute trials'));
        xlabel('Time post-cue 3 onset');
        tidyfig(gcf,14);box off
        line([-100 400],[5 5],'Color','k','LineStyle','--','LineWidth',0.5);
    end
end
