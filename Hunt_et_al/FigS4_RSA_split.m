clear all; close all;
[bd] = get_basedir;

addpath(genpath(fullfile(bd)));

%% load in info about units
load('all_units_info.mat');
nUnits = length(all_units);

for i = 1:nUnits
    area(i) = all_units{i}.UnitInfo.area;
    average_firing_rate(i) = all_units{i}.UnitInfo.average_firing_rate;
    isfrank(i) = all_units{i}.UnitInfo.isfrank;
    ismiles(i) = all_units{i}.UnitInfo.ismiles;
end

%% define areas. indices 1 = ACC; 2 = LPFC 3 = OFC
area_index = [1 2 3 5];
color_index = {'r' 'g' 'b' 'c' 'k'};
area_name = {'ACC' 'DLPFC' 'OFC' '' 'Other'};

area([87 find(average_firing_rate<1)]) = nan; %bad cells

ACCcells = find(area==1);
LPFCcells = find(area==2);
OFCcells = find(area==3);

tin = 100;
tout = 500;

%compute unitmat, which has dimensions nConditions*nUnits, for each brain region:
[unitmat{1},unitmats{1}] = calculate_RSA_matrices_updown(ACCcells,tin,tout);
[unitmat{2},unitmats{2}] = calculate_RSA_matrices_updown(LPFCcells,tin,tout);
[unitmat{3},unitmats{3}] = calculate_RSA_matrices_updown(OFCcells,tin,tout);

%% draw figure S4
figure(1);clf;
set(gcf,'Position',[-1744         688        1695         559]);
colormap hot
for i = 1:3
    subplot(1,3,mod(i,3)+1);
    ok = ~all(unitmat{i}==0); %exclude any cells that don't spike
    RSAmat{i} = corrcoef(normalise(unitmat{i}(:,ok),1)');
    imagesc(RSAmat{i});caxis([-0.3 0.3]);
    title(area_name{i});
    axis square; tidyfig;
    set(gca,'YTickLabel',[1:5],'YTick',1:40,'XTickLabel',[1:5],'XTick',1:40,'FontSize',10);
    line([10.5 10.5],[-0.5 40.5],'LineWidth',3,'Color','k');
    line([5.5 5.5],[-0.5 40.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    line([15.5 15.5],[-0.5 40.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    
    line([-0.5 40.5],[10.5 10.5],'LineWidth',3,'Color','k')
    line([-0.5 40.5],[5.5 5.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    line([-0.5 40.5],[15.5 15.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    
    line([30.5 30.5],[-0.5 40.5],'LineWidth',3,'Color','k');
    line([25.5 25.5],[-0.5 40.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    line([35.5 35.5],[-0.5 40.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    
    line([-0.5 40.5],[30.5 30.5],'LineWidth',3,'Color','k')
    line([-0.5 40.5],[25.5 25.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    line([-0.5 40.5],[35.5 35.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    
    line([-0.5 40.5],[20.5 20.5],'LineWidth',4.5,'Color','k');
    line([20.5 20.5],[-0.5 40.5],'LineWidth',4.5,'Color','k');
    xlabel(sprintf(['Prob.   Mag.   Prob.   Mag.   Prob.   Mag.   Prob.   Mag.\n' ...
        'Left               Right               Left               Right\n' ... 
        'Top                                    Bottom']))
    set(get(gca,'XLabel'),'FontSize',14);
    ylabel(sprintf(['Bottom                                 Top\n'...
        'Right              Left                Right             Left\n' ... 
        'Mag.   Prob.   Mag.   Prob.   Mag.   Prob.   Mag.   Prob.']))
    %ylabel('Left                              Right')
    set(get(gca,'YLabel'),'FontSize',14);
end

%% to compute sliding RSA matrix:
%{
for i = 1:3
    for ttt = 1:800
        dat = squeeze(unitmats{i}(:,:,ttt)');
        ok = ~all(dat==0); %exclude any cells that don't spike
        RSAmatt{i}(:,:,ttt) = corrcoef(normalise(dat(:,ok),1)');
    end
end
%}