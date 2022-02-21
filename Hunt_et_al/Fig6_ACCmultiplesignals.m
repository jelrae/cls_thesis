clear all; close all;
[bd] = get_basedir;

addpath(genpath(bd));

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
area_name = {'ACC' 'DLPFC' 'OFC' '' 'Other' ''};

area([87 find(average_firing_rate<1)]) = nan; %bad cells

ACCcells = find(area==1);
LPFCcells = find(area==2);
OFCcells = find(area==3);


%% load in neuronal analyses

load(fullfile(bd,'neuronal_regression_results','info_gathering_cue123_alltrials_attentional_final.mat'),...
    't_stats*','timebins','subspace*proj');

% n.b. these results can also be obtained by running
% out = InfoGathering_Neurons_Cue123_alltrials_attentional_value_final(0);
% fields_to_workspace(out);

%% grab t-statistics from these analyses, timelocked to Cue 1, Cue2, Cue 3 and response

for i = 1:length(t_stats1); tt1(:,:,i) = t_stats1{i}; tt1_odd(:,:,i) = t_stats1_odd{i}; tt1_even(:,:,i) = t_stats1_even{i}; end
for i = 1:length(t_stats2); tt2(:,:,i) = t_stats2{i}; tt2_odd(:,:,i) = t_stats2_odd{i}; tt2_even(:,:,i) = t_stats2_even{i};   end
for i = 1:length(t_stats3); tt3(:,:,i) = t_stats3{i}; tt3_odd(:,:,i) = t_stats3_odd{i}; tt3_even(:,:,i) = t_stats3_even{i};   end
for i = 1:length(t_statsR); ttR(:,:,i) = t_statsR{i}; ttR_odd(:,:,i) = t_statsR_odd{i}; ttR_even(:,:,i) = t_statsR_even{i};   end

%% define which brain region to look at, and which timebin (post cue-onset) to use for defining subspaces

aoi = 1; %area of interest - 1 = ACC, 2 = LPFC, 3= OFC
thisbin = find(timebins==300); % find bin corresponding to t=300ms


%% is there a belief confirmation subspace? try correlating the relevant regressors in all the different subregions

figure; set(gcf,'Position',[25 1202 1135 410]);
for this_aoi = 1:3
    subplot(1,3,mod(this_aoi,3)+1);
    
    %different T-stats for the different 'belief confirmation regressors:
    tmp = [squeeze(tt1(thisbin,19,area==this_aoi)) squeeze(tt2(thisbin,13,area==this_aoi)) squeeze(tt2(thisbin,14,area==this_aoi)) ...
        squeeze(tt3(thisbin,15,area==this_aoi)) squeeze(tt3(thisbin,16,area==this_aoi))];

    [r,p] = corrcoef(tmp);
    p = (p<0.005) + (p<0.01) + (p<0.05);
    p = 3 - (tril(p)-eye(5)*3);
    p = p(2:5,1:4);
    r = r(2:5,1:4);
    imagesc(p);caxis([-3 3]);colormap bone; box off;axis square
    title(area_name(this_aoi));
    tidyfig;
    
    for i = 1:4
        for j = 1:4
            if i<=j
                text(i-0.25,j,sprintf('%0.2f',r(j,i)));
            end
        end
    end
    
    set(gca,'XTick',1:4,'XTickLabel',{'EV1/EV5, Cue1' 'EV13, Cue2' 'EV14, Cue2' 'EV15, Cue3'},'XTickLabelRotation',90,'FontSize',9);
    set(gca,'YTick',1:4,'YTickLabel',{'EV13, Cue2' 'EV14, Cue2' 'EV15, Cue3' 'EV16, Cue3'});
end

%% now do Fisher's r-to-Z transformation on this, to compare across regions.

for this_aoi = 1:3
    %different T-stats for the different 'belief confirmation regressors:
    tmp = [squeeze(tt1(thisbin,19,area==this_aoi)) squeeze(tt2(thisbin,13,area==this_aoi)) squeeze(tt2(thisbin,14,area==this_aoi)) ...
        squeeze(tt3(thisbin,15,area==this_aoi)) squeeze(tt3(thisbin,16,area==this_aoi))];

    rr(:,:,this_aoi) = corrcoef(tmp); %correlation matrix
    nC(this_aoi) = sum(area==this_aoi); %number of cells
end

zACC_vs_OFC = (atanh(rr(:,:,1))-atanh(rr(:,:,3)))./sqrt(1/(sum(area==1)-3) + 1/(sum(area==3)-3));
pACC_vs_OFC = (1-normcdf(abs(zACC_vs_OFC),0,1))*2;
zACC_vs_DLPFC = (atanh(rr(:,:,1))-atanh(rr(:,:,2)))./sqrt(1/(sum(area==1)-3) + 1/(sum(area==2)-3));
pACC_vs_DLPFC = (1-normcdf(abs(zACC_vs_DLPFC),0,1))*2;

figure;

subplot(2,2,1);
p = (pACC_vs_OFC<0.005) + (pACC_vs_OFC<0.01) + (pACC_vs_OFC<0.05);
p = 3 - (tril(p)-eye(5)*3);
p = p(2:5,1:4);
z = zACC_vs_OFC(2:5,1:4);
z=tril(z); z(z==0)=nan; %just looks at lower left corner
imagesc(p);caxis([-3 3]);colormap bone; box off;axis square
imagesc(p);caxis([-3 3]);colormap bone; box off;axis square
title('ACC vs. OFC (Fisher r-to-Z transformation)');
set(gca,'XTickLabel','','YTickLabel','');
tidyfig;

for i = 1:4
    for j = 1:4
        if i<=j
            text(i-0.25,j,sprintf('%0.2f',z(j,i)));
        end
    end
end

subplot(2,2,2);
hist(z(:));xlim([-6 6]);box off; l = line([0 0],[0 3],'Color','k','LineStyle','--');
xlabel('Z-score (ACC > OFC)');

subplot(2,2,3);
p = (pACC_vs_DLPFC<0.005) + (pACC_vs_DLPFC<0.01) + (pACC_vs_DLPFC<0.05);
p = 3 - (tril(p)-eye(5)*3);
p = p(2:5,1:4);
z = zACC_vs_DLPFC(2:5,1:4); 
z=tril(z); z(z==0)=nan; %just looks at lower left corner
imagesc(p);caxis([-3 3]);colormap bone; box off;axis square
imagesc(p);caxis([-3 3]);colormap bone; box off;axis square
title(area_name(this_aoi));
title('ACC vs. DLPFC (Fisher r-to-Z transformation)');
set(gca,'XTickLabel','','YTickLabel','');
tidyfig;

for i = 1:4
    for j = 1:4
        if i<=j
            text(i-0.25,j,sprintf('%0.2f',z(j,i)));
        end
    end
end


subplot(2,2,4);
hist(z(:));xlim([-6 6]);box off; l = line([0 0],[0 3],'Color','k','LineStyle','--');
xlabel('Z-score (ACC > DLPFC)');


%% now load in data that has been sorted by reaction time:

load(fullfile(bd,'neuronal_regression_results','info_gathering_ramp_to_choice.mat'),'out');

% n.b. these results can also be obtained by running:
% out = InfoGathering_Neurons_Ramp_to_choice(0);

timebins_Cue1Long = out.timebins_Cue1Long;
timebins_ResponseLong = out.timebins_ResponseLong;

%% define different subspaces

%...on all trials:
%average of the four 'belief confirmation' regressors - 13 and 14 at cue 2, 15 and 16 at cue 3
subspace_belief = squeeze((tt3(thisbin,15,area==aoi)+tt3(thisbin,16,area==aoi)+...
    tt2(thisbin,14,area==aoi)+tt2(thisbin,13,area==aoi))/4); 
% average of the two regressors for 'choose left' minus 'choose right'
subspace_RL = (squeeze(ttR(30,17,area==aoi))+squeeze(ttR(30,18,area==aoi)))/2;
% average response on option trials:
subspace_O = squeeze(ttR(30,11,area==aoi));
% average response on attribute trials:
subspace_A = squeeze(ttR(30,12,area==aoi));

%...as above, but for even trials only:
subspace_belief_even = squeeze((tt3_even(thisbin,15,area==aoi)+tt3_even(thisbin,16,area==aoi)+...
    tt2_even(thisbin,14,area==aoi)+tt2_even(thisbin,13,area==aoi))/4);
subspace_RL_even = (squeeze(ttR_even(30,17,area==aoi))+squeeze(ttR_even(30,18,area==aoi)))/2;
subspace_O_even = squeeze(ttR_even(30,11,area==aoi));
subspace_A_even = squeeze(ttR_even(30,12,area==aoi));

%...as above, but for odd trials only:
subspace_belief_odd = squeeze((tt3_odd(thisbin,15,area==aoi)+tt3_odd(thisbin,16,area==aoi)+...
    tt2_odd(thisbin,14,area==aoi)+tt2_odd(thisbin,13,area==aoi))/4);
subspace_RL_odd = (squeeze(ttR_odd(30,17,area==aoi))+squeeze(ttR_odd(30,18,area==aoi)))/2;
subspace_O_odd = squeeze(ttR_odd(30,11,area==aoi));
subspace_A_odd = squeeze(ttR_odd(30,12,area==aoi));

%% now project the data (with different average RT lengths) onto these subspaces using OLS

%use subspaces from even trials, data from odd trials
dm_even = [demean(subspace_belief_even) demean(subspace_RL_even) ...
    ones(size(subspace_belief_even)) (subspace_A_even+subspace_O_even)/2];
nSub = size(dm_even,2); %number of Subspaces
data = out.av_Cue1L_odd(area==aoi,:,:); %Cue1 locked data from odd trials

%now do projection:
datr = squeeze(data(:,:)); 
cc_even = ols(datr,dm_even,eye(nSub));
cc_even = reshape(cc_even,[nSub size(data,2) size(data,3)]);

%use subspaces from odd trials, data from even trials:
dm_odd = [demean(subspace_belief_odd) demean(subspace_RL_odd) ...
    ones(size(subspace_belief_odd)) (subspace_A_odd+subspace_O_odd)/2];
data = out.av_Cue1L_even(area==aoi,:,:); %Cue1 locked data from even trials

%now do projection:
datr = squeeze(data(:,:));
cc_odd = ols(datr,dm_odd,eye(nSub));
cc_odd = reshape(cc_odd,[nSub size(data,2) size(data,3)]); 

%average these two projections:
cc_stim = (cc_even+cc_odd)/2;

%% do the same thing for stimulus locked data split by left/right final response

%even trial data projected onto odd trial subspace
data = out.av_Cue1RT_LR_even(area==aoi,:,:,:);
datr = squeeze(data(:,:));
badcells = any(isnan(datr),2); %a few cells contain nans due to lack of appropriate trials - trim these
ccLR_even = ols(datr(~badcells,:),dm_odd(~badcells,:),eye(nSub));
ccLR_even = reshape(ccLR_even,[nSub size(data,2) size(data,3) size(data,4)]);

%odd trial data projected onto even trial subspace
data = out.av_Cue1RT_LR_odd(area==aoi,:,:,:);
datr = squeeze(data(:,:));
badcells = any(isnan(datr),2); %a few cells contain nans due to lack of appropriate trials - trim these
ccLR_odd = ols(datr(~badcells,:),dm_even(~badcells,:),eye(nSub));
ccLR_odd = reshape(ccLR_odd,[nSub size(data,2) size(data,3) size(data,4)]);

%average the two
ccLR_stim = (ccLR_even+ccLR_odd)/2;

%% do the same thing for response locked data split by left/right final response

%even trial data projected onto odd trial subspace
data = out.av_RespRT_LR_even(area==aoi,:,:,:);
datr = squeeze(data(:,:));
badcells = any(isnan(datr),2); %a few cells contain nans due to lack of appropriate trials - trim these
ccLR_resp_even = ols(datr(~badcells,:),dm_odd(~badcells,:),eye(nSub));
ccLR_resp_even = reshape(ccLR_resp_even,[nSub size(data,2) size(data,3) size(data,4)]);

%odd trial data projected onto even trial subspace
data = out.av_RespRT_LR_odd(area==aoi,:,:,:);
datr = squeeze(data(:,:));
badcells = any(isnan(datr),2); %a few cells contain nans due to lack of appropriate trials - trim these
ccLR_resp_odd = ols(datr(~badcells,:),dm_even(~badcells,:),eye(nSub));
ccLR_resp_odd = reshape(ccLR_resp_odd,[nSub size(data,2) size(data,3) size(data,4)]);

%average the two
ccLR_resp = (ccLR_resp_even+ccLR_resp_odd)/2;

%% timebins of interest for plotting
tboi = timebins_Cue1Long>=-100&timebins_Cue1Long<=2200; 
tboir = timebins_ResponseLong>=-500&timebins_ResponseLong<=100;

%% now plot activity in belief confirmation subspace, timelocked to Cue1 onset

figure;
set(gcf,'Position',[-1346 709 899 661]);

%fill cc with nans so that we only plot data *up to time of choice*
tmax = [1200 1400 1600 1800 50000];
for i = 1:5
    cc_stim(:,timebins_Cue1Long>tmax(i),i) = nan; 
end
    
subplot(6,4,sort([5:4:24 6:4:24 7:4:24]))
plot(timebins_Cue1Long(tboi),squeeze(cc_stim(1,tboi,:)-...
    repmat(mean(cc_stim(1,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:),2),[1 sum(tboi) 1])),... %baseline correction
    'LineWidth',4);
ylabel(sprintf('%s Belief confirmation subspace (a.u.)',area_name{aoi}));tidyfig;box off;
xlabel('Time after Cue 1 onset (ms)');

l = legend({'Response time < 1.2s' 'Response time 1.2-1.4s' 'Response time 1.4-1.6s' ...
    'Response time 1.6-1.8s' 'Response time >1.8s'},'Location','NorthWest');
set(l,'Box','off');
xlim([0 2100]);
ylim([-0.08 0.06]);
tidyfig(gcf,12);

%% now plot the boxplot of Cue2, Cue3, Cue4 onsets
subplot(6,4,1:3);
boxplot([out.Cue2_ON_list-out.Cue1_ON_list;out.Cue3_ON_list-out.Cue1_ON_list;...
    out.Cue4_ON_list-out.Cue1_ON_list]',...
    'symbol','','orientation','horizontal','plotstyle','compact')
axis off;
xlim([-100 2100]);

%% now plot activity in belief confirmation subspace, timelocked to response
tmax = [1200 1400 1600 1800 50000];
for i = 1:5;
    ccLR_stim(:,timebins_Cue1Long>tmax(i),i,:) = nan;
    ccLR_resp(:,timebins_ResponseLong<-tmax(i),i,:) = nan;
end

vecx = timebins_ResponseLong(tboir);
%averaged over left/right responses, and debaselined with respect to pre-stimulus period as before
vecz = squeeze(mean(ccLR_resp(1,tboir,:,:)-...
    repmat(mean(ccLR_stim(1,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,:)),[1 sum(tboir) 1 1]),4)); 


subplot(6,4,8:4:24);
plot(vecx,vecz,'LineWidth',4);
xlim([-500 100]);
ylim([-0.08 0.06]);
box off;
xlabel('Time until response (ms)');
tidyfig(gcf,12);

%% now plot in left/right subspace, timelocked to response

figure;hold on;
set(gcf,'Position',[-608   720   598   507]);
tboir = timebins_ResponseLong>=-2000&timebins_ResponseLong<=100;
vecx = timebins_ResponseLong(tboir);
vecy = squeeze(ccLR_resp(2,tboir,:,1)-repmat(mean(ccLR_stim(2,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,1)),[1 sum(tboir) 1 1]));
pp = plot(vecx,vecy,'LineWidth',4);
vecy = squeeze(ccLR_resp(2,tboir,:,2)-repmat(mean(ccLR_stim(2,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,2)),[1 sum(tboir) 1 1]));
p = plot(vecx,vecy,'LineWidth',4,'LineStyle','-.');
xlim([-1500 100]);
ylim([-0.12 0.12]);
for i = 1:5
   set(p(i),'Color',pp(i).Color);
end
xlabel('Time until response (ms)');
ylabel(sprintf('%s Chosen action subspace (a.u.)',area_name{aoi}))
tidyfig(gcf,12);

%% supplementary figure - individual correlations of the different regressors that make up the belief confirmation subspace

figure;

tmp = [squeeze(tt1(thisbin,19,area==aoi))  squeeze(tt2(thisbin,13,area==aoi)) squeeze(tt2(thisbin,14,area==aoi))...
    squeeze(tt3(thisbin,15,area==aoi)) squeeze(tt3(thisbin,16,area==aoi))];

labels = {sprintf('Value at Cue 1 \n(EV 1/EV 5 average)')  sprintf('Belief confirmation at \nCue 2, option (EV 13)') ...
    sprintf('Belief confirmation at \nCue 2, attribute (EV 14)')...
    sprintf('Belief confirmation at \nCue 3, option (EV 15)') ...
    sprintf('Belief confirmation at \nCue 3, attribute (EV 16)')};

set(gcf,'Position',[-1413 374  1218 1067]);
for i = 1:5
    for j = 1:5
        if i>j
            subplot(4,4,(i-2)*4+j);
            scatter_plus_fit(tmp(:,i),tmp(:,j),{'MarkerSize' 10 'Color' 'k'},{'LineWidth' 1 'Color' 'b'});

            if j==1
                ylabel(labels{i});
            else
                ylabel('');
            end
            if i==5
                xlabel(labels{j});
            else
                xlabel('');
            end
            tidyfig(gcf,10);box off;
        end
    end
end

%% 
