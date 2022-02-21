clear all; close all;
[bd] = get_basedir;

addpath(genpath(bd));
rmpath(genpath(fullfile(bd,'nish_scripts')));
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
    't_stats*','timebins');

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

%% timebins_Cue1Long of interest for plotting
tbois = timebins_Cue1Long>=0&timebins_Cue1Long<=2000; 
tboir = timebins_ResponseLong>=-2000&timebins_ResponseLong<=0;


%%

tmax = [1200 1400 1600 1800 50000];
for i = 1:5;
    ccLR_stim(:,timebins_Cue1Long>tmax(i),i,:) = nan;
    ccLR_resp(:,timebins_ResponseLong<-tmax(i),i,:) = nan;
end

%% set up movie

drawmov = 1;

if drawmov
    fc = 1;
    myMovie = VideoWriter('subspace_movie_final.avi');
    myMovie.FrameRate = 30;
    myMovie.Quality = 100;
    myMovie.open;
end

%% plot stimulus locked, shift to response locked

vecx = timebins_Cue1Long(tbois);
vecz = squeeze(mean(ccLR_resp(1,tboir,:,:)-repmat(mean(ccLR_stim(1,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,:)),[1 sum(tboir) 1 1]),4));

movediff = sum(isnan(vecz));
shiftmove = [0 1./(1+exp(-(-6:0.02:6))) 1];
shiftmove = ceil(shiftmove'*movediff);

figure;set(gcf,'Position',[440   252   723   546]);
title('Stimulus locked');

for j = 1:length(movediff)
    veczs(:,j) = circshift(vecz(:,j),-shiftmove(end,j),1);
end
pp = plot(vecx,veczs,'LineWidth',3);cmap = reshape([pp(:).Color],[3 5])';
xlabel('Time after stimulus onset(ms)');
ylabel('Belief confirmation subspace');box off;
set(gca,'XTick',-500:500:2500,'Xlim',[-500 2500],...
    'YTick',-0.15:0.05:0.1,'YLim',[-0.15 0.1]);
tidyfig(gcf,14);
l = legend({'Response time < 1.2s' 'Response time 1.2-1.4s' 'Response time 1.4-1.6s' 'Response time 1.6-1.8s' 'Response time >1.8s'},'Location','SouthWest');
set(l,'Box','off');
if drawmov;
    for i = 1:100;
        if mod(fc,3)==1;drawnow;ff = getframe(gcf);writeVideo(myMovie,ff);end;fc=fc+1;
    end
end

vecx = timebins_ResponseLong(tboir);
if drawmov;
    for i = length(shiftmove):-1:1;
        for j = 1:length(movediff)
            veczs(:,j) = circshift(vecz(:,j),-shiftmove(i,j),1);
        end
        
        pp = plot(vecx,veczs,'LineWidth',3);
        
        xlabel('Time before response(ms)');
        ylabel('Belief confirmation subspace');box off;
        set(gca,'XTick',-2500:500:500,'Xlim',[-2500 500],...
            'YTick',-0.15:0.05:0.1,'YLim',[-0.15 0.1]);
        set(get(gca,'XAxis'),'Visible','off');tidyfig(gcf,14);
        if mod(fc,3)==1;drawnow;ff = getframe(gcf);writeVideo(myMovie,ff);end;fc=fc+1;
    end
else
    i = 1;
            for j = 1:length(movediff)
            veczs(:,j) = circshift(vecz(:,j),-shiftmove(i,j),1);
        end
        
        pp = plot(vecx,veczs,'LineWidth',3);
        
        xlabel('Time before response(ms)');
        ylabel('Belief confirmation subspace');box off;
        set(gca,'XTick',-2500:500:500,'Xlim',[-2500 500],...
            'YTick',-0.15:0.05:0.1,'YLim',[-0.15 0.2]);
        set(get(gca,'XAxis'),'Visible','off');tidyfig(gcf,14);
end

set(get(gca,'XAxis'),'Visible','on');
title('Response locked');

%% fade in split by R/L 

colshift = [0 1./(1+exp(-(-6:0.02:6))) 1];
vecz = squeeze(ccLR_resp(1,tboir,:,1)-repmat(mean(ccLR_stim(1,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,1)),[1 sum(tboir) 1 1]));
hold on;
pp2 = plot(vecx,vecz,'LineWidth',2);
vecz = squeeze(ccLR_resp(1,tboir,:,2)-repmat(mean(ccLR_stim(1,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,2)),[1 sum(tboir) 1 1]));
pp3 = plot(vecx,vecz);


for i = 1:length(colshift)
    for j = 1:5
        pp(j).Color = [cmap(j,1:3) 1-colshift(i)];
        pp2(j).Color = [cmap(j,1:3) colshift(i)];
        pp3(j).Color = [cmap(j,1:3) colshift(i)];
    end
    if i>(length(colshift)/2)
        title('Split by chose left/right');
    end
    if drawmov;if mod(fc,3)==1;drawnow;ff = getframe(gcf);writeVideo(myMovie,ff);end;fc=fc+1;end
end

%% make 3d plot

clf;
spin = diff(([0 1./(1+exp(-(-6:0.02:6))) 1]));
spinlast = [spin(1:146) ...
    diff(sum(spin(1:146)):0.5/length(spin):(1-sum(spin(1:146)))) ...
    spin(457:end)];

vecx = timebins_ResponseLong(tboir);
vecz = squeeze(ccLR_resp(1,tboir,:,1)-repmat(mean(ccLR_stim(1,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,1)),[1 sum(tboir) 1 1]));
vecy = squeeze(ccLR_resp(2,tboir,:,1)-repmat(mean(ccLR_stim(2,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,1)),[1 sum(tboir) 1 1]));
pp = plot3(vecx,vecy,vecz,'LineWidth',2);
for j = 1:5
    pp(j).Color = [cmap(j,1:3)];
end

xlabel('Time before response(ms)');
zlabel('Belief confirmation subspace');
ylabel('Chosen action subspace');

hold on;

%axis off;
% vecx = timebins_Cue1Long(tbois);
% vecz = squeeze(ccLR_stim(1,tbois,:,2)-repmat(mean(ccLR_stim(1,timebins_Cue1Long>=-100&timebins_Cue1Long<=0,:,2)),[1 sum(tbois) 1 1]));
% vecy = squeeze(ccLR_stim(3,tbois,:,2)-repmat(mean(ccLR_stim(3,timebins_Cue1Long>=-100&timebins_Cue1Long<=0,:,2)),[1 sum(tbois) 1 1]));
% plot3(vecx,vecy,vecz,'LineWidth',2);tidyfig;

vecx = timebins_ResponseLong(tboir);
vecz = squeeze(ccLR_resp(1,tboir,:,2)-repmat(mean(ccLR_stim(1,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,2)),[1 sum(tboir) 1 1]));
vecy = squeeze(ccLR_resp(2,tboir,:,2)-repmat(mean(ccLR_stim(2,timebins_Cue1Long>=-50&timebins_Cue1Long<=50,:,2)),[1 sum(tboir) 1 1]));
pp = plot3(vecx,vecy,vecz);
for j = 1:5
    pp(j).Color = [cmap(j,1:3)];
end

%% rotations
vav = [0 0 0 0 45 45+180 ];
vev = [0 0 90 90 60 30];
axoff = [0 1 0 0 1];

tidyfig(gcf,14);

va = vav(1);ve = vev(1);
for t = 1:(length(vav)-2)
    if t == 1
        title('Split by chose left/right');
    elseif t==2
        title('Rotate into action subspace');
    else
        title('');
    end

    if axoff(t)
        if(t==2)
            set(get(gca,'YAxis'),'Visible','off');
            set(get(gca,'ZAxis'),'Visible','off');
        else
            axis off;
        end
    end
    sc = 0;
    for i = spin;
        view(va,ve);drawnow;if drawmov;if mod(fc,3)==1;ff = getframe(gcf);writeVideo(myMovie,ff);end;fc=fc+1;end
        va = va + (vav(t+1)-vav(t))*i;
        ve = ve + (vev(t+1)-vev(t))*i;
        sc = sc+1;
        if ~axoff(t)
        set(gca,'YAxisLocation','right',...
            'XTick',-2500:500:500,'Xlim',[-2500 500],...
            'YTick',-0.2:0.05:0.2,'YLim',[-0.18 0.18],...
            'ZTick',-0.15:0.05:0.1,'ZLim',[-0.15 0.1]);
        end
    end
    if axoff(t)
        if(t==2)
            set(get(gca,'YAxis'),'Visible','on');
            set(get(gca,'ZAxis'),'Visible','on');
        else
            axis on;
        end
    end
end

t = t+1;
for i = spinlast;
    if axoff(t)
        set(get(gca,'YAxis'),'Visible','off');
        set(get(gca,'ZAxis'),'Visible','off');
    end
    view(va,ve);drawnow;if drawmov;if mod(fc,3)==1;ff = getframe(gcf);writeVideo(myMovie,ff);end;fc=fc+1;end
    va = va + (vav(t+1)-vav(t))*i;
    ve = ve + (vev(t+1)-vev(t))*i;
    if axoff(t)
        set(get(gca,'YAxis'),'Visible','on');
        set(get(gca,'ZAxis'),'Visible','on');
    end
end
grid on;
if drawmov
drawnow;ff = getframe(gcf);writeVideo(myMovie,ff);
myMovie.close;
end