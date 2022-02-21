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
area_name = {'ACC' 'DLPFC' 'OFC' '' 'Other' ''};

area([87 find(average_firing_rate<1)]) = nan; %bad cells
%area(isfrank) = nan; %to look only at Miles's RSA matrix, see Fig S2
%area(ismiles) = nan; %to look only at Frank's RSA matrix, see Fig S2

ACCcells = find(area==1);
LPFCcells = find(area==2);
OFCcells = find(area==3);

%% return the data for calculating the RSA matrices for each area
tin = 100;
tout = 500;

%compute unitmat, which has dimensions nConditions*nUnits, for each brain region:
[unitmat{1},unitmats{1}] = calculate_RSA_matrices(ACCcells,tin,tout);
[unitmat{2},unitmats{2}] = calculate_RSA_matrices(LPFCcells,tin,tout);
[unitmat{3},unitmats{3}] = calculate_RSA_matrices(OFCcells,tin,tout);

%% plot the RSA matrices - figure 4
figure(1);clf;
set(gcf,'Position',[100 100 1695 559]);
colormap hot
axorder = [3 2 1];
for i = 1:3
    subplot(1,3,axorder(i));
    ok = ~all(unitmat{i}==0); %exclude any cells that don't spike
    RSAmat{i} = corrcoef(normalise(unitmat{i}(:,ok),1)'); %this line of code calculates the RSA matrix
    RSAmat_pdist{i} = pdist(normalise(unitmat{i}(:,ok),1),'Correlation');
    imagesc(RSAmat{i});caxis([-0.3 0.3]);
    title(area_name{i});
    axis square; tidyfig;
    set(gca,'YTickLabel',[1:5],'YTick',1:20,'XTickLabel',[1:5],'XTick',1:20,'FontSize',10);
    line([10.5 10.5],[-0.5 20.5],'LineWidth',3,'Color','k');
    line([5.5 5.5],[-0.5 20.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    line([15.5 15.5],[-0.5 20.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    
    line([-0.5 20.5],[10.5 10.5],'LineWidth',3,'Color','k')
    line([-0.5 20.5],[5.5 5.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    line([-0.5 20.5],[15.5 15.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
    xlabel(sprintf(['Prob.              Mag.              Prob.             Mag.\n' ... 
        'Left                                      Right']))
    set(get(gca,'XLabel'),'FontSize',14);
    ylabel(sprintf(['Right                                      Left\n' ...
        'Mag.             Prob.               Mag.             Prob.']))
    %ylabel('Left                              Right')
    set(get(gca,'YLabel'),'FontSize',14);
end

%% construct sliding RSA matrices for many different points in time
for i = 1:3
    for ttt = 1:800
        dat = squeeze(unitmats{i}(:,:,ttt)');
        ok = ~all(dat==0); %exclude any cells that don't spike
        RSAmatt{i}(:,:,ttt) = corrcoef(normalise(dat(:,ok),1)');
    end
end

%% design matrix to be regressed onto RSA matrix (each area separately)

m1 = repmat(-2:2,5,1);
m2 = repmat([-2:2]',1,5);

%build regressors
r1 = ones(20,20); rt{1} = 'mean';
r2 = eye(20);  rt{2} = 'mean diagonal';
r3 = [zeros(10) eye(10); eye(10) zeros(10)]; rt{3} = 'stimulus identity';
r4 = repmat(m1.*m2,4,4); rt{4} = 'attended value'; 
%alternative version of attended value regressor, for Fig. S11: 
%r4 = -repmat(abs(m1-m2),4,4); r4 = r4 - mean(r4(:)); 
r5 = repmat((m1<0&m2<0)|(m1>0&m2>0),4,4) - repmat((m1>0&m2<0)|(m1<0&m2>0),4,4); rt{5} = 'accept/reject';
r6 = [ones(10,10) zeros(10,10); zeros(10,10) ones(10,10)]; rt{6} = 'spatial attention';
r7  = (r6>0).*r4; rt{7} = 'left/right value';

%% perform the regression
for i = 1:3
    RSAdat(:,i) = RSAmat{i}(:);
end

%build design matrix:
dm = [r1(:) r2(:) r3(:) r4(:) r5(:) r6(:) r7(:)];
dm(:,3:7) = dm(:,3:7) ./ repmat(max(abs(dm(:,3:7))),[size(dm,1) 1]);
cmat = eye(size(dm,2));
%cmat(8,:) = cmat(4,:)+cmat(7,:);rt{8} = rt{4};

%run Ordinary Least Squares regression
[c,v,t] = ols(RSAdat,dm,cmat);

%reshape to nRegressors*20*20, for plotting below
dmr = (reshape(dm',size(dm,2),20,20)); 

%%  design matrix (all areas combined), for F-tests between regions

dmS = zeros(1200,21);
dmS(1:400,1:7) = dm;
dmS(401:800,8:14) = dm;
dmS(801:1200,15:21) = dm;
for i = 1:7
   Ftests{i} = [i i+7]; 
end

%% perform F-test
cmatS = [eye(14) -[eye(7); eye(7)]];
[cS,vS,tS,FS,Fdof] = ols(RSAdat(:),dmS,cmatS,Ftests');

%% run permutation tests, to calculate null distribution for both T-stats and F-stats

nP = 100; %number of permutations - in the paper this was set to 10000.

for p = 1:nP
    pp = randperm(20); %randomly permute the 20 conditions
    
    %recompute the RSA matrices for each region:
    for i = 1:3
        ok = ~all(unitmat{i}==0); %exclude any cells that don't spike
        RSAmatp{i} = corrcoef(normalise(unitmat{i}(pp,ok),1)');
        RSAdatp(:,i) = RSAmatp{i}(:);
    end
    
    %calculate the test statistics in this permuted RSA matrix:
    [cp(p,:,:),vp(p,:,:),tp(p,:,:)] = ols(RSAdatp,dm,cmat);
    [~,~,~,FSp(p,:)] = ols(RSAdatp(:),dmS,cmatS,Ftests');
    if mod(p,100)==0
        fprintf('%0.0f of %0.0f permutations complete',p,nP);
    end
end

%% repeat regression on sliding analysis, and calculate sliding CPD

timeind = -199:600;

for k=1:3
    for j = 1:800
        RSAdatt(:,j,k) = squash(RSAmatt{k}(:,:,j));
    end
end

[cpd_out,RX] = cpd(RSAdatt(:,:),dm);cpd_out = cpd_out.*100;

cpds = reshape(cpd_out,size(dm,2),800,3);

%% recompute with different permutations of the residuals

Xb = dm*pinv(dm)*RSAdatt(:,:); %original Xb from GLM
for p = 1:nP
    pmat = randperm(400);
    new_data = Xb + RX(pmat,:); %new simulated data, with the residualals from the model fit
    
    [cpd_out_perm] = cpd(new_data,dm);cpd_out_perm = cpd_out_perm.*100; %re-estimate CPDout
    cpds_perm(:,:,:,p) = reshape(cpd_out_perm,size(dm,2),800,3);
    
    %estimate t75s:
    for i = 3:7 %regressor
        for r = 1:3 %region
            [t75max_perm(i,r,p),tmax_perm(i,r,p),maxval_perm(i,r,p)] = calculate_t75max(squeeze(cpds_perm(i,:,r,p)),timeind);
        end
    end
end

%% plot these - supplementary figure 5C

figure; 
set(gcf,'Position',[1000         172         482        1166]);
stack_order = [6 1; 6 3; 6 2; 4 3; 5 3; 5 1; 3 3;  7 2; 7 1;]; 
for i = 1:length(stack_order)
    t75_perm_stack(:,i) = t75max_perm(stack_order(i,1), stack_order(i,2), :); 
    stack_label{i} = sprintf('%s: %s',area_name{stack_order(i,2)},rt{stack_order(i,1)});
end
   
hold on;
for i = 1:size(t75_perm_stack,2)
    plot(i+randn(nP,1)*0.05,t75_perm_stack(:,i),'.','Color',[0.75 0.75 0.75]);
end

bph = boxplot(t75_perm_stack,stack_label,'symbol','','Notch','on');
set(bph,{'linew'},{2})
box off


set(gca,'XTickLabelRotation',90,'YAxisLocation','Right','YTickLabelRotation',90);
ylabel('t_7_5 (ms after cue 1 onset)');
ylim([0 550]);
tidyfig;

%%

%% alternative permutation test, to calculate null distribution of CPD - here we permute the design matrix rather than the data
%{
for p = 1:nP
    pp = randperm(20); %randomly permute the 20 conditions

    dmr_perm = dmr(:,pp,pp); 
    dm_perm = dmr_perm(:,:)';
    cpd_out_perm = cpd(RSAdatt(:,:),dm_perm).*100;

    cpds_perm(:,:,:,p) = reshape(cpd_out_perm,size(dm_perm,2),800,3);

end

max_perm = squeeze(max(cpds_perm,[],2));
cpd_threshold = prctile(max_perm,99.5,3); %p<0.01 threshold (corrected for multiple comparisons)
%}


%% plot results for different regressors

axorder = [3 2 1];
figure(2);clf;set(gcf,'Position',[395  145 1344 607]);
ebars = squeeze(std(tp,[],1));

for i = 3:7
    %subplot(5,3,(i-3)*3+1);
    subplot(3,5,i-2);
    imagesc(squeeze(dmr(i,:,:)));title(rt{i});tidyfig;
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    switch i
        case 3
            caxis([-1.5 1.5]);
        case 4
            caxis([-1.5 1.5]);
        case 5
            caxis([-1.5 1.5]);
        case 6
            caxis([-1.5 1.5]);
        case 7
            caxis([-1.5 1.5]);
    end

    colormap hot;freezeColors;
    %ylabel(rt{i});
    axis square
    line([10.5 10.5],[-0.5 20.5],'LineWidth',1,'Color','k');
    line([-0.5 20.5],[10.5 10.5],'LineWidth',1,'Color','k');
    

    if (i==4|i==5)
    line([5.5 5.5],[-0.5 20.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    line([15.5 15.5],[-0.5 20.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    
    line([-0.5 20.5],[5.5 5.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    line([-0.5 20.5],[15.5 15.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    end
    
    %subplot(5,3,(i-3)*3+2);
    subplot(3,5,i+3);
    if i==3
        bb=bar([t(i,axorder); zeros(size(t(i,axorder)))],0.5,'LineWidth',2);
        xlim([0.5 1.5]); box off; set(gca,'Box','off','XTickLabel',[],'LineWidth',2);
        l = legend(area_name(axorder));set(l,'Box','off');
    else
        bb=bar([t(i,axorder); zeros(size(t(i,axorder)))],0.5,'LineWidth',2);
        xlim([0.5 1.5]); box off; set(gca,'Box','off','XTickLabel',[],'LineWidth',2);
    end
    ylabel(sprintf('T-statistic\n(%s)',rt{i}));
    set(gca,'FontSize',13);set(get(gca,'YLabel'),'FontSize',13);
    axis square
    colormap parula;freezeColors
    
    %subplot(5,3,(i-3)*3+3);hold on;
    subplot(3,5,i+8);hold on;
    ccvec = [237 228 22;62 181 156;53 45 132]/255;
    for rr = 1:3
        plot(timeind, squeeze(cpds(i,:,axorder(rr))), 'LineWidth',2,'Color',ccvec(axorder(rr),:));
        tmp = get(gca,'YLim'); ymax = tmp(2); clear tmp;
        %significant_timepoints = cpds(i,:,axorder(rr))>cpd_threshold(i,axorder(rr));
        %plot(timeind(significant_timepoints),ymax)
    end
    xlabel('Time after stimulus onset (ms)');
    ylabel(sprintf('CPD (%%)\n(%s)',rt{i}));
    set(gca,'FontSize',13);set(get(gca,'YLabel'),'FontSize',13);
    axis square
    xlim([-200 500]);box off
end

%% plot each of the different regressors
%{
figure;
for i = 2:7
    %subplot(5,3,(i-3)*3+1);
    subplot(1,6,i-1);
    imagesc(squeeze(dmr(i,:,:)));title(rt{i});tidyfig;
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    switch i
        case 2
            caxis([-max(abs(dmr(i,:))) max(abs(dmr(i,:)))]);
        case 3
            caxis([-max(abs(dmr(i,:))) max(abs(dmr(i,:)))]);
        case 4
            caxis([-max(abs(dmr(i,:))) max(abs(dmr(i,:)))]*1.5);
        case 5
            caxis([-max(abs(dmr(i,:))) max(abs(dmr(i,:)))]*1.5);
        case 6
            caxis([-max(abs(dmr(i,:))) max(abs(dmr(i,:)))]*1.5);
        case 7
            caxis([-max(abs(dmr(i,:))) max(abs(dmr(i,:)))]*1.5);
    end

    colormap hot;freezeColors;
    %ylabel(rt{i});
    axis square
    line([10.5 10.5],[-0.5 20.5],'LineWidth',1,'Color','k');
    line([-0.5 20.5],[10.5 10.5],'LineWidth',1,'Color','k');
    

    if (i==4|i==5)
    line([5.5 5.5],[-0.5 20.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    line([15.5 15.5],[-0.5 20.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    
    line([-0.5 20.5],[5.5 5.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    line([-0.5 20.5],[15.5 15.5],'LineWidth',1,'Color',[0.5 0.5 0.5]);
    end
end
%}

%% supplementary figure - relative timing in different regions

figure;
set(gcf,'Position',[1386 886 1074 420]);tidyfig;
for i = 1:3; 
    subplot(1,3,i); 
    plot(timeind,100*cpds(3:7,:,i)./repmat(max(cpds(3:7,:,i),[],2),1,800),'LineWidth',2); 
    xlim([-100 500]);ylim([0 100]);
    xlabel('Time after Cue 1 onset (ms)');
    if i==1
        ylabel('Coefficient of partial determination (%)');
        l =legend(rt(3:7));
        set(l,'Location','NorthWest','Box','off');
    end
    box off;
    title(area_name{i});
    tidyfig(gcf,12);
end

%% statistics on timing - what is time to maximum/half-maximum of CPD?

for i = 3:7 %regressor
    for r = 1:3 %region
        [t75max(i,r),tmax(i,r),maxval(i,r)] = calculate_t75max(squeeze(cpds(i,:,r)),timeind);
    end
end

tmax(maxval<10) = nan;
fprintf('======\n')
for i = [6 3 4 5 7]
    fprintf('Regressor: %s\n',rt{i});
    for r = 1:3
        if maxval(i,r)>10
            fprintf('%s Time of maximum: %0.0f ms\tPeak CPD: %0.2f%%\n',area_name{r},t75max(i,r),maxval(i,r));
        end
    end
    fprintf('\n');
end
fprintf('======\n')

%% show stats from permutation tests (region comparisons)

fprintf('\n***************\n')
R1 = 3; R2 = 2; reg = 3; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);
R1 = 3; R2 = 1; reg = 3; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);

R1 = 3; R2 = 2; reg = 4; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);
R1 = 3; R2 = 1; reg = 4; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);

R1 = 1; R2 = 3; reg = 5; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);
R1 = 1; R2 = 2; reg = 5; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);

R1 = 2; R2 = 1; reg = 6; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);
R1 = 2; R2 = 3; reg = 6; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);

R1 = 1; R2 = 2; reg = 7; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);
R1 = 1; R2 = 3; reg = 7; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);
R1 = 2; R2 = 3; reg = 7; pval = 1-mean((t(reg,R1)-t(reg,R2))>(tp(:,reg,R1)-tp(:,reg,R2)));
fprintf('Comparison for %s > %s for %s: p=%0.4f\n',area_name{R1},area_name{R2},...
    rt{reg},pval);


%% show stats from permutation tests (significant difference from 0

fprintf('\n***************\n')
sigmat = 1-squeeze(mean(permute(repmat(t,1,1,nP),[3 1 2])>tp,1)); %signficance for all
for j = [3 1 2]
    fprintf('Significance for %s\n',area_name{j})
    for i = 3:7
        if i ==5
            fprintf('%s: \t\tT=%0.4f, p=%0.4f\n',rt{i},t(i,j),sigmat(i,j));
        else
            fprintf('%s: \tT=%0.4f, p=%0.4f\n',rt{i},t(i,j),sigmat(i,j));
        end
    end
    fprintf('\n')
end

%%

sigFstats = 1-mean(repmat(FS',nP,1)>FSp);
for i = 3:7
    fprintf('Significance for F-stat, regressor %s: F = %0.2f, p=%0.4f\n',rt{i},FS(i),sigFstats(i))
end


%% supplementary video

drawmov = 0;

if drawmov
    fc = 1;
    myMovie = VideoWriter('RSA_movie_all_coregressors.avi');
    myMovie.FrameRate = 8;
    myMovie.Quality = 100;
    myMovie.open;
    
    figure(2);clf;
    set(gcf,'Position',[100  100        1738        1026]);
    colormap hot
    for ttt = 5:5:700
        
        for i = 1:3
            subplot(2,3,axorder(i));
            dat = squeeze(unitmats{i}(:,:,ttt)');
            ok = ~all(dat==0); %exclude any cells that don't spike
            tmp = corrcoef(normalise(dat(:,ok),1)');
            imagesc(tmp);caxis([-0.3 0.3]);
            title(area_name{i});
            axis square; tidyfig;
            set(gca,'YTickLabel',[1:5],'YTick',1:20,'XTickLabel',[1:5],'XTick',1:20,'FontSize',10);
            line([10.5 10.5],[-0.5 20.5],'LineWidth',3,'Color','k');
            line([5.5 5.5],[-0.5 20.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
            line([15.5 15.5],[-0.5 20.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
            
            line([-0.5 20.5],[10.5 10.5],'LineWidth',3,'Color','k')
            line([-0.5 20.5],[5.5 5.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
            line([-0.5 20.5],[15.5 15.5],'LineWidth',1.5,'Color',[0.5 0.5 0.5]);
            xlabel(sprintf(['Prob.              Mag.              Prob.             Mag.\n' ...
                'Left                                      Right']))
            set(get(gca,'XLabel'),'FontSize',14);
            ylabel(sprintf(['Right                                      Left\n' ...
                'Mag.             Prob.               Mag.             Prob.']))
            %ylabel('Left                              Right')
            set(get(gca,'YLabel'),'FontSize',14);
        end
        if ttt<find(timeind==0)
            suptitle(sprintf('Time = %0.0fms: pre-stimulus',timeind(ttt)));
        elseif ttt<find(timeind==300)
            tit = suptitle(sprintf('Time = %0.0fms: stimulus on',timeind(ttt)));
            set(tit,'Color','red');
        else
            suptitle(sprintf('Time = %0.0fms: stimulus off',timeind(ttt)));
        end
        
        plotcount = 1;
        for p = [6 3 4 7 5]
            subplot(4,5,plotcount+10);imagesc(squeeze(dmr(p,:,:)));title(rt{p});
            if p>3
                caxis([-max(abs(dmr(p,:))) max(abs(dmr(p,:)))]*1.5);
            else
                caxis([-max(abs(dmr(p,:))) max(abs(dmr(p,:)))]);
            end
            
            subplot(4,5,plotcount+15);
            hold on;
            for rr = [1:3]
                plot(timeind(1:ttt),squeeze(cpds(p,1:ttt,axorder(rr))),'LineWidth',2,'Color',ccvec(axorder(rr),:));
            end
            xlim([-200 500]);ylim([0 90]);
            box off;
            xlabel('Time after stimulus onset (ms)');
            if p == 3
                set(gca,'YLim',[0 60]);
            elseif p ==4
                set(gca,'YLim',[0 15]);
            elseif p ==5
                set(gca,'YLim',[0 40]);
                l = legend(area_name([3 2 1]));
                set(l,'Box','off','FontSize',11,'Location','NorthWest');
            elseif p ==6
                ylabel('Coefficient of Partial Determination (%)');
                set(gca,'YLim',[0 100]);
            elseif p == 7
                set(gca,'YLim',[0 30]);
            end
            plotcount = plotcount + 1;
        end
        
        drawnow;ff = getframe(gcf);writeVideo(myMovie,ff);
        clf
    end
    myMovie.close;
end