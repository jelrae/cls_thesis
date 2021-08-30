%% Set up paths
addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\cls_thesis\gc_hierarchies\')
addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\cls_thesis\helper_functions\')
%addpath('D:/MATLAB/mvgc_v1.0')
startup
addpath('C:\Users\Jordan\Documents\cls_thesis\matlab\fieldtrip-20210411')
ft_defaults

format short;
clear all;
close all;clc;
ptic('starting\n')

load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_all_AttIn.mat')
monkey = 'kurt';

%% Parameters

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

model_order    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
% morder    = 'AIC'; 
momax     = 30;     % maximum model order for model order estimation (<= 30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1000;    % sample rate (Hz)

% fres      = fs/nobs; % for FFT defined as fs/N samples
fres      = 1000;
fnq       = fs/2;

seed      = 0;      % random seed (0 for unseeded)
verb = false; 

%% Variables to help with counting and data locations.  Not sure if needed
sizes = [];
errors = [];

%% Storage for all GC

gc_forward = {};
gc_backward = {};

%% Identify the regions
fig_6_ROIS;

%% Do GC
% Get all the regions
roi_cfg = [];
roi_cfg.channel = all_channels;
all_regions = ft_selectdata(roi_cfg, all_AttIn);
all_regions = all_regions.trial;
all_regions = cat(3,all_regions{:});
ntrials   = size(all_regions,3);     % no. trials
nobs      = size(all_regions,2);   % no. obs per trial

figure(1);
regions(6) = [];
% for each region
for region = 1 : length(regions)
    %Reset the gc storing arrays
    gc_one    = [];
    gc_two    = [];
    % Get the regions 
    roi_cfg = [];
    roi_cfg.channel = regions{region};
    current_region = ft_selectdata(roi_cfg, all_AttIn);
    current_region = current_region.trial;
    current_region = cat(3,current_region{:});
    % Cycle over the channels for each of the two (all regions and current)
    % for testing purposes we are just going to do one of each and try to
    % plot
    region_comp = cat(1, current_region(1,:,:), all_regions(1,:,:));
    
    %% Start the calcs of granger cause
    % VAR - channel combinations 
    [aic,bic,moaic,mobic] = tsdata_to_infocrit(region_comp,momax,icregmode, verb);

    if strcmpi(model_order,'AIC')
        morder = moaic;
        fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(model_order,'BIC')
        morder = mobic;
        fprintf('\nusing BIC best model order = %d\n',morder);
    else
        fprintf('\nusing specified model order = %d\n',morder);
        morder = model_order;
    end
    
    ptic('\n*** tsdata_to_var... ');
    [A,SIG] = tsdata_to_var(region_comp,morder,regmode);

    ptoc;
    % Check for failed regression
    assert(~isbad(A),'VAR estimation failed');

    % Autocovariance calculation
    ptic('*** var_to_autocov... ');
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    ptoc;
    
    try
        var_info(info,true); % report results (and bail out on error)

        ptic('\n*** autocov_to_spwcgc... ');
        f= autocov_to_spwcgc(G,fres);
        ptoc;
        % Check for failed spectral GC calculation
        assert(~isbad(f,false),'spectral GC calculation failed');

        if ~isreal(f)
            fprintf('\nCurrent channel combination is complex: %d, %d, didnt work!\n', chan1, chan4)
        end

        % Only collect f of size maz_len_f
        gc_one(end+1,:) = squeeze(f(1,2,:)); 
        gc_two(end+1,:) = squeeze(f(2,1,:));

        sizes(end+1) = size(f,3);

    catch % for intstances with unstable VAR root

        fprintf('\nCurrent channel combination is of region %d: %d, %d, didnt work!\n', region, chan1, chan4)
        errors(end+1,:) = [region 1 1];
        continue   
    end
        
    gc_forward{end+1} = gc_one;
    gc_backward{end+1} = gc_two;
    
    %% start plotting
    
    % average the results and plot
    gc_1_ave = mean(gc_one, 1);
    gc_2_ave = mean(gc_two, 1);
    % f_all = vertcat(gc_1_ave, gc_2_ave);
    % gc_plot(f_all, fs);

    %% do some plotting
    x_range = (1:1:length(gc_1_ave))./length(gc_1_ave);
    x_range = x_range.*fnq;
    blue=[.3 .8 .9];red=[1 .1 .1];
    
    subplot(length(regions), 4, (4*region) - 3)
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%     set(gca,'FontSize',24,'LineWidth',5,'TickLength',[0.02 0.02])
%     set(gca,'box','off');legend('FF','FB');legend boxoff;
    legend('FF','FB');legend boxoff;
    xlim([0 10]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);title('V1-V4');set(gca,'Layer','top');
    ylabel('Granger causality');xlabel('Frequency (Hz)');
    
    subplot(length(regions), 4, (4*region) - 2)
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%     set(gca,'FontSize',24,'LineWidth',5,'TickLength',[0.02 0.02])
%     set(gca,'box','off');legend('FF','FB');legend boxoff;
    legend('FF','FB');legend boxoff;
    xlim([10 30]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);title('V1-V4');set(gca,'Layer','top');
    ylabel('Granger causality');xlabel('Frequency (Hz)');
    
    subplot(length(regions), 4, (4*region) - 1)
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%     set(gca,'FontSize',24,'LineWidth',5,'TickLength',[0.02 0.02])
%     set(gca,'box','off');legend('FF','FB');legend boxoff;
    legend('FF','FB');legend boxoff;
    xlim([30 100]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);title('V1-V4');set(gca,'Layer','top');
    ylabel('Granger causality');xlabel('Frequency (Hz)');
    
    subplot(length(regions), 4, (region*4))
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%     set(gca,'FontSize',24,'LineWidth',5,'TickLength',[0.02 0.02])
%     set(gca,'box','off');legend('FF','FB');legend boxoff;
    legend('FF','FB');legend boxoff;
    xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);title('V1-V4');set(gca,'Layer','top');
    ylabel('Granger causality');xlabel('Frequency (Hz)');
    
end