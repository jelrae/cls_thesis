% adding imports.

% addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
% set_paths_roi;

addpath('/home/jordan/common/matlab/fieldtrip-20210411')

clear all
clc

startup
ft_defaults

load('monkey_sessions_F.mat')

monkey = 'F';

all_regions = {'ACC', 'DLPFC', 'OFC'};

pre_cfg = [];
pre_cfg.demean = 'yes';
% pre_cfg.detrend = 'yes';
pre_cfg.dftfilter = 'yes';
pre_cfg.dftfreq = [50];

%% Parameters

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default) OG LWR

model_order    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
% morder    = 'AIC'; 
momax     = 30;     % maximum model order for model order estimation (<= 30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1000;    % sample rate (Hz)

% fres      = fs/nobs; % for FFT defined as fs/N samples
% I think these might be wrong?
fres      = 1000;
fnq       = fs/2;

seed      = 0;      % random seed (0 for unseeded)
verb = false; 

% change for monkey
sessions = sessions_F;

figure(1);

for sn = 1:length(sessions)
    session_number = sessions{sn};
    disp(session_number)
    file_name = sprintf('%s_%s_data.mat',monkey, session_number);
    base_name = sprintf('%s_%s_data',monkey, session_number);
    data_name = append(base_name,'_1');
    load(file_name);
    
    %reset the session regions
    session_regions = {};
    
    % figure out which regions there are data for
    for region = all_regions
        region_name = append(base_name, '_', region{1}, '_channels');
        eval(sprintf('region_channels = %s;',region_name));
        if ~isempty(region_channels)
            session_regions(end+1) = {region{1}};
        end
    end
    % now check if we have at least 2 regions with things in it
    if size(session_regions,2) >= 2
        a_cfg = [];
        a_cfg.method = 'summary';
        a_cfg.channel = {'all'};
%         eval(sprintf('%s = ft_rejectvisual(a_cfg, %s);', data_name, data_name));
        eval(sprintf('%s = ft_preprocessing(pre_cfg, %s);', data_name, data_name));
%         data = ft_rejectvisual(a_cfg,data);
        for r1 = session_regions
            roi_cfg = [];
            eval(sprintf('roi_cfg.channel = %s_%s_channels;', base_name, r1{1}));
            region1 = ft_selectdata(roi_cfg, eval(data_name));
            region1 = region1.trial;
            region1 = cat(3,region1{:});
            for r2 = session_regions
                if ~strcmp(r1{1}, r2{1})
                    roi_cfg = [];
                    gc_one    = [];
                    gc_two    = [];
                    eval(sprintf('roi_cfg.channel = %s_%s_channels;', base_name, r2{1}));
                    region2 = ft_selectdata(roi_cfg, eval(data_name));
                    region2 = region2.trial;
                    region2 = cat(3,region2{:});
                    % add the GC loops
                    for chan1 = 1 : size(region1, 1);
                        for chan2 = 1 : size(region2, 1); 
                            try
                                fprintf('\nCurrent channel combination for regions %s %s is: %d, %d\n', r1{1}, r2{1}, chan1, chan2)
                                region_comp = cat(1, region1(chan1,:,:), region2(chan2,:,:));

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

                                [A,SIG] = tsdata_to_var(region_comp,morder,regmode);

                                % Check for failed regression
                                assert(~isbad(A),'VAR estimation failed');

                                % Autocovariance calculation
                                [G,info] = var_to_autocov(A,SIG,acmaxlags);
                            catch
                                disp('\nCurrent channel combination is bad pre gc calc')
                            end
                            
                            try
%                                 var_info(info,true); % report results (and bail out on error)

                                f= autocov_to_spwcgc(G,fres);
                                % Check for failed spectral GC calculation
                                assert(~isbad(f,false),'spectral GC calculation failed');

                                if ~isreal(f)
                                    fprintf('\nCurrent channel combination is complex %s, %s: %d, %d, was not pos def for chol!\n', region_names(r1), region_names(r2), chan1, chan2)
                                end

                                % Updated from source code first dim is to, second is
                                % from
                                gc_one(end+1,:) = squeeze(f(2,1,:)); 
                                gc_two(end+1,:) = squeeze(f(1,2,:));
                            catch % for intstances with unstable VAR root
                                disp('\nCurrent channel combination is bad gc calc')
%                                 fprintf('\nCurrent channel combination of regions %s, %s: %d, %d, didnt work!\n', region_names(r1), region_names(r2), chan1, chan2)
%                                 errors(end+1,:) = [r1 r2 chan1 chan2];
                                continue   
                            end 
                        end
                    end
                    % Done with region v.s. region compare
                    % Do test plotting
                    gc_1_ave = mean(gc_one,1);
                    gc_2_ave = mean(gc_two,1);
                    
                    x_range = (1:1:length(gc_1_ave))./length(gc_1_ave);
                    x_range = x_range.*fnq;
                    blue=[.3 .8 .9];red=[1 .1 .1];

                    if strcmp(r1{1}, 'ACC') && strcmp(r2{1}, 'DLPFC')
                        subplot(3,3,2)
                        plot(x_range,gc_2_ave,'Color',blue,'LineWidth',1);hold on;
                        plot(x_range,gc_1_ave,'Color',red,'LineWidth',1);hold on;
                        xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                         ylim([0 zgc]);                        
                        title('ACC v.s. DLPFC');
                        set(gca,'Layer','top');
                    elseif strcmp(r1{1}, 'ACC') && strcmp(r2{1}, 'OFC')
                        subplot(3,3,3)
                        plot(x_range,gc_2_ave,'Color',blue,'LineWidth',1);hold on;
                        plot(x_range,gc_1_ave,'Color',red,'LineWidth',1);hold on;
                        xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                         ylim([0 zgc]);                        
                        title('ACC v.s. OFC');
                        set(gca,'Layer','top');
                    elseif strcmp(r1{1}, 'DLPFC') && strcmp(r2{1}, 'OFC')
                        subplot(3,3,6)
                        plot(x_range,gc_2_ave,'Color',blue,'LineWidth',1);hold on;
                        plot(x_range,gc_1_ave,'Color',red,'LineWidth',1);hold on;
                        xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                         ylim([0 zgc]);                        
                        title('DLPFC v.s. OFC');
                        set(gca,'Layer','top');
                    elseif strcmp(r1{1}, 'DLPFC') && strcmp(r2{1}, 'ACC')
                        subplot(3,3,4)
                        plot(x_range,gc_2_ave,'Color',blue,'LineWidth',1);hold on;
                        plot(x_range,gc_1_ave,'Color',red,'LineWidth',1);hold on;
                        xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                         ylim([0 zgc]);                        
                        title('DLPFC v.s. ACC');
                        set(gca,'Layer','top');
                    elseif strcmp(r1{1}, 'OFC') && strcmp(r2{1}, 'ACC')
                        subplot(3,3,7)
                        plot(x_range,gc_2_ave,'Color',blue,'LineWidth',1);hold on;
                        plot(x_range,gc_1_ave,'Color',red,'LineWidth',1);hold on;
                        xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                         ylim([0 zgc]);
                        title('OFC v.s. ACC');
                        set(gca,'Layer','top');
                    elseif strcmp(r1{1}, 'OFC') && strcmp(r2{1}, 'DLPFC')
                        subplot(3,3,8)
                        plot(x_range,gc_2_ave,'Color',blue,'LineWidth',1);hold on;
                        plot(x_range,gc_1_ave,'Color',red,'LineWidth',1);hold on;
                        xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                         ylim([0 zgc]);
                        title('OFC v.s. DLPFC');
                        set(gca,'Layer','top');
                    else 
                        subplot(3,3,1)
                        plot(x_range,gc_2_ave,'Color',blue,'LineWidth',1);hold on;
                        plot(x_range,gc_1_ave,'Color',red,'LineWidth',1);hold on;
                        xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                         ylim([0 zgc]);
                        title('What are these?');
                        set(gca,'Layer','top');
                    end
%                     disp(session_number)
%                     plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
%                     plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%                     legend('FF','FB');legend boxoff;
%                     xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%                     ylim([0 zgc]);
%                     title('GC Plot F');
%                     set(gca,'Layer','top');
%                      disp('--------------------');
                end
            end
        end
    end
end

han = axes(figure(1), 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han,'Granger Causality');xlabel(han, 'Frequency (Hz)');
title(han,sprintf('%s Granger Causeality \n', monkey));




% if ~isempty(region1_channels)
%             for region2 = all_regions
%                 if region1 ~= region2
%                     region2_channels = append(base_name, '_', region2, '_channels');
%                     if ~isempty(region1_channels)
% 
%                     end
%                 end
%             end
%         end