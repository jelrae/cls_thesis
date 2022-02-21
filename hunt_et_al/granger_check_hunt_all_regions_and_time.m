%% Initialization

clear all
clc

addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
set_paths_roi;

ft_defaults

%% Load in basic things
% addpath('/home/jordan/neuro_thesis/data/hunt_et_al/scripts/merged_F')
% load('monkey_sessions_F.mat')
% monkey = 'F';
% sessions = sessions_F;

% addpath('/home/jordan/neuro_thesis/data/hunt_et_al/scripts/merged_M')
load('monkey_sessions_M.mat')
monkey = 'M';
sessions = sessions_M;

save_dir = 'power_check';

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
fres      = 1000;
fnq       = fs/2;

seed      = 0;      % random seed (0 for unseeded)
verb = false; 

%% Preprocessing cfg
pre_cfg = [];
pre_cfg.demean = 'yes';
pre_cfg.detrend = 'yes';
pre_cfg.dftfilter = 'yes';

%% Variables to help with counting and data locations.  Not sure if needed
sizes = [];
errors = [];
bad_combinations = [];

%% Region checking points

region_list = {'DLPFC', 'ACC', 'OFC'};
% combination_list = {{'ACC','DLPFC'}, {'ACC', 'OFC'}, {'OFC', 'DLPFC'}};
region_list = {'OFC','ACC'};
num_regions = length(combination_list);

num_sessions = length(sessions);
% skip_sessions_F = [13,23,24,30,31,33,34];
% skip_sessions_M = [10,12,22,36,41,44,46,48];
include_sessions = [07,11,12,14,15,21,23,27,30,41,42];
plotting_colors = distinguishable_colors(length(include_sessions));
% plotting_colors = distinguishable_colors(num_sessions);
% plotting_colors = distinguishable_colors(length(skip_sessions));
color_num = 1;

plot_titles = {};
plot_titles(end+1) = {'Pre Fixation'};
plot_titles(end+1) = {'During Fixation'};
plot_titles(end+1) = {'First Image'};
plot_titles(end+1) = {'Second Image'};
plot_titles(end+1) = {'Pre Choice'};
plot_titles(end+1) = {'Pre Feedback'};

for sn = 1:length(sessions)
    session_number = sessions{sn};
    disp(session_number)
    s_n = str2num(session_number(2:end));
    if ismember(s_n,include_sessions)
%     if ismember(s_n,skip_sessions)
        file_name = sprintf('%s_%s_data.mat',monkey, session_number);
        base_name = sprintf('%s_%s_data',monkey, session_number);
        load(file_name);

        data_sets = {};
        data_sets(end+1) = {append(base_name,'_first_fix')};
        data_sets(end+1) = {append(base_name,'_finish_fix')};
        data_sets(end+1) = {append(base_name,'_1')};
        data_sets(end+1) = {append(base_name,'_2')};
        data_sets(end+1) = {append(base_name,'_pre_choice_300')};
        data_sets(end+1) = {append(base_name,'_pre_feedback')};

        for time = 1:length(data_sets)
            data_name = data_sets{time};
            for c_n = 1:length(combination_list)
                gc_one    = [];
                gc_two    = [];
                interest_regions = combination_list{c_n};
                interest_region_1 = interest_regions{1};
                interest_region_2 = interest_regions{2};
                % select the channels
                eval(sprintf('channels_1 = %s_%s_channels;',base_name, interest_region_1));
                eval(sprintf('channels_2 = %s_%s_channels;',base_name, interest_region_2));
                if size(channels_1,1) > 0 && size(channels_2,1) > 0
                    eval(sprintf('%s = ft_preprocessing(pre_cfg, %s);', data_name, data_name));
                    % get the things for the channels
                    for chan1 = channels_1
                        roi_cfg = [];
                        roi_cfg.channel = chan1;
                        eval(sprintf('region1 = ft_selectdata(roi_cfg, %s);', data_name))
                        region1 = region1.trial;
                        region1 = cat(3,region1{:});
                        for chan2 = channels_2
                            roi_cfg = [];
                            roi_cfg.channel = chan2;
                            eval(sprintf('region2 = ft_selectdata(roi_cfg, %s);', data_name))
                            region2 = region2.trial;
                            region2 = cat(3,region2{:});
                            region_comp = cat(1, region1, region2);
                            
                            fprintf('\nCurrent channel combination for regions %s %s is: %s, %s\n', interest_region_1, interest_region_2, chan1{1}, chan2{1})
                            
                            try
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
                                fprintf('\nCurrent channel combination is of regions %s %s: %d, %d, was not pos def for chol!\n', region_names(r1), region_names(r2), chan1, chan2)
                                bad_combinations(end+1,:) = [r1 r2 chan1 chan2];
                            end

                            try
                %                 var_info(info,true); % report results (and bail out on error)

                                f= autocov_to_spwcgc(G,fres);
                                % Check for failed spectral GC calculation
                                assert(~isbad(f,false),'spectral GC calculation failed');

                                if ~isreal(f)
                                    fprintf('\nCurrent channel combination is complex %s, %s: %d, %d, was not pos def for chol!\n', region_names(r1), region_names(r2), chan1, chan2)
                                end

                                % Updated from source code first dim is to, second is
                                gc_one(end+1,:) = squeeze(f(2,1,:)); 
                                gc_two(end+1,:) = squeeze(f(1,2,:));

                                sizes(end+1) = size(f,3);

                            catch % for intstances with unstable VAR root

                                fprintf('\nCurrent channel combination of regions %s, %s: %d, %d, didnt work!\n', region_names(r1), region_names(r2), chan1, chan2)
                                errors(end+1,:) = [r1 r2 chan1 chan2];
                                continue   
                            end
                        end
                    end
                end
                gc_1_ave = mean(gc_one, 1);
                gc_2_ave = mean(gc_two, 1);
                x_range = (1:1:length(gc_1_ave))./length(gc_1_ave);
                x_range = x_range.*fnq;
                % Feed Forward
                hold on
                ff_pn = (2*num_regions*(time-1)+(2*c_n)-1);
                subplot(length(data_sets),length(combination_list)*2, ff_pn)
                %loglog
                plot(x_range,gc_2_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('%s - %s during %s', interest_region_1, interest_region_2, plot_titles{time}))
                xlim([0 140]);
                ylabel('Granger Causality')
                xlabel('Frequency (Hz)')
                % Feedback
                hold on
                subplot(length(data_sets),length(combination_list)*2, ff_pn+1)
                %loglog
                plot(x_range,gc_1_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('%s - %s during %s', interest_region_2, interest_region_1, plot_titles{time}))
                xlim([0 140]);z
                ylabel('Granger Causality')
                xlabel('Frequency (Hz)')
%                     set(gca,'xscale','log')
%                     set(gca,'yscale','log')                
%                 legend
            end 
        end
        color_num = color_num + 1;
    end
end
subplot(length(data_sets),length(combination_list)*2, 31)
legend
subplot(length(data_sets),length(combination_list)*2, 33)
legend
subplot(length(data_sets),length(combination_list)*2, 35)
legend
sgtitle(sprintf('Granger Causality Analysis %s', monkey))
