%% Initialization

clear all
clc

addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
set_paths_roi;

ft_defaults

%% Load in basic things
% addpath('/home/jordan/neuro_thesis/data/hunt_et_al/scripts/merged_F')
load('monkey_sessions_F.mat')
monkey = 'F';
sessions = sessions_F;

% addpath('/home/jordan/neuro_thesis/data/hunt_et_al/scripts/merged_M')
% load('monkey_sessions_M.mat')
% monkey = 'M';
% sessions = sessions_M;

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
fres      =  [];
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
combination_list = {{'OFC', 'DLPFC'}};
region_list = {'OFC','ACC'};
num_regions = length(combination_list);

num_sessions = length(sessions);
% skip_sessions_F = [13,23,24,30,31,33,34];
% skip_sessions_M = [10,12,22,36,41,44,46,48];
include_sessions = [07,11,12,14,15,21,23,27,30,41,42];
% include_sessions = sessions;
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
%     if 1 
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
                                f_res_bucket_width = fs/size(region_comp,2);
                                fres = fnq/f_res_bucket_width;
                                %% Start the calcs of granger cause
                                % Find the cpsd
                                ptic('\n*** tsdata_to_cpsd... ');
                                % Added stupid windows to try to smooth
%                                 cpsd_trial = tsdata_to_cpsd(region_comp, floor(size(region_comp, 2)/2),'MT', size(region_comp, 2), floor(size(region_comp, 2)/2),3, 4);
                                cpsd_trial = tsdata_to_cpsd(region_comp, fres,'MT', size(region_comp, 2), floor(size(region_comp, 2)/2),3, 4);
                                %% Failed attempt to use the fieldtrip information instead of mvgc doesnt match output/wrong shapes
%                                 cfg_field = [];   
%                                 cfg_field.method = 'mtmfft';
%                                 cfg_field.output = 'powandcsd';
%                                 cfg_field.channel = {chan1{1}; chan2{1}};
% %                                 cfg_field.keeptapers = 'yes';
% %                                 cfg_field.keeptrials = 'yes';%'no';
%                                 cfg_field.pad='maxperlen'; % Padding: not adding zeros
%                                 cfg_field.flag = 0;
% 
%                                 % have hanning and 1 taper as Conrado suggested
%                                 cfg_field.taper = 'dpss';
%                                 cfg_field.tapsmofrq = 7;
%                                 cfg_field.foilim = [1 140];
                                
                                ptoc;
                                % find the autocov
                                ptic('\n*** cpsd_to_autocov... ');
                                [G,q] = cpsd_to_autocov(cpsd_trial);
                                ptoc;
                            catch 
                                fprintf('\nCurrent channel combination is of regions %s %s: %d, %d, was not pos def for chol!\n', region_names(r1), region_names(r2), chan1, chan2)
                                bad_combinations(end+1,:) = [r1 r2 chan1 chan2];
                            end

                            try

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
%                 x_range = (1:1:length(gc_1_ave))./length(gc_1_ave);
%                 x_range = x_range.*fnq;
                x_range = sfreqs(fres,fs);
                % Feed Forward
                hold on
                ff_pn = ((time-1)*6);
                
                %Feedforward theta
                hold on
                subplot(length(data_sets),6, ff_pn+1)
                plot(x_range,gc_2_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('Theta FF during %s', plot_titles{time}))
                xlim([0 10]);
%                 ylabel('Granger Causality')
%                 xlabel('Frequency (Hz)')
                
                % Feedback theta
                hold on
                subplot(length(data_sets), 6, ff_pn+2)
                %loglog
                plot(x_range,gc_1_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('Theta FB during %s', plot_titles{time}))
                xlim([0 10]);
%                 ylabel('Granger Causality')
%                 xlabel('Frequency (Hz)')
                
                %Feedforward Beta
                hold on
                subplot(length(data_sets), 6, ff_pn+3)
                plot(x_range,gc_2_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('Beta FF during %s', plot_titles{time}))
                xlim([10 30]);
%                 ylabel('Granger Causality')
%                 xlabel('Frequency (Hz)')
                
                % Feedback Beta
                hold on
                subplot(length(data_sets),6, ff_pn+4)
                %loglog
                plot(x_range,gc_1_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('Beta FB during %s', plot_titles{time}))
                xlim([10 30]);
%                 ylabel('Granger Causality')
%                 xlabel('Frequency (Hz)')
                
                %Feedforward gamma
                hold on
                subplot(length(data_sets),6, ff_pn+5)
                plot(x_range,gc_2_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('Gamma FF during %s', plot_titles{time}))
                xlim([30 100]);
%                 ylabel('Granger Causality')
%                 xlabel('Frequency (Hz)')
                
                % Feedback gamma
                hold on
                subplot(length(data_sets),6, ff_pn+6)
                %loglog
                plot(x_range,gc_1_ave, 'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                title(sprintf('Gamma FF during %s', plot_titles{time}))
                xlim([30 100]);
%                 ylabel('Granger Causality')
%                 xlabel('Frequency (Hz)')

            end 
        end
        color_num = color_num + 1;
    end
end
subplot(length(data_sets),6, 31)
legend('NumColumns',4)
subplot(length(data_sets),6, 33)
legend('NumColumns',4)
subplot(length(data_sets),6, 35)
legend('NumColumns',4)
sgtitle(sprintf('%s Granger Causality Analysis %s to %s', monkey, interest_region_1,interest_region_2))
han = axes(figure(1), 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han,'Granger Causality');xlabel(han, 'Frequency (Hz)');
% title(han,sprintf('%s Granger Causeality of Regions Against all Other \n', monkey));
