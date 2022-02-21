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

% Fourier transfrom vars

pow_cfg = [];
pow_cfg.average = 'yes';
pow_cfg.method = 'mtmfft';
pow_cfg.keeptapers = 'yes';
pow_cfg.output = 'fourier';
pow_cfg.keeptrials = 'yes';%'no';
pow_cfg.pad='maxperlen'; % Padding: not adding zeros
pow_cfg.flag = 0;

% have hanning and 1 taper as Conrado suggested
pow_cfg.taper = 'hanning';
pow_cfg.tapsmofrq = 1;
pow_cfg.foilim = [5 80];

desc_cfg = [];
desc_cfg.psi = 'no'; % Phase slope Index
desc_cfg.jackknife = 'yes';
desc_cfg.avgChOI = 'yes';

pre_cfg = [];
pre_cfg.demean = 'yes';
pre_cfg.detrend = 'yes';
pre_cfg.dftfilter = 'yes';

% Region checking points

region_list = {'DLPFC', 'ACC', 'OFC'};
combination_list = {{'ACC','DLPFC'}, {'ACC', 'OFC'}, {'OFC', 'DLPFC'}};
% region_list = {'OFC','ACC'};
num_regions = length(combination_list);

num_sessions = length(sessions);
% skip_sessions_F = [13,23,24,30,31,33,34];
% skip_sessions_M = [10,12,22,36,41,44,46,48];
include_sessions = [07,11,12,14,15,21,23,27,30,41,42];
% 40, 37
% skip_sessions = skip_sessions_M;
skip_sessions = [];

% plotting_colors = distinguishable_colors(num_sessions - length(skip_sessions));
plotting_colors = distinguishable_colors(num_sessions);
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
%     if ~ismember(s_n,skip_sessions)
    if ismember(s_n,include_sessions)
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
                interest_regions = combination_list{c_n};
                interest_region_1 = interest_regions{1};
                interest_region_2 = interest_regions{2};
                % select the channels
                eval(sprintf('channels_1 = %s_%s_channels;',base_name, interest_region_1));
                eval(sprintf('channels_2 = %s_%s_channels;',base_name, interest_region_2));
                if size(channels_1,1) > 0 && size(channels_2,1) > 0
                    pow_cfg.channel = [channels_1,channels_2];
                    desc_cfg.channel = pow_cfg.channel;

                    % do some preproc
                    disp(session_number)
                    eval(sprintf('%s = ft_preprocessing(pre_cfg, %s);', data_name, data_name));

                    % do the Fouier

                    eval(sprintf('freqanal = ft_freqanalysis(pow_cfg,%s);', data_name));
                    freqanal = ft_selectdata(pow_cfg, freqanal);
                    
                    %% connectivity analysis
                    cfg = [];
                    cfg.method = 'coh';
                    coh = ft_connectivityanalysis(cfg, freqanal);

                    %% Average over channels
                    cfg = [];
                    cfg.avgoverchan = 'yes';
                    avgcoh = ft_selectdata(cfg, coh);
                    ff_pn = (2*num_regions*(time-1)+(2*c_n)-1);

                    hold on
                    subplot(length(data_sets),length(combination_list)*2, ff_pn)
%                         subplot(length(data_sets),length(region_list), (time))
                    loglog(avgcoh.freq, squeeze(avgcoh.cohspctrm), ...
                        'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
%                     subplot(length(region_list),length(data_sets), (3*(time-1)+r_n), freqanal.freq,freqanal.powspctrm(:,:), 'DisplayName', base_name(1:6))
                    title(sprintf('%s - %s during %s', interest_region_1, interest_region_2, plot_titles{time}))
                    set(gca,'xscale','log')
                    set(gca,'yscale','log')
                    ylabel('Coherence')
                    xlabel('Frequency (Hz)')
%                     set(gca, 'XScale', 'linear')
                    xlim([5 20])
                    
                    hold on
                    subplot(length(data_sets),length(combination_list)*2, ff_pn+1)
%                         subplot(length(data_sets),length(region_list), (time))
                    loglog(avgcoh.freq, squeeze(avgcoh.cohspctrm), ...
                        'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
%                     subplot(length(region_list),length(data_sets), (3*(time-1)+r_n), freqanal.freq,freqanal.powspctrm(:,:), 'DisplayName', base_name(1:6))
                    title(sprintf('%s - %s during %s', interest_region_1, interest_region_2, plot_titles{time}))
                    set(gca,'xscale','log')
                    set(gca,'yscale','log')
                    ylabel('Coherence')
                    xlabel('Frequency (Hz)')
%                     set(gca, 'XScale', 'linear')
                    xlim([20 80])

                end
            end 
        end
        color_num = color_num + 1;
    end
end
subplot(length(data_sets),length(combination_list)*2, 31)
legend('NumColumns',4)
subplot(length(data_sets),length(combination_list)*2, 33)
legend('NumColumns',4)
subplot(length(data_sets),length(combination_list)*2, 35)
legend('NumColumns',4)
sgtitle(sprintf('Coherence Analysis %s', monkey))
