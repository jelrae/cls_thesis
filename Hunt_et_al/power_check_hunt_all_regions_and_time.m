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
% pre_cfg.dftfreq = [50];

% Region checking points

region_list = {'DLPFC', 'ACC', 'OFC'};
% region_list = {'ACC','OFC'};
num_regions = length(region_list);

num_sessions = length(sessions);
skip_sessions_F = [13,23,24,30,31,33,34];
% skip_sessions_M = [14,15,19,21,29,30,40,42,44,45];
include_sessions = [07,11,12,14,15,21,23,27,30,41,42];
skip_sessions_M = [];
maybe_skip_M = [21,25,16];
% 40, 37
skip_sessions = skip_sessions_M;
plotting_colors = distinguishable_colors(num_sessions - length(skip_sessions));
% plotting_colors = distinguishable_colors(num_sessions);
% plotting_colors = distinguishable_colors(length(skip_sessions));
color_num = 1;

for sn = 1:length(sessions)/2
% for sn = length(sessions)/2:length(sessions)
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
    %     data_sets(end+1) = {append(base_name,'_pre_choice_500')};
        data_sets(end+1) = {append(base_name,'_pre_feedback')};

        plot_titles = {};
        plot_titles(end+1) = {'Pre Fixation'};
        plot_titles(end+1) = {'During Fixation'};
        plot_titles(end+1) = {'First Image'};
        plot_titles(end+1) = {'Second Image'};
        plot_titles(end+1) = {'Pre Choice'};
    %     plot_titles(end+1) = {'Pre Choice 500'};
        plot_titles(end+1) = {'Pre Feedback'};

        for time = 1:length(data_sets)
            data_name = data_sets{time};
            for r_n = 1:length(region_list)
                interest_region = region_list{r_n};
                % select the channels
                eval(sprintf('pow_cfg.channel = %s_%s_channels;',base_name,interest_region));
                desc_cfg.channel = pow_cfg.channel;

                % the things if the region has things in it
                if size(desc_cfg.channel,1) ~= 0 
                    % do some preproc
                    disp(session_number)
                    eval(sprintf('%s = ft_preprocessing(pre_cfg, %s);', data_name, data_name));

                    % do the Fouier

                    eval(sprintf('freqanal = ft_freqanalysis(pow_cfg,%s);', data_name));

                    %% Freqdescriptives
                    freqanal = ft_freqdescriptives(desc_cfg, freqanal);
                    ff_pn = (2*num_regions*(time-1)+(2*r_n)-1);
                    if size(freqanal.powspctrm,1) > 1
                        hold on
                        subplot(length(data_sets),length(region_list)*2, ff_pn)
%                         subplot(length(data_sets),length(region_list), (time))
                        loglog(freqanal.freq,(mean(freqanal.powspctrm(:,:))), ...
                            'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                        title(sprintf('%s during %s', interest_region, plot_titles{time}))
                        ylabel('Power')
                        xlabel('Frequency (Hz)')
%                         legend('NumColumns',2)
%                         set(gca, 'XScale', 'linear')
                        xlim([5 20])
                        
                        hold on
                        subplot(length(data_sets),length(region_list)*2, ff_pn+1)
%                         subplot(length(data_sets),length(region_list), (time))
                        loglog(freqanal.freq,(mean(freqanal.powspctrm(:,:))), ...
                            'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                        title(sprintf('%s during %s', interest_region, plot_titles{time}))
                        ylabel('Power')
                        xlabel('Frequency (Hz)')
%                         legend('NumColumns',2)
%                         set(gca, 'XScale', 'linear')
                        xlim([20 80])
                    else
                        hold on
                        subplot(length(data_sets),length(region_list)*2, ff_pn)
%                         subplot(length(data_sets),length(region_list), (time))
                        loglog(freqanal.freq,freqanal.powspctrm(:,:), ...
                            'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                        title(sprintf('%s during %s', interest_region, plot_titles{time}))
    %                     set(gca,'xscale','log')
    %                     set(gca,'yscale','log')
                        ylabel('Power')
                        xlabel('Frequency (Hz)')
    %                     label(base_name(1:6))
%                         legend('NumColumns',2)
%                         set(gca, 'XScale', 'linear')
                        
                        hold on
                        subplot(length(data_sets),length(region_list)*2, ff_pn+1)
%                         subplot(length(data_sets),length(region_list), (time))
                        loglog(freqanal.freq,freqanal.powspctrm(:,:), ...
                            'DisplayName', append('S', base_name(5:6)), 'Color', plotting_colors(color_num,:), 'LineWidth',1.5)
                        title(sprintf('%s during %s', interest_region, plot_titles{time}))
                        ylabel('Power')
                        xlabel('Frequency (Hz)')
%                         legend('NumColumns',2)
%                         set(gca, 'XScale', 'linear')
                        xlim([20 80])
                    end
                end
            end 
        end
        color_num = color_num + 1;
    end
end
% legend
sgtitle(sprintf('Power Analysis %s', monkey))
subplot(length(data_sets),length(region_list)*2, 31)
legend('NumColumns',4)
subplot(length(data_sets),length(region_list)*2, 33)
legend('NumColumns',4)
subplot(length(data_sets),length(region_list)*2, 35)
legend('NumColumns',4)