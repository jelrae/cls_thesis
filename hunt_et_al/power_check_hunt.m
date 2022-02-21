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
pow_cfg.foilim = [20 140];

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

% interest_region = 'DLPFC';
% interest_region = 'ACC';
interest_region = 'OFC';
skip_sessions = [13,23,24,30,31,33,34];
colors = distinguishable_colors(16);
color = 1;

for sn = 1:length(sessions)
    session_number = sessions{sn};
    disp(session_number)
    if ~ismember(str2num(session_number(2:end)),skip_sessions)
        file_name = sprintf('%s_%s_data.mat',monkey, session_number);
        base_name = sprintf('%s_%s_data',monkey, session_number);
    %     data_name = append(base_name,'_1');
    %     data_name = append(base_name,'_2');
        data_name = append(base_name,'_pre_choice_300');
    %     data_name = append(base_name,'_pre_choice_500');
    %     data_name = append(base_name,'_pre_feedback');
        load(file_name);

        % select the channels
    %     eval(sprintf('pow_cfg.channel = %s.label;',data_name));
    %     desc_cfg.channel = pow_cfg.channel;
        eval(sprintf('pow_cfg.channel = %s_%s_channels;',base_name,interest_region));
        desc_cfg.channel = pow_cfg.channel;

        % the things if the region has things in it
        if size(desc_cfg.channel,1) ~= 0 
            % do some preproc
            eval(sprintf('%s = ft_preprocessing(pre_cfg, %s);', data_name, data_name));

            % do the Fouier

            eval(sprintf('freqanal = ft_freqanalysis(pow_cfg,%s);', data_name));

            %% Freqdescriptives
            freqanal = ft_freqdescriptives(desc_cfg, freqanal);

            hold on
            if size(freqanal.powspctrm,1) > 1
                loglog(freqanal.freq,(mean(freqanal.powspctrm(:,:))), ...
                    'DisplayName', base_name(1:6), 'Color', colors(color,:), 'LineWidth',2)
                color = color+1;
    %             set(gca, 'LineWidth',3);
            else
                loglog(freqanal.freq,freqanal.powspctrm(:,:), ...
                    'DisplayName', base_name(1:6), 'Color', colors(color,:), 'LineWidth',2)
                color = color+1;
    %             set(gca, 'LineWidth',3);
            end
            legend
        %     clearvars -except monkey save_dir session_M session_F pow_cfg desc_cfg
        end
    end
end

ylabel('Power')
xlabel('Frequency (Hz)')
title(sprintf('Power analysis %s %s end', monkey, interest_region))
legend