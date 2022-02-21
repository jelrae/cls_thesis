%% Initialization

clear all
clc

addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
set_paths_roi;

ft_defaults

%% Load in basic things
addpath('/home/jordan/neuro_thesis/data/hunt_et_al/scripts/merged_F')

load('monkey_sessions_F.mat')
monkey = 'F';
save_dir = 'power_check';

% Fourier transfrom vars

pow_cfg = [];
pow_cfg.average = 'yes';
pow_cfg.method = 'mtmfft';
pow_cfg.keeptapers = 'yes';
pow_cfg.output = 'fourier';
pow_cfg.channel = 'all';
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
pre_cfg.dftfreq = [50,100,150];

info_store = [];

for sn = 1:length(sessions_F)
    session_number = sessions_F{sn};
    disp(session_number)
    file_name = sprintf('%s_%s_data.mat',monkey, session_number);
    base_name = sprintf('%s_%s_data',monkey, session_number);
    data_name = append(base_name,'_1');
    load(file_name);
    
    % select the channels
    eval(sprintf('pow_cfg.channel = %s.label;',data_name));
    desc_cfg.channel = pow_cfg.channel;
    
    % do some preproc
    
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
    
    info_store(end+1,:) = avgcoh.freq;
    
    hold on
    A_in = squeeze(avgcoh.cohspctrm);
    plot(avgcoh.freq, A_in, 'DisplayName', base_name(1:6));
%     clearvars -except monkey save_dir session_M session_F pow_cfg desc_cfg
end
hold on
% A_in = squeeze(avgcoh.cohspctrm);
% plot(avgcoh.freq, A_in, 'DisplayName', 'Average coherence');

ylabel('Power')
xlabel('Frequency (Hz)')
title('Coherence analysis F')
legend