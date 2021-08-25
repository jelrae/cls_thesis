addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\gc_hierarchies\helper_functions\')
set_paths_roi;

clear all
clc
ft_defaults

%% Select the monkey and define region interest (hyperparamters)
% monkey = 'kurt';
monkey = 'kurt';

load_data = true;

% Define and set analysis objectives and specs:

des_data = 'post';
exp_desc = 'Coherence';
roi_string = 'V1 and V4'; % Set for plotting

% Run script to get ROIs specific to the monkey
vbm_rois;
% Define region of interest 1
roi = V1;
%Define region of interest 2
roi_2 = V4;
roi_combine = [roi roi_2];

%% Create / load data

if ~load_data
    % Configuration parameters - the following are consistent
    cfg = config_setup(roi, monkey, des_data, 'all', false);
    mcfg = cfg;
    [all_AttIn, all_AttOut] = merge_data2(mcfg);
    
else
   load(sprintf('%s_p_all_AttIn.mat',monkey));
   load(sprintf('%s_p_all_AttOut.mat',monkey));
end

% should I save off this fouier data so I dont have to recreate each time? 
% seems like a good idea

%% Calculate Fourier on ALL channels
cnd = {'AttIn' 'AttOut'}; % conditions
tapsmofrq = [8 4 2]; %last with 12 instead 0f 9
foilim = {[20 140],[8 52],[2 15]}; %this is old{[20 140] [5 30] [2 10]};
coh_cfg = [];
coh_cfg.average = 'yes';
coh_cfg.method = 'mtmfft';
coh_cfg.keeptapers = 'yes';
coh_cfg.output = 'fourier';
coh_cfg.channel = roi_combine;
coh_cfg.keeptrials = 'yes';%'no';
coh_cfg.pad='maxperlen'; % Padding: not adding zeros
coh_cfg.flag = 0;

for z=1:length(tapsmofrq)
    coh_cfg.tapsmofrq = tapsmofrq(z);
    coh_cfg.frequency = foilim{z};
    if z==3
        coh_cfg.taper = 'hanning';
    else
        coh_cfg.taper = 'dpss';
    end

    for n = 1:length(cnd)
        tmp_angle = [];
        tmp = [];
        tmp_pow = [];
        sprintf('Condition %s',cnd{n})
        % Here is where attmin and out are used.
        eval(sprintf('tmp_dat = all_%s;',cnd{n}))
        
        tmp_freq = ft_freqanalysis(coh_cfg,tmp_dat);
        tmp = ft_selectdata(coh_cfg, tmp_freq);
        eval(sprintf('freq_tap%d_%s=tmp;',tapsmofrq(z),cnd{n}));
    end
end

clear all_AttIn all_AttOut

%% Connectivity (coherence) analysis Tap8
cfg = [];
cfg.method = 'coh';
coh_in = ft_connectivityanalysis(cfg, freq_tap8_AttIn);
coh_out = ft_connectivityanalysis(cfg, freq_tap8_AttOut);

%% Select data - average over channels
cfg = [];
cfg.avgoverchan = 'yes';

avgcoh_in = ft_selectdata(cfg, coh_in);
avgcoh_out = ft_selectdata(cfg, coh_out);
%% Plot coherence for ALL channels
% --TAP8--
tap_num = 'tap8';
% AttIn tap8
AttIn_label = 'AttIn %s';
A_in = squeeze(avgcoh_in.cohspctrm);
axes('Position', [ 0.6218    0.1100    0.3629    0.8150]);

plot(avgcoh_in.freq, A_in, 'DisplayName', sprintf(AttIn_label, tap_num));
if strcmp(monkey,'kurt')
    xlim([30 120]);
    ylim([0.06 0.22]);
elseif strcmp(monkey,'pele')
    xlim([40 100]);
    ylim([0 0.14]);
end
hold on
% AttOut 
AttOut_label = 'AttOut %s';
A_out = squeeze(avgcoh_out.cohspctrm);
% axes('Position', [ 0.6218    0.1100    0.3629    0.8150]);
plot(avgcoh_out.freq, A_out, 'DisplayName', sprintf(AttOut_label, tap_num));
% Plotting labels
ylabel('Coherence')
xlabel('Frequency (Hz)')
formatSpec = 'Coherence analysis for Monkey %s: %s (Areas: %s)';
title(sprintf(formatSpec, monkey, exp_desc, roi_string))
 
legend

%% Tap2
cfg = [];
cfg.method = 'coh';
coh_in2 = ft_connectivityanalysis(cfg, freq_tap2_AttIn);
coh_out2 = ft_connectivityanalysis(cfg, freq_tap2_AttOut);

% Select data - average over channels
cfg = [];
cfg.avgoverchan = 'yes';

avgcoh_in_tap2 = ft_selectdata(cfg, coh_in2);
avgcoh_out_tap2 = ft_selectdata(cfg, coh_out2);

% Plot 
tap_num2 = 'tap2';
hold on
% AttIn tap2
AttIn_label = 'AttIn';
A_in = squeeze(avgcoh_in_tap2.cohspctrm);
axes('Position', [0.1300    0.0967    0.1575    0.8283])
% yticklabels({'0','0.2','0.4','0.6'})
plot(avgcoh_in_tap2.freq, A_in, 'DisplayName', sprintf(AttIn_label));
if strcmp(monkey,'kurt')
    xlim([0 16]);
    ylim([0.06 0.22]);
elseif strcmp(monkey,'pele')
    xlim([0 20]);
    ylim([0 0.14]);
end
hold on
% AttOut
AttOut_label = 'AttOut';
A_out = squeeze(avgcoh_out_tap2.cohspctrm);
% axes('Position', [0.1300    0.0967    0.1575    0.8283])
plot(avgcoh_in_tap2.freq, A_out, 'DisplayName', sprintf(AttOut_label));
ylabel('Coherence')
% xlabel('Frequency (Hz)')
formatSpec = 'Freq %s';
title(sprintf(formatSpec, tap_num2))
legend

%% Tap4
cfg = [];
cfg.method = 'coh';
coh_in4 = ft_connectivityanalysis(cfg, freq_tap4_AttIn);
coh_out4 = ft_connectivityanalysis(cfg, freq_tap4_AttOut);

% Select data - average over channels
cfg = [];
cfg.avgoverchan = 'yes';

avgcoh_in_tap4 = ft_selectdata(cfg, coh_in4);
avgcoh_out_tap4 = ft_selectdata(cfg, coh_out4);

% Plot 
tap_num4 = 'tap4';
hold on
% AttIn tap2
AttIn_label = 'AttIn';
A_in = squeeze(avgcoh_in_tap4.cohspctrm);
axes('Position', [0.3189    0.1062    0.2771    0.8188]);

% yticklabels({'0','0.3','0.6','0.9'})
plot(avgcoh_in_tap4.freq, A_in, 'DisplayName', sprintf(AttIn_label));
if strcmp(monkey,'kurt')
    xlim([15 35]);
    ylim([0.06 0.14]);
elseif strcmp(monkey,'pele')
    xlim([5 50]);
    ylim([0 0.14]);
end
hold on
% AttOut
AttOut_label = 'AttOut';
A_out = squeeze(avgcoh_out_tap4.cohspctrm);
% axes('Position', [0.3189    0.1062    0.2771    0.8188]);
plot(avgcoh_in_tap4.freq, A_out, 'DisplayName', sprintf(AttOut_label));

xlabel('Frequency (Hz)')
formatSpec = 'Freq %s';
title(sprintf(formatSpec, tap_num4))
legend