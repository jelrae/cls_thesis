%% Initalize

addpath('/home/jordan/neuro_thesis/gc_hierarchies/helper_functions/')
set_paths_roi;

clear all
clc
ft_defaults

%Load previously made data
load_data = true; 

% Declare your Monkey!
% monkey = 'kurt';
monkey = 'pele';

% Define and set analysis objectives and specifics:
exp_desc = 'V1 compared to baseline (loglog)';
roi_string = 'V1'; % Set for plotting
plot_loglog = 'yes'; % Set to yes to plot loglog

% Run script to get ROIs specific to the monkey
vbm_rois;
% Define region of interest
roi = V1;

%% Create/load data post data

if ~load_data
    % Configuration parameters - the following are consistent 
    cfg = [];
    cfg.equal = 'yes';
    cfg.bsln = 0;
    cfg.norm = 'yes';
    cfg.montage = 'yes'; %here you need to add the name of the montage structure you want to use... alongfinger_v2
    % Do monkey things and run merge_data2  Should put the creation of
    % these cfg into a 
else
   load(sprintf('%s_p_all_AttIn.mat',monkey));
   load(sprintf('%s_p_all_AttOut.mat',monkey));
end

%% Fourier transfrom post
cnd = {'AttIn' 'AttOut'}; % conditions
tapsmofrq = [8 4 2]; %last with 12 instead 0f 9
foilim = {[20 140],[8 52],[2 15]}; %this is old{[20 140] [5 30] [2 10]};
pow_cfg = [];
pow_cfg.average = 'yes';
pow_cfg.method = 'mtmfft';
pow_cfg.keeptapers = 'yes';
pow_cfg.output = 'fourier';
pow_cfg.channel = roi;

pow_cfg.keeptrials = 'yes';%'no';
pow_cfg.pad='maxperlen'; % Padding: not adding zeros
pow_cfg.flag = 0;

for z=1:length(tapsmofrq)
    pow_cfg.tapsmofrq = tapsmofrq(z);
    pow_cfg.foilim = foilim{z};
    if z==3
        pow_cfg.taper = 'hanning';
    else
        pow_cfg.taper = 'dpss';
    end

    for n = 1:length(cnd)
        tmp_angle = [];
        tmp = [];
        tmp_pow = [];
        sprintf('Condition %s',cnd{n})
        eval(sprintf('tmp_dat = all_%s;',cnd{n}))
        
        tmp = ft_freqanalysis(pow_cfg,tmp_dat);
        eval(sprintf('freq_tap%d_%s=tmp;',tapsmofrq(z),cnd{n}));
    end
end

clear all_AttIn all_AttOut

%% Freqdescriptives
cnd = {'AttIn' 'AttOut'};
tapsmofrq = [8 4 2];
centerfrqs = {[40:8:120],[12:6:48],[4:2:12]};
freqwin = [8 6 2];

pow2_cfg = [];
pow2_cfg.psi = 'no'; % Phase slope Index
pow2_cfg.channel = roi;
pow2_cfg.jackknife = 'yes';
pow2_cfg.avgChOI = 'yes';

%% Obtaining and plotting the power analysis

% Frequency of interest limits
freq_tap8_AttIn.freq;
freq_tap8_AttOut.freq;

% Power
powtap8_AttIn = ft_freqdescriptives(pow2_cfg, freq_tap8_AttIn);
powtap8_AttOut = ft_freqdescriptives(pow2_cfg, freq_tap8_AttOut);
hold on
% Plot graph AttIn
AttIn_label = 'AttIn %s';
% 
% loglog
if strcmp(plot_loglog,'yes')
    loglog(powtap8_AttIn.freq,(mean(powtap8_AttIn.powspctrm(:,:))), ...
        'DisplayName', sprintf(AttIn_label, roi_string))
elseif strcmp(plot_loglog,'no')
    plot(powtap8_AttIn.freq, (mean(powtap8_AttIn.powspctrm(:,:))),...
    'DisplayName', sprintf(AttIn_label, roi_string))
end

% Plot graph AttOut
AttOut_label = 'AttOut %s';

% loglog
if strcmp(plot_loglog,'yes')
    loglog(powtap8_AttOut.freq, (mean(powtap8_AttOut.powspctrm(:,:))),...
        'DisplayName', sprintf(AttOut_label, roi_string))
elseif strcmp(plot_loglog,'no')
    plot(powtap8_AttOut.freq, (mean(powtap8_AttOut.powspctrm(:,:))),...
    'DisplayName', sprintf(AttOut_label, roi_string))
end
hold on

%% Baseline load data

if ~load_data
    % Configuration parameters - the following are consistent 
    cfg = [];
    cfg.equal = 'yes';
    cfg.bsln = 0;
    cfg.norm = 'yes';
    cfg.montage = 'yes'; %here you need to add the name of the montage structure you want to use... alongfinger_v2
    % Do monkey things and run merge_data2  Should put the creation of
    % these cfg into a 
else
   load(sprintf('%s_b_all_AttIn.mat',monkey));
end

baseline = all_AttIn;
clear all_AttIn

%% Frequency descriptives

cnd = {'baseline'};
tapsmofrq = [8 4 2]; %last with 12 instead 0f 9
foilim = {[20 140],[8 52],[2 15]}; %this is old{[20 140] [5 30] [2 10]};
pow_cfg = [];
pow_cfg.method = 'mtmfft';
pow_cfg.keeptapers = 'yes';
pow_cfg.output = 'fourier';
pow_cfg.channel = 'all';
pow_cfg.keeptrials = 'yes';%'no';
pow_cfg.pad='maxperlen'; % Padding: not adding zeros
pow_cfg.flag = 0;


for z=1:length(tapsmofrq)
    pow_cfg.tapsmofrq = tapsmofrq(z);
    pow_cfg.foilim = foilim{z};
    if z==3
        pow_cfg.taper = 'hanning';
    else
        pow_cfg.taper = 'dpss';
    end

    for n = 1:length(cnd)
        tmp_angle = [];
        tmp = [];
        tmp_pow = [];
        tmp = ft_freqanalysis(pow_cfg, baseline);
        eval(sprintf('freq_tap%d_%s=tmp;',tapsmofrq(z),cnd{n}));
    end
end

pow2_cfg = [];
pow2_cfg.psi = 'no'; % Phase slope Index
pow2_cfg.channel = 'all';%{'H26-H25' 'J17-J16'}; 
pow2_cfg.jackknife = 'yes';
pow2_cfg.avgChOI = 'yes';

tapsmofrq = [8 4 2];
centerfrqs = {[40:8:120],[12:6:48],[4:2:12]};
freqwin = [8 6 2];

%% Obtaining and plotting the power analysis

% Frequency of interest limits
freq_tap8_baseline.freq;

% Powerclear all
powtap8_baseline = ft_freqdescriptives(pow2_cfg, freq_tap8_baseline);

% Specific bipolarderivation that is known to present well: J17-J16
match_str(powtap8_baseline.label,'J17-J16');

% Plot graph Baseline
baseline_label = 'Baseline';

if strcmp(plot_loglog,'yes')
    loglog(powtap8_baseline.freq, mean(powtap8_baseline.powspctrm(:,:)), 'DisplayName', sprintf(baseline_label));
elseif strcmp(plot_loglog,'no')
    plot(powtap8_baseline.freq, mean(powtap8_baseline.powspctrm(:,:)), 'DisplayName', sprintf(baseline_label));

end

% Plotting labels
ylabel('Power')
xlabel('Frequency (Hz)')
formatSpec = 'Power analysis for Monkey %s: %s';
title(sprintf(formatSpec, monkey, exp_desc))
legend