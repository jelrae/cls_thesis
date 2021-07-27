%% Baseline power

%% Add paths
addpath('/home/jordan/common/matlab/fieldtrip-20210411')
addpath('/home/jordan/neuro_thesis/gc_hierarchies/helper_functions/') 
% addpath('/home/jordan/neuro_thesis/data/monkey_kurt_data/bipolar_baseline_data/') 
addpath('/home/jordan/neuro_thesis/data/monkey_kurt_data/raw_k_baseline_data/') 
addpath('/home/jordan/neuro_thesis/data/monkey_kurt_data/kurt_montage_layout/')
addpath('/home/jordan/neuro_thesis/data/monkey_pele_data/raw_p_baseline_data/') 
addpath('/home/jordan/neuro_thesis/data/monkey_pele_data/pele_montage_layout/') 

%% Get sessions (Normalization is in merge_data)
clear all
clc
ft_defaults
cfg = [];
cfg.bsln = 1;
cfg.norm = 'yes';
cfg.montage = 'yes';
cfg.equal = 'yes';
%% Kurt
% cfg.monkey = 'ku';
% cfg.path = '/home/jordan/neuro_thesis/data/monkey_kurt_data/raw_k_baseline_data/'; %% Does not exist
% cfg.bipolarder = 'home/jordan/neuro_thesis/data/monkey_kurt_data/kurt_montage_layout/kurt_montage_bipolar_alongfinger_v2.mat';  
% cfg.v_session = [67 68 69];
% plot_loglog = 'yes'; % Set to yes to plot loglog
%% Pele
cfg.monkey = 'pe';
cfg.path = '/home/jordan/neuro_thesis/data/monkey_pele_data/raw_p_baseline_data/';
cfg.bipolarder = 'home/jordan/neuro_thesis/data/monkey_pele_data/pele_montage_layout/pele_montage_bipolar_alongfinger_v2.mat';
cfg.v_session = [39 40 42 87 88 90];
plot_loglog = 'no'; 
%% Frequency descriptives
mcfg = cfg;
[baseline] = merge_data2(mcfg);

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

% Plot graph AttIn
plot(powtap8_baseline.freq, mean(powtap8_baseline.powspctrm(:,:)))
% hold on
% 
% Plot graph AttOut
baseline_label = 'Baseline';
plot(powtap8_baseline.freq, mean(powtap8_baseline.powspctrm(:,:)), 'DisplayName', sprintf(baseline_label));

% loglog
if strcmp(plot_loglog,'yes')
    loglog(powtap8_baseline.freq, mean(powtap8_baseline.powspctrm(:,:)), 'DisplayName', sprintf(baseline_label));
elseif strcmp(plot_loglog,'no')
    plot(powtap8_baseline.freq, mean(powtap8_baseline.powspctrm(:,:)), 'DisplayName', sprintf(baseline_label));

end
% labels
ylabel('Power')
xlabel('Frequency (Hz)')
legend

% clear all;
% 
