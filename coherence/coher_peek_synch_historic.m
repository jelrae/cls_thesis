addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\cls_thesis\helper_functions\')
set_paths_roi;

clear all
clc
ft_defaults

%% Select the monkey and define region interest (hyperparamters)
monkey = 'kurt';
% monkey = 'pele';

load_data = true;

% Define and set analysis objectives and specs:

des_data = 'post';
exp_desc = 'Coherence';
roi_string = 'V1 and V4'; % Set for plotting

% Run script to get ROIs specific to the monkey
vbm_rois;

%% load data
load(sprintf('%s_p_all_AttIn.mat',monkey));
load(sprintf('%s_p_all_AttOut.mat',monkey));

%Check for region being in data

% region_1 = V1;
% region_2 = V4;
% 
% for i = region_1
%     if ~ismember(i, all_AttIn.label)
%     region_1(strcmp(region_1,i)) = [];
%     end
% end
% 
% for i = region_2
%     if ~ismember(i, all_AttIn.label)
%     region_2(strcmp(region_2,i)) = [];
%     end
% end

% Define region of interest 1
roi = V1;
%Define region of interest 2
roi_2 = DP;
roi_combine = [roi roi_2];

roi_list = cell(0,2);

for i = V1
    if ismember(i, all_AttIn.label)
        for j = DP
            if ismember(j, all_AttIn.label)
                roi_list(end+1,:) = {i{1} j{1}};
            end
        end
    end
end

%% Calculate Fourier on ALL channels
cnd = {'AttIn' 'AttOut'}; % conditions
% cnd = {'AttIn'}; % conditions
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


coh_cfg.tapsmofrq = 8;
coh_cfg.frequency = [2 140];

% coh_cfg.taper = 'hanning';
coh_cfg.taper = 'dpss';

for n = 1:length(cnd)
    tmp_angle = [];
    tmp = [];
    tmp_pow = [];
    sprintf('Condition %s',cnd{n})
    % Here is where attmin and out are used.
    eval(sprintf('tmp_dat = all_%s;',cnd{n}))

    tmp_freq = ft_freqanalysis(coh_cfg,tmp_dat);
    tmp = ft_selectdata(coh_cfg, tmp_freq);
    eval(sprintf('freq_%s=tmp;',cnd{n}));
end


clear all_AttIn all_AttOut

%% Connectivity (coherence) analysis Tap8
cfg = [];
cfg.method = 'coh';
cfg.channelcmb = roi_list;
coh_in = ft_connectivityanalysis(cfg, freq_AttIn);
coh_out = ft_connectivityanalysis(cfg, freq_AttOut);

%% Calculate the peak and align

% gp_info = [];
% gp_info.location = [];
% gp_info.value = [];
% gp_info.width = [];
% gp_info.no_peak = [];
% 
% gp_start = 30;
% gp_end = 60;
% 
% for i =1:size(coh_in.cohspctrm,1)
% 	[pk, pl, pw, pp] = findpeaks(coh_in.cohspctrm(i,gp_start:gp_end));
%     if ~isempty(pk)
%         [vm, lm] = max(pk);
%         gp_info.value(end+1) = vm;
%         gp_info.location(end+1) = gp_start - 1 + pl(lm);
%         gp_info.width(end+1) = pw(lm);
%     else
%         gp_info.value(end+1) = 0;
%         gp_info.location(end+1) = 0;
%         gp_info.width(end+1) = 0;
%         gp_info.no_peak(end+1) = i;
%     end
% end
% 
% gp_info.ave_loc = round(sum(gp_info.location)/nnz(gp_info.location));
% gp_info.location(gp_info.location==0) = gp_info.ave_loc;
% 
% w = 10;
% 
% g_peak = zeros(1,2*w+1);
% 
% for i = 1:size(gp_info.location,2)
%     start_loc = gp_info.location(i) - w;
%     end_loc = gp_info.location(i) + w;
%     g_peak = g_peak + coh_in.cohspctrm(i,start_loc:end_loc);
% end
% g_peak = g_peak./size(coh_in.cohspctrm,1);
% 
% g_min = gp_info.ave_loc - w;
% g_max = gp_info.ave_loc + w;
% plot(coh_in.freq(g_min:g_max), g_peak);

peak_start = 30;
peak_end = 60;
w = 10;

p_info = PeakAlignSpectrum(peak_start,peak_end,w,coh_in.cohspctrm);

g_min = p_info.ave_loc - w;
g_max = p_info.ave_loc + w;
plot(coh_in.freq(g_min:g_max), p_info.peak);

%% Select data - average over channels
cfg = [];
% cfg.avgoverchan = 'yes';
cfg.avgoverchancmb = 'yes';

avgcoh_in = ft_selectdata(cfg, coh_in);
avgcoh_out = ft_selectdata(cfg, coh_out);



%% Plot coherence for ALL channels
% --TAP8--
% AttIn tap8
AttIn_label = 'Att. Inside';
A_in = squeeze(avgcoh_in.cohspctrm);
% axes('Position', [ 0.6218    0.1100    0.3629    0.8150]);

plot(avgcoh_in.freq, A_in, 'DisplayName', AttIn_label);
% if strcmp(monkey,'kurt')
%     xlim([30 120]);
%     ylim([0.06 0.22]);
% elseif strcmp(monkey,'pele')
%     xlim([40 100]);
%     ylim([0 0.14]);

hold on
% AttOut 
AttOut_label = 'Att. Outside';
A_out = squeeze(avgcoh_out.cohspctrm);
% axes('Position', [ 0.6218    0.1100    0.3629    0.8150]);
plot(avgcoh_out.freq, A_out, 'DisplayName', AttOut_label);
% Plotting labels
ylabel('Coherence')
xlabel('Frequency (Hz)')
formatSpec = 'Coherence analysis V1 - DP for %s';
title('Coherence Analysis Kurt V1 - DP')
 
legend