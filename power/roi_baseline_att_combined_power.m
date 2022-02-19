%% Initalize

addpath('C:/Users/Jordan/Documents/cls_thesis/neuro_thesis/cls_thesis/helper_functions/')
set_paths_roi;

clear all
clc
ft_defaults

%Load previously made data
load_data = true; 

% Declare your Monkey!
monkey = 'kurt';
monkey_caps = 'Kurt';
% monkey = 'pele';
% monkey_caps = 'Pele';

plot_loglog = 'no'; % Set to yes to plot loglog

lineColor = linspecer(3,'qualitative');

% Run script to get ROIs specific to the monkey
fig_6_ROIS;

monkeys = {'kurt', 'pele'};
for m = monkeys
    
    for r = region_names

    end
end 
% Define region of interest
roi = V4;
roi_string = "V4"; % Set for plotting

%% Load data

load(sprintf('%s_p_all_AttIn.mat',monkey));
%commented out for the noise check test
% load(sprintf('%s_p_all_AttOut.mat',monkey));

%% for the checking of if noise changes the power 
% all_AttOut = all_AttIn;
% 
% for i = 1:length(all_AttOut)
%     all_AttOut.trial{i} = all_AttOut.trial{i} + (.25 .* randn(size(all_AttOut.trial{i})));
% end


%% Power Analysis Attention In
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

% for the beta and gamma range
pow_cfg.tapsmofrq = 4;
pow_cfg.foilim = [8 140];
pow_cfg.taper = 'dpss';

b_g_power_analysis = ft_freqanalysis(pow_cfg,all_AttIn);

% for the theta  range
pow_cfg.tapsmofrq = 1;
pow_cfg.foilim = [2 8];
pow_cfg.taper = 'hanning';

t_power_analysis = ft_freqanalysis(pow_cfg,all_AttIn);

clear all_AttIn all_AttOut

%% Extract power and plot

pow2_cfg = [];
pow2_cfg.psi = 'no'; % Phase slope Index
pow2_cfg.channel = roi;
pow2_cfg.jackknife = 'yes';
pow2_cfg.avgChOI = 'yes';

% Obtaining and plotting the power analysis

% Power
b_g_power_analysis = ft_freqdescriptives(pow2_cfg, b_g_power_analysis);
t_power_analysis = ft_freqdescriptives(pow2_cfg, t_power_analysis);

fig1 = figure(1);

hold on

% Plot graph AttIn Theta
AttIn_label = 'AttIn %s';


subplot(1,2,1)
if strcmp(plot_loglog,'yes')
    loglog(t_power_analysis.freq,(mean(t_power_analysis.powspctrm(:,:))), ...
        'DisplayName', 'AttIn', 'Color', lineColor(1,:), 'LineWidth', 1.25)
elseif strcmp(plot_loglog,'no')
    plot(t_power_analysis.freq, (mean(t_power_analysis.powspctrm(:,:))),...
    'DisplayName', 'AttIn', 'Color', lineColor(1,:), 'LineWidth', 1.25)
end
hold on

% Plot graph AttIn Beta Gamma

subplot(1,2,2)
if strcmp(plot_loglog,'yes')
    loglog(b_g_power_analysis.freq, (mean(b_g_power_analysis.powspctrm(:,:))),...
        'DisplayName', 'AttIn', 'Color', lineColor(1,:), 'LineWidth', 1.25)
elseif strcmp(plot_loglog,'no')
    plot(b_g_power_analysis.freq, (mean(b_g_power_analysis.powspctrm(:,:))),...
    'DisplayName', 'AttIn', 'Color', lineColor(1,:), 'LineWidth', 1.25)
end

%% Adding in Attention Out

load(sprintf('%s_p_all_AttOut.mat',monkey));

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

% for the beta and gamma range
pow_cfg.tapsmofrq = 4;
pow_cfg.foilim = [8 140];
pow_cfg.taper = 'dpss';

b_g_power_analysis = ft_freqanalysis(pow_cfg,all_AttOut);

% for the theta  range
pow_cfg.tapsmofrq = 1;
pow_cfg.foilim = [2 8];
pow_cfg.taper = 'hanning';

t_power_analysis = ft_freqanalysis(pow_cfg,all_AttOut);

clear all_AttOut

%% Extract power and plot

pow2_cfg = [];
pow2_cfg.psi = 'no'; % Phase slope Index
pow2_cfg.channel = roi;
pow2_cfg.jackknife = 'yes';
pow2_cfg.avgChOI = 'yes';

% Obtaining and plotting the power analysis

% Power
b_g_power_analysis = ft_freqdescriptives(pow2_cfg, b_g_power_analysis);
t_power_analysis = ft_freqdescriptives(pow2_cfg, t_power_analysis);

hold on

% Plot graph AttOut Theta

subplot(1,2,1)
if strcmp(plot_loglog,'yes')
    loglog(t_power_analysis.freq,(mean(t_power_analysis.powspctrm(:,:))), ...
        'DisplayName', 'AttOut', 'Color', lineColor(2,:), 'LineWidth', 1.25)
elseif strcmp(plot_loglog,'no')
    plot(t_power_analysis.freq, (mean(t_power_analysis.powspctrm(:,:))),...
    'DisplayName', 'AttOut', 'Color', lineColor(2,:), 'LineWidth', 1.25)
end
hold on

% Plot graph AttOut Beta Gamma

subplot(1,2,2)
if strcmp(plot_loglog,'yes')
    loglog(b_g_power_analysis.freq, (mean(b_g_power_analysis.powspctrm(:,:))),...
        'DisplayName', 'AttOut', 'Color', lineColor(2,:), 'LineWidth', 1.25)
elseif strcmp(plot_loglog,'no')
    plot(b_g_power_analysis.freq, (mean(b_g_power_analysis.powspctrm(:,:))),...
    'DisplayName', 'AttOut', 'Color', lineColor(2,:), 'LineWidth', 1.25)
end


%% adding in baselines

load(sprintf('%s_b_all_AttIn.mat',monkey));

cnd = {'Baseline'};

% Fourier transfrom
cnd = {'AttIn' 'AttOut'}; % conditions
pow_cfg = [];
pow_cfg.average = 'yes';
pow_cfg.method = 'mtmfft';
pow_cfg.keeptapers = 'yes';
pow_cfg.output = 'fourier';
pow_cfg.channel = roi;

pow_cfg.keeptrials = 'yes';%'no';
pow_cfg.pad='maxperlen'; % Padding: not adding zeros
pow_cfg.flag = 0;

% for the beta and gamma range
pow_cfg.tapsmofrq = 4;
pow_cfg.foilim = [8 140];
pow_cfg.taper = 'dpss';

b_g_power_analysis = ft_freqanalysis(pow_cfg,all_AttIn);

% for the theta  range
pow_cfg.tapsmofrq = 1;
pow_cfg.foilim = [2 8];
pow_cfg.taper = 'hanning';

t_power_analysis = ft_freqanalysis(pow_cfg,all_AttIn);

clear all_AttIn

% Freqdescriptives
cnd = {'Baseline'};

pow2_cfg = [];
pow2_cfg.psi = 'no'; % Phase slope Index
pow2_cfg.channel = roi;
pow2_cfg.jackknife = 'yes';
pow2_cfg.avgChOI = 'yes';

% Obtaining and plotting the power analysis

% Power analysis
b_g_power_analysis = ft_freqdescriptives(pow2_cfg, b_g_power_analysis);
t_power_analysis = ft_freqdescriptives(pow2_cfg, t_power_analysis);

hold on

% Plot Baseline theta range

subplot(1,2,1)
if strcmp(plot_loglog,'yes')
    loglog(t_power_analysis.freq,(mean(t_power_analysis.powspctrm(:,:))), ...
        'DisplayName', 'Baseline', 'Color', lineColor(3,:), 'LineWidth', 1.25)
elseif strcmp(plot_loglog,'no')
    plot(t_power_analysis.freq, (mean(t_power_analysis.powspctrm(:,:))),...
    'DisplayName', 'Baseline', 'Color', lineColor(3,:), 'LineWidth', 1.25)
end
hold on

% Plot Baseline beta gamma range
subplot(1,2,2)
if strcmp(plot_loglog,'yes')
    loglog(b_g_power_analysis.freq, (mean(b_g_power_analysis.powspctrm(:,:))),...
        'DisplayName', 'Baseline', 'Color', lineColor(3,:), 'LineWidth', 1.25)
elseif strcmp(plot_loglog,'no')
    plot(b_g_power_analysis.freq, (mean(b_g_power_analysis.powspctrm(:,:))),...
    'DisplayName', 'Baseline', 'Color', lineColor(3,:), 'LineWidth', 1.25)
end

legend

sfh1 = subplot(1,2,1, 'Parent', gcf);
sfh2 = subplot(1,2,2, 'Parent', gcf);

init_width = 0.3347;

% sfh1.Position = sfh1.Position - [0 0 init_width-0.05 0];
sfh1.Position = [0.11 0.11 0.075 0.815];
sfh2.Position = [0.225 0.11 0.725 0.815];
% sfh2.Position = sfh2.Position + [0 0 1.175-init_width 0];
% sfh2.Position = sfh2.Position - [init_width-0.05 0 0 0];

% Plotting labels
han = axes(figure(1), 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylh = ylabel(han,'Power');xlabel(han, 'Frequency (Hz)');
ylh.Position(1) = -.075;
formatSpec = 'Power Analysis %s %s';
title(han, sprintf(formatSpec, monkey_caps, roi_string))

saveas(fig1, sprintf('results/tbg_%s_%s.jpg', monkey_caps, roi_string))
close(fig1)