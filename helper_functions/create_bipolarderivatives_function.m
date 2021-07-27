function create_bipolarderivatives_function(monkey, des_data, data_used)
%% Create the Polpolar derivatives Preprocessing

%% Start Clear and hparams

% clear all
% clc

%     % Select the monkey and define region interest (hyperparamters) 
%     % monkey = 'kurt';
%     monkey = 'pele';
% 
%     % Specify which data you want to run this on
%     % des_data = 'baseline';
%     des_data = 'post';
%     data_used = 'all';

% Define and set analysis objectives and specifics:
exp_desc = 'Coherence';
roi_string = 'V1 and V4'; % Set for plotting


%% Add paths
set_paths(des_data);

ft_defaults

%% Run script to get ROIs specific to the monkey
vbm_rois;
% Define region of interest 1
roi = V1;
%Define region of interest 2
roi_2 = V4;
roi_combine = [roi roi_2];

%% Configuration parameters and bipolar merging
cfg = config_setup(roi, monkey, des_data, data_used, 1);

% cfg = config_setup(roi, monkey, des_data, data_used);
mcfg = cfg;

[dataAttIn,dataAttOut] = merge_data2(mcfg, true);

%% Select data for respective ROIs
% ROI One
roi_cfg = [];
roi_cfg.channel = roi;
% ROI Two
roi_cfg2 = [];
roi_cfg2.channel = roi_2;
% All data in one file + indices for two ROIs
all_cfg = [];
all_cfg.channel = roi_combine;

% Changing shit for non bilpolar version

v1_AttIn =  ft_selectdata(roi_cfg, dataAttIn);
v4_AttIn =  ft_selectdata(roi_cfg2, dataAttIn);
all_AttIn =  ft_selectdata(all_cfg, dataAttIn);
    

clear dataAttIn

if ~cfg.bsln
    v1_AttOut = ft_selectdata(roi_cfg, dataAttOut);
    v4_AttOut = ft_selectdata(roi_cfg2, dataAttOut);
    all_AttOut = ft_selectdata(all_cfg, dataAttOut);
end

clear dataAttOut

disp('-------------Data Selection done---------------')
disp('-------------Save Begin---------------')

%% Save

save_dir = sprintf('C:/Users/Jordan/Documents/cls_thesis/neuro_thesis/data/monkey_%s_data/bipolar_lowpass_%s_data', monkey, des_data);

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% mk = monkey(1:end-2);

save(append(save_dir, sprintf('/%s_%s_all_AttIn', monkey, des_data(1))), 'all_AttIn')
% Save data for input to model
save(append(save_dir, sprintf('/%s_%s_v1_AttIn', monkey, des_data(1))), 'v1_AttIn')
save(append(save_dir, sprintf('/%s_%s_v4_AttIn', monkey, des_data(1))), 'v4_AttIn')

clear all_AttIn all_AttIn all_AttIn

if ~cfg.bsln
    save(append(save_dir, sprintf('/%s_%s_all_AttOut', monkey, des_data(1))), 'all_AttOut')
    save(append(save_dir, sprintf('/%s_%s_v1_AttOut', monkey, des_data(1))), 'v1_AttOut')
    save(append(save_dir, sprintf('/%s_%s_v4_AttOut', monkey, des_data(1))), 'v4_AttOut')
end

clear all_AttOut all_AttOut all_AttOut

disp('-------------Done---------------')