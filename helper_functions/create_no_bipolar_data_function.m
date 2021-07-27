function create_no_bipolar_data_function(monkey, des_data, data_used, bipolar)
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
cfg = config_setup(roi, monkey, des_data, data_used);
cfg.COI = load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\kurt_montage_layout\kurt_COI_alongfiner2.mat');
mcfg = cfg;

[dataAttIn,dataAttOut] = merge_data2(mcfg, bipolar);

%% Select data for respective ROIs
% ROI One
roi_cfg = [];
roi_cfg.channel = roi;
% ROI Two
roi_cfg2 = [];
roi_cfg2.channel = roi_2;
% All data in one file + indices for two ROIs

COI = cfg.COI;
COI = struct2cell(COI);
COI = COI{1};
COI = COI(1:50);
all_cfg = [];
all_cfg.channel = permute(COI,[2,1]);
% Changing shit for non bilpolar version

all_AttIn =  ft_selectdata(all_cfg, dataAttIn);

clear dataAttIn

if ~cfg.bsln
    all_AttOut = ft_selectdata(all_cfg, dataAttOut);
end

clear dataAttOut

disp('-------------Data Selection done---------------')
disp('-------------Save Begin---------------')

%% Save
if bipolar
    save_dir = sprintf('C:/Users/Jordan/Documents/cls_thesis/neuro_thesis/data/monkey_%s_data/no_bipolar_%s_data', monkey, des_data);
else
    save_dir = sprintf('C:/Users/Jordan/Documents/cls_thesis/neuro_thesis/data/monkey_%s_data/no_bipolar_%s_data', monkey, des_data);
end
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% mk = monkey(1:end-2);

save(append(save_dir, sprintf('/%s_%s_all_AttIn', monkey, des_data(1))), 'all_AttIn')

clear all_AttIn all_AttIn all_AttIn

if ~cfg.bsln
    save(append(save_dir, sprintf('/%s_%s_all_AttOut', monkey, des_data(1))), 'all_AttOut')
end

clear all_AttOut all_AttOut all_AttOut

disp('-------------Done---------------')