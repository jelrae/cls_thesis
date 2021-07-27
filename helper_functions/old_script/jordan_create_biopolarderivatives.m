%% Create the Polpolar derivatives Preprocessing

%% Start Clear and hparams

clear all
clc

% Select the monkey and define region interest (hyperparamters) 
% monkey = 'kurt';
monkey = 'pele';

% Specify which data you want to run this on
% des_data = 'baseline';
des_data = 'post';
data_used = 'all';

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

%% Configuration parameters and bioplar merging
cfg = config_setup(roi, monkey, des_data, data_used);
mcfg = cfg;

[dataAttIn,dataAttOut] = merge_data2(mcfg);

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

v1_AttIn =  ft_selectdata(roi_cfg, dataAttIn);
v4_AttIn =  ft_selectdata(roi_cfg2, dataAttIn);
all_AttIn =  ft_selectdata(all_cfg, dataAttIn);

if ~cfg.bsln
    v1_AttOut = ft_selectdata(roi_cfg, dataAttOut);
    v4_AttOut = ft_selectdata(roi_cfg2, dataAttOut);
    all_AttOut = ft_selectdata(all_cfg, dataAttOut);
end

disp('-------------Data Selection done---------------')
disp('-------------Save Begin---------------')

%% Save
save_dir = sprintf('/home/jordan/neuro_thesis/data/monkey_%s_data/bioplar_%s_data', monkey, des_data);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% mk = monkey(1:end-2);

save(append(save_dir, sprintf('/%s_all_AttIn', monkey)), 'all_AttIn')
% Save data for input to model
save(append(save_dir, sprintf('/%s_v1_AttIn', monkey)), 'v1_AttIn')
save(append(save_dir, sprintf('/%s_v4_AttIn', monkey)), 'v4_AttIn')

if ~cfg.bsln
    save(append(save_dir, sprintf('/%s_all_AttOut', monkey)), 'all_AttOut')
    save(append(save_dir, sprintf('/%s_v1_AttOut', monkey)), 'v1_AttOut')
    save(append(save_dir, sprintf('/%s_v4_AttOut', monkey)), 'v4_AttOut')
end


%% Save data for input to model
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v1_attin', 'ku_v1_attin')
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v1_attout', 'ku_v1_attout')
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v4_attin', 'ku_v1_attin')
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v4_attout', 'ku_v1_attout')

disp('-------------Done---------------')