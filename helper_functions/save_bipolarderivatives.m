%addpath('D:/MATLAB/forSiobhan/UvA_internship/granger_02112020.m')
addpath('/home/jordan/common/matlab/fieldtrip-20210411')
addpath('/home/jordan/neuro_thesis/gc_hierarchies/helper_functions/')
addpath('D:/MATLAB/forSiobhan/Kurt/Data/')
addpath('/home/jordan/neuro_thesis/data/kurt_montage_layout')
addpath('D:/MATLAB/forSiobhan/UvA_internship/Pele/Data/')
addpath('/home/jordan/neuro_thesis/data/pele_montage_layout')
clear all
clc
ft_defaults

%% Select the monkey and define region interest (hyperparamters)
monkey = 'kurt';
% monkey = 'pele';
%% Define and set analysis objectives and specifics:
exp_desc = 'Coherence';
roi_string = 'V1 and V4'; % Set for plotting
%% Run script to get ROIs specific to the monkey
vbm_rois;
% Define region of interest 1
roi = V1;
%Define region of interest 2
roi_2 = V4;
roi_combine = [roi roi_2];

%% Configuration parameters - the following are consistent
cfg = [];
cfg.equal = 'yes';
cfg.bsln = 0;
cfg.norm = 'yes';
cfg.montage = 'yes'; %here you need to add the name of the montage structure you want to use... alongfinger_v2
cfg.channel = roi;
%% Assigning parameter 'monkey' above sets up the cfg with monkey specific info
if strcmp(monkey,'kurt')
    % Monkey specific information
    monkey_code = 'K';
    cfg.path = 'D:\MATLAB\forSiobhan\UvA_internship\Kurt\Data\';
    cfg.monkey = 'ku';
    cfg.v_session = [40 41 42 43 44 45 46 48]; % Vector referring to number of sessions
%     cfg.v_session = [40 41];    
    cfg.bipolarder = 'D:\MATLAB\forSiobhan\UvA_internship\Kurt\kurt_montage_layout\kurt_montage_bipolar_alongfinger_v2.mat';
    mcfg = cfg;
    [dataAttIn,dataAttOut] = merge_data2(mcfg);
elseif strcmp(monkey,'pele')    
    monkey_code = 'P';
    cfg.path = 'D:\MATLAB\forSiobhan\UvA_internship\Pele\Data\';
    cfg.monkey = 'pe';
    cfg.v_session = [28 29 30 31 32 33 34 35 36 37 38 39 40 42]; % Vector referring to number of sessions
    cfg.bipolarder = 'D:\MATLAB\forSiobhan\UvA_internship\Pele\pele_montage_layout\pele_montage_layout\pele_montage_bipolar_alongfinger_v2.mat';
    mcfg = cfg;
    [dataAttIn,dataAttOut] = merge_data2(mcfg); 
end

%% Select data for respective ROIs
if strcmp(monkey,'kurt')
    % ROI One
    roi_cfg = [];
    roi_cfg.channel = roi;
    ku_v1_attin =  ft_selectdata(roi_cfg, dataAttIn);
    ku_v1_attout = ft_selectdata(roi_cfg, dataAttOut);
    % ROI Two
    roi_cfg2 = [];
    roi_cfg2.channel = roi_2;
    ku_v4_attin =  ft_selectdata(roi_cfg2, dataAttIn);
    ku_v4_attout = ft_selectdata(roi_cfg2, dataAttOut);
    % All data in one file + indices for two ROIs
    all_cfg = [];
    all_cfg.channel = roi_combine;
    ku_all_attin =  ft_selectdata(all_cfg, dataAttIn);
    ku_all_attout = ft_selectdata(all_cfg, dataAttOut);
elseif strcmp(monkey,'pele')
     % ROI One
    roi_cfg = [];
    roi_cfg.channel = roi;
    pe_v1_attin =  ft_selectdata(roi_cfg, dataAttIn);
    pe_v1_attout = ft_selectdata(roi_cfg, dataAttOut);
    % ROI Two
    roi_cfg2 = [];
    roi_cfg2.channel = roi_2;
    pe_v4_attin =  ft_selectdata(roi_cfg2, dataAttIn);
    pe_v4_attout = ft_selectdata(roi_cfg2, dataAttOut);
    % All data in one file + indices for two ROIs
    all_cfg = [];
    all_cfg.channel = roi_combine;
    pe_all_attin =  ft_selectdata(all_cfg, dataAttIn);
    pe_all_attout = ft_selectdata(all_cfg, dataAttOut);
end
%% Save
if strcmp(monkey,'kurt')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_all_attin', 'ku_all_attin')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_all_attout', 'ku_all_attout')
    % Save data for input to model
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v1_attin', 'ku_v1_attin')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v1_attout', 'ku_v1_attout')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v4_attin', 'ku_v4_attin')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v4_attout', 'ku_v4_attout')
elseif strcmp(monkey,'pele')  
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/pe_all_attin', 'pe_all_attin')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/pe_all_attout', 'pe_all_attout')
    % Save data for input to model
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/pe_v1_attin', 'pe_v1_attin')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/pe_v1_attout', 'pe_v1_attout')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/pe_v4_attin', 'pe_v4_attin')
    save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/pe_v4_attout', 'pe_v4_attout')
end
%% Save data for input to model
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v1_attin', 'ku_v1_attin')
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v1_attout', 'ku_v1_attout')
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v4_attin', 'ku_v1_attin')
% save('D:/MATLAB/forSiobhan/UvA_internship/GC-model/GC_practise/ku_v4_attout', 'ku_v1_attout')

