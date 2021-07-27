%% Set up paths
addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\gc_hierarchies\')
%addpath('D:/MATLAB/mvgc_v1.0')
startup
addpath('C:\Users\Jordan\Documents\cls_thesis\matlab\fieldtrip-20210411')
ft_defaults

format short;
clear all;
close all;clc;
ptic('starting\n')

%% Load data and restructure
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_baseline_data\kurt_b_v1_AttIn.mat')
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_baseline_data\kurt_b_v4_AttIn.mat')

% Load the full pele data, this seems to work fine
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v1_AttIn.mat')
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v4_AttIn.mat')

% Load the full kurt data
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_v1_AttIn.mat')
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_v4_AttIn.mat')

% Get the trial part of the structure with the data
v1_in = v1_AttIn.trial;
v4_in = v4_AttIn.trial;

% Concat them along the 3rd dim so we have the data in the structure we
% want.  Structure needed is (#regions x #obs x #trials)
v1_in = cat(3,v1_in{:});
v4_in = cat(3,v4_in{:});

[max_val, loc] = max(v1_in, [], [1,2], 'linear');
[min_val, loc] = min(v1_in, [], [1,2], 'linear');

max_val = squeeze(max_val);
min_val = squeeze(min_val);

examine = [];

for i = 1:length(max_val)
    if (max_val(i) > 100) | (min_val(i) < -100)
        examine(:,:,end+1) = v1_in(:,:,i)
    end
end
