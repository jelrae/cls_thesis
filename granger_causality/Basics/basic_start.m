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
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_baseline_data\kurt_b_v1_AttIn.mat')
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_baseline_data\kurt_b_v4_AttIn.mat')

% Get the trial part of the structure with the data
v1_all = v1_AttIn.trial;
v4_all = v4_AttIn.trial;

% Concat them along the 3rd dim so we have the data in the structure we
% want.  Structure needed is (#regions x #obs x #trials)
v1_all = cat(3,v1_all{:});
v4_all = cat(3,v4_all{:});

% Way to permute them if we need to do this at some point
% v1_all = permute(v1_all, [3,1,2]);
% v4_all = permute(v4_all, [3,1,2]);

% Try to get the regions and concat them

v1_1 = v1_all(1,:,:);
v4_1 = v4_all(1,:,:);
region_comp = cat(1, v1_1, v4_1);
