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

load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_all_AttIn.mat')

vbm_rois;