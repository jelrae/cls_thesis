%% Load data
addpath('D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/f_baseline/area_1');
lfp = load('D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/f_baseline/area_1/baseline_F_3_17_3_4.mat');

%% Plot LFP data
% Load the raw lfp data. 

% The following line is based off the data structure as saved by `get_lfps`
lfp = lfp.baseline_lfp;

fs = 1000; % Hz
T = 1/fs;

tx = (0:length(lfp)-1)/fs;

figure; plot(tx,lfp); ...
    xlabel('Time (s)'), ylabel('Amplitude (uV)'), title('LFP');
