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

% Load the full AttIn pele data, this seems to work fine
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v1_AttIn.mat')
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v4_AttIn.mat')

% Load the full AttOut pele data, not tested yet
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v1_AttOut.mat')
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v4_AttOut.mat')

% Load in pele low pass filtered
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_lowpass_post_data\pele_p_v1_AttIn.mat')
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_lowpass_post_data\pele_p_v4_AttIn.mat')

% Load the full kurt data - Problem child
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_v1_AttIn.mat')
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_v4_AttIn.mat')

% Get the trial part of the structure with the data AttIn
v1_in_preclean = v1_AttIn.trial;
v4_in_preclean = v4_AttIn.trial;

% Get the trial part of the structure with the data AttOut
% v1_in_preclean = v1_AttOut.trial;
% v4_in_preclean = v4_AttOut.trial;

% Concat them along the 3rd dim so we have the data in the structure we
% want.  Structure needed is (#regions x #obs x #trials)

%For when we dont preclean
v1_in = cat(3,v1_in_preclean{:});
v4_in = cat(3,v4_in_preclean{:});

%For when we preclean
% v1_in_preclean = cat(3,v1_in_preclean{:});
% v4_in_preclean = cat(3,v4_in_preclean{:});

% Dirty remove outlier

% [max_val_1, ~] = max(v1_in_preclean, [], [1,2], 'linear');
% [min_val_1, ~] = min(v1_in_preclean, [], [1,2], 'linear');
% 
% max_val_1 = squeeze(max_val_1);
% min_val_1 = squeeze(min_val_1);
% 
% [max_val_4, ~] = max(v4_in_preclean, [], [1,2], 'linear');
% [min_val_4, ~] = min(v4_in_preclean, [], [1,2], 'linear');
% 
% max_val_4 = squeeze(max_val_4);
% min_val_4 = squeeze(min_val_4);
% 
% v1_in = [];
% v4_in = [];
% 
% thresh = 130;
% 
% for i = 1:length(max_val_1)
%     if ((max_val_1(i) < thresh) & (min_val_1(i) > -thresh)) & ((max_val_4(i) < thresh) & (min_val_4(i) > -thresh))
%         v1_in(:,:,end+1) = v1_in_preclean(:,:,i);
%         v4_in(:,:,end+1) = v4_in_preclean(:,:,i);
%     end
% end

% Add some random noise to hope to help...
% v1_in = v1_in +(.25 .* randn(size(v1_in)));
% v4_in = v4_in +(.25 .* randn(size(v4_in)));

clear v1_in_preclean v4_in_preclean
clear max_val_1 max_val_4 min_val_1 min_val_4
clear v1_AttIn v4_AttIn

% v1_in(:,:,1) = [];
% v4_in(:,:,1) = [];

% Reduce number of trials

% v1_in = v1_in(:,:,1:1235);
% v4_in = v4_in(:,:,1:1235);


%% Parameters
ntrials   = size(v1_in,3);     % no. trials
nobs      = size(v1_in,2);   % no. obs per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 30;     % maximum model order for model order estimation (<= 30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1000;    % sample rate (Hz)

% fres      = fs/nobs; % for FFT defined as fs/N samples
fres      = 1000;
fnq       = fs/2;

gc_one    = [];
gc_two    = [];

seed      = 0;      % random seed (0 for unseeded)
verb = false; 

%% Variables to help with counting and data locations.  Not sure if needed
sizes = [];
errors = [];

% get the number of channels
num_chan_v1 = size(v1_in, 1);
num_chan_v4 = size(v4_in, 1);

v1_chan = v1_in(1,:,:);
v4_chan = v4_in(1,:,:);
        
% Concatanate the channels
region_comp = cat(1, v1_chan, v4_chan);
        
%% Start the calcs of granger cause
% Find the cpsd
cpsd_trial_no_noise = tsdata_to_cpsd(region_comp,fres);
region_comp_noise = region_comp + (.25 .* randn(size(region_comp)));
cpsd_trial_noise = tsdata_to_cpsd(region_comp_noise,fres);