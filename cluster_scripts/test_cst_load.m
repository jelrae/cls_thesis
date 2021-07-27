%% test control systems

addpath('/home/shall/matlab/mvgc_v1.0');
addpath('/home/shall/matlab/control');
addpath('/home/shall/scripts/gc_hierarchies/cluster_scripts/');
addpath('/home/cbosman1/matlab/fieldtrip/');
addpath('/home/shall/data/bipolarderivatives');


startup
ft_defaults

format short;
rng(938197);
ptic('/nstarting/n')

load('/home/shall/data/bipolarderivatives/ku_v1_attin.mat')
load('/home/shall/data/bipolarderivatives/ku_v4_attin.mat')

choose_resampling = 'yes';

   
resample_cfg = [];
resample_cfg.resamplefs = 250;
resample_cfg.detrend = 'yes';
resample_cfg.demean = 'yes';
[ku_v1_attin] = ft_resampledata(resample_cfg, ku_v1_attin);
[ku_v4_attin] = ft_resampledata(resample_cfg, ku_v4_attin);
v1_in_all = ku_v1_attin.trial;
v1_in_tmp = cell2mat(v1_in_all);
v1_in = v1_in_tmp(1:2,1:1261);

v4_in_all = ku_v4_attin.trial;
v4_in_tmp = cell2mat(v4_in_all);
v4_in = v4_in_tmp(1:2,1:1261);

%% Parameters
ntrials   = 10;     % no. trials
nobs      = 126;   % no. obs/ trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

% morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 30;     % maximum model order for model order estimation (<= 30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 250;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

max_len_f = 800; % Set up matrix to collect f values
len_long = 0;
len_short = 2 * max_len_f;
gc_temp_one    = zeros(1,1,max_len_f);
gc_temp_two    = zeros(1,1,max_len_f);

chan_v1 = size(v1_in, 1);
chan_v4 = size(v4_in, 1);
v1_in_combo = zeros((ntrials*chan_v1*chan_v4),nobs);
v4_in_combo = zeros((ntrials*chan_v1*chan_v4),nobs);
gc_one    = zeros(1,1,max_len_f);
gc_two    = zeros(1,1,max_len_f);

seed      = 0;      % random seed (0 for unseeded)
verb = true; 
%% Variables to help with counting & keeping track of loops, matrix rows etc
count_v4v1 = 0;
skipped_trials = []; 
loop_cnt = 0;
trial_count = 0;
row_count = 0;
too_short = 0;
division_mean = 0;

len_trials = size(ntrials, 1);
dim_ = chan_v1*chan_v4*ntrials;

%% Create matrix of all trial + channel combinations
for v1=1:size(chan_v1,1)
    % Set temp variables
    
    
    for v4=1:chan_v4
        for t=1:ntrials
            trial_count = trial_count  + 1;
            trial_index = t * nobs;
            
            % Create datasets in both directions
            v1_in_combo(t,:) = v1_in(v1,(trial_index - nobs + 1): trial_index);
            v4_in_combo(t,:) = v4_in(v4,(trial_index - nobs + 1): trial_index);
        end
    end
end
disp(v1_in_combo);
%% Loop through channels of roi 1 and 2 and all trials
  

for row=1:size(v1_in_combo,1)
    
    morder    = 'AIC';
    % Create datasets in both directions
    tmp_v1 = v1_in_combo(row,:);
    tmp_v4 = v4_in_combo(row, :);
    v4v1 = vertcat(tmp_v1, tmp_v4);
    input = cat(3, v4v1, 1*ones(size(v4v1)));
    count_v4v1 = count_v4v1 + 1;
    
    disp(input);
    % VAR - channel combinations 
    [aic,bic,moaic,mobic] = tsdata_to_infocrit(input,momax,icregmode, verb);
end
ptoc('/end/n')