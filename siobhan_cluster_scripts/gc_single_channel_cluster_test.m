%% Set up and load data
addpath('/home/shall/matlab/mvgc_v1.0')
    
startup

format short;
clear all;
close all;clc;
rng(938197);
%% Load data
load('/home/shall/data/bipolarderivatives/ku_v1_attin.mat')
load('/home/shall/data/bipolarderivatives/ku_v4_attin.mat')
Nareas2=2;

%% Parameters
ntrials   = 1;     % no. trials
nobs      = 501;   % no. obs/ trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 3;     % maximum model order for model order estimation (30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1000;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

% gc_one    = zeros(50, 50, []);
% gc_two    = zeros(50, 50, []);

seed      = 0;      % random seed (0 for unseeded)
verb = true; % Confirm
%% Set up input matrices
v1_in_tmp = ku_v1_attin.trial;
v1_in = cell2mat(v1_in_tmp);
v1_channelone = v1_in(1, :);
v4_in_tmp = ku_v4_attin.trial;
v4_in = cell2mat(v4_in_tmp);
v4_channelone = v4_in(1, :);
%% Loop 
len_trials = size(ntrials, 1);
chan_v1 = size(v1_in, 1);
chan_v4 = size(v4_in, 1);
count_v4v1 = 0;
% input = zeros(2,2505,2);
trial_count = 0;


for v1=1:chan_v1
    for v4=1:chan_v4
        % Create datasets in both directions
        tmp_v1 = v1_in(v1,1:501);
        tmp_v4 = v4_in(v4,1:501);
        input = vertcat(tmp_v4, tmp_v1);
        count_v4v1 = count_v4v1 + 1;

%             VAR - channel combinations 
        [aic,bic,moaic,mobic] = tsdata_to_infocrit(input,momax,icregmode, verb);
        fprintf('\nbest model order (AIC) = %d\n',moaic);
        fprintf('best model order (BIC) = %d\n',mobic);

        if strcmpi(morder,'AIC')
            morder = moaic;
            fprintf('\nusing AIC best model order = %d\n',morder);
        elseif strcmpi(morder,'BIC')
            morder = mobic;
            fprintf('\nusing BIC best model order = %d\n',morder);
        else
            fprintf('\nusing specified model order = %d\n',morder);
        end

        ptic('\n*** tsdata_to_var... ');
        [A,SIG] = tsdata_to_var(input,morder,regmode);

        ptoc;
        % Check for failed regression
        assert(~isbad(A),'VAR estimation failed');

        % Autocovariance calculation
        ptic('*** var_to_autocov... ');
        [G,info] = var_to_autocov(A,SIG,acmaxlags);
        ptoc;
        var_info(info,true); % report results (and bail out on error)


        ptic('\n*** autocov_to_spwcgc... ');
        f= autocov_to_spwcgc(G,fres);
        ptoc;
        % Check for failed spectral GC calculation
        assert(~isbad(f,false),'spectral GC calculation failed');

        save('/home/shall/results/f.mat', 'f')

        figure(3)
        plot_spw(f,fs);
        hold on
    end
end

        