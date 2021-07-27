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
v1_in = v1_AttIn.trial;
v4_in = v4_AttIn.trial;

% Concat them along the 3rd dim so we have the data in the structure we
% want.  Structure needed is (#regions x #obs x #trials)
v1_in = cat(3,v1_in{:});
v4_in = cat(3,v4_in{:});

%% Parameters
ntrials   = 417;     % no. trials
nobs      = 502;   % no. obs/ trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 30;     % maximum model order for model order estimation (<= 30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1000;    % sample rate (Hz)
fres      = fs*2;     % frequency resolution (empty for automatic calculation)

max_len_f = 700; % Set up matrix to collect f values
len_long = 0;
len_short = 2 * max_len_f;
% gc_one    = zeros(1,1,max_len_f);
% gc_two    = zeros(1,1,max_len_f);
gc_one    = [];
gc_two    = [];

seed      = 0;      % random seed (0 for unseeded)
verb = false; 

%% Variables to help with counting and data locations.  Not sure if needed
trial_count = 1;
too_short = 0;
sizes = [];

% get the number of channels
num_chan_v1 = size(v1_in, 1);
num_chan_v4 = size(v4_in, 1);

for chan1=1:num_chan_v1
    
%     fprintf('\nCurrent V1 channel: %d\n', chan1)
    
    % Access the selected channel for region 1
    v1_chan = v1_in(chan1,:,:);
    
    for chan4=1:num_chan_v4
        % show where we are
        
        fprintf('\nCurrent channel combination is: %d, %d\n', chan1, chan4)
        % Access the selected channel for region 4
        v4_chan = v4_in(chan4,:,:);
        
        % Concatanate the channels
        region_comp = cat(1, v1_chan, v4_chan);
        
        %% Start the calcs of granger cause
        % VAR - channel combinations 
        [aic,bic,moaic,mobic] = tsdata_to_infocrit(region_comp,momax,icregmode, verb);

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
        [A,SIG] = tsdata_to_var(region_comp,morder,regmode);

        ptoc;
        % Check for failed regression
        assert(~isbad(A),'VAR estimation failed');

        % Autocovariance calculation
        ptic('*** var_to_autocov... ');
        [G,info] = var_to_autocov(A,SIG,acmaxlags);
        ptoc;
        
        try
            var_info(info,true); % report results (and bail out on error)

            ptic('\n*** autocov_to_spwcgc... ');
            f= autocov_to_spwcgc(G,fres);
            ptoc;
            % Check for failed spectral GC calculation
            assert(~isbad(f,false),'spectral GC calculation failed');
            
            % Only collect f of size maz_len_f
            if size(f, 3) >= max_len_f
%                 trial_count = trial_count + 1;
%                 gc_one(trial_count,:) = squeeze(f(1,2,1:max_len_f)); 
%                 gc_two(trial_count,:) = squeeze(f(2,1,1:max_len_f));
                gc_one(end+1,:) = squeeze(f(1,2,1:max_len_f)); 
                gc_two(end+1,:) = squeeze(f(2,1,1:max_len_f));

            else
                too_short = too_short+1;
            end
            
            sizes(end+1) = size(f,3);

        catch % for intstances with unstable VAR root
%             skipped_trials(end+1) = t;
%             fprintf('\nskipping trial %d channel combo V1[%d] V4[%d]',t,v1, v4);
%             fid = fopen(skipped_trials_path,'at');
%             fprintf(fid, '\nskipping trial %d channel combo V1[%d] V4[%d]',t,v1, v4);
%             fclose(fid);
            continue    
        end
    end
end

% average the results and plot
gc_1_ave = mean(gc_one, 1);
gc_2_ave = mean(gc_two, 1);
f_all = vertcat(gc_1_ave, gc_2_ave);
gc_plot(f_all, fs);