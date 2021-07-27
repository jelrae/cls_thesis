%% Set up and load data
addpath('/home/shall/matlab/mvgc_v1.0');
addpath('/home/cbosman1/matlab/fieldtrip/');
addpath('/home/shall/scripts/gc_hierarchies/cluster_scripts/Copy_of_gc_plot.m');
addpath('/home/shall/matlab/control');

ft_defaults
startup

format short;
clear all;
close all;clc;
rng(938197);
ptic('/nstarting/n')
%% Load data
load('/home/shall/data/bipolarderivatives/ku_v1_attin.mat')
load('/home/shall/data/bipolarderivatives/ku_v4_attin.mat')

%% Create paths to collect trial information
% Experiment description
roi_1 = 'V1';
roi_2 = 'V4';
data_code = 'ku'; % 'pe'
exp_desc = strcat('Data: ', data_code, ' Areas: ', roi_1, ' to ', roi_2);
date = regexprep(char(datetime('now')),':',''); 
folder_path = '/home/shall/results/trial_info/';

% path to csv file to collect skipped trials
length_f_path = strcat(folder_path, 'length_f_', date, '.txt');
fid = fopen(length_f_path,'w+');
fclose(fid);
% path to text file to collect skipped trials + channel combinations
skipped_trials_path = strcat(folder_path,'skipped_',roi_1, roi_2, '_', date, '.txt');
fid = fopen(skipped_trials_path,'w+');
fclose(fid);
skipped_f_path = strcat(folder_path,'skipped_f_',roi_1, roi_2, '_', date, '.txt');
fid = fopen(skipped_f_path,'w+');
fclose(fid);

% Create downsampling cfg - otherwise use matrix as loaded
choose_resampling = 'yes'; % 'no'

if strcmpi(choose_resampling,'yes')    
    resample_cfg = [];
    resample_cfg.resamplefs = 250;
    resample_cfg.detrend = 'yes';
    resample_cfg.demean = 'yes';
    [ku_v1_attin] = ft_resampledata(resample_cfg, ku_v1_attin);
    [ku_v4_attin] = ft_resampledata(resample_cfg, ku_v4_attin);
    v1_in_all = ku_v1_attin.trial;
    v1_in_tmp = cell2mat(v1_in_all);
%     v1_in = v1_in_tmp(1,1:103);
    v1_in = v1_in_tmp;
    v4_in_all = ku_v4_attin.trial;
    v4_in_tmp = cell2mat(v4_in_all);
    v4_in = v4_in_tmp;
%     v4_in = v4_in_tmp(1,1:103);
end

%% Parameters
ntrials   = 1235;     % no. trials
nobs      = 126;   % no. obs/ trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 30;     % maximum model order for model order estimation (<= 30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 250;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

max_len_f = 350; % Set up matrix to collect f values
len_long = 0;
len_short = 2 * max_len_f;
gc_temp_one    = zeros(1,1,max_len_f);
gc_temp_two    = zeros(1,1,max_len_f);

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
chan_v1 = size(v1_in, 1);
chan_v4 = size(v4_in, 1);
len_trials = size(ntrials, 1);
dim_ = chan_v1*chan_v4*ntrials;

%% Loop through channels of roi 1 and 2 and all trials
for v1=1:chan_v1
    for v4=1:chan_v4
        for t=1:ntrials
            trial_count = trial_count  + 1;
            trial_index = t * nobs;
            
            % Create datasets in both directions
            tmp_v1 = v1_in(v1,(trial_index - nobs + 1): trial_index);
            tmp_v4 = v4_in(v4,(trial_index - nobs + 1): trial_index);
            v4v1 = vertcat(tmp_v1, tmp_v4);
            input = cat(3, v4v1, 1*ones(size(v4v1)));
            count_v4v1 = count_v4v1 + 1;
            
            % VAR - channel combinations 
            [aic,bic,moaic,mobic] = tsdata_to_infocrit(input,momax,icregmode, verb);
            morder    = 'AIC';
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
            try
                var_info(info,true); % report results (and bail out on error)


                ptic('\n*** autocov_to_spwcgc... ');
                f= autocov_to_spwcgc(G,fres);
                ptoc;
                % Check for failed spectral GC calculation
                assert(~isbad(f,false),'spectral GC calculation failed');
                
                % Only collect f of size maz_len_f
                if size(f, 3) >= max_len_f
                    gc_one = gc_one + f(1,2,1:max_len_f);
                    gc_two = gc_two + f(2,1,1:max_len_f);

                    division_mean = division_mean + 1;
                    
                elseif size(f, 3) < max_len_f % Too short to be saved
                    fid = fopen(skipped_f_path,'at');
                    fprintf(fid, '\nNot saving f from trial %d channel combo V1[%d] V4[%d]',t,v1, v4);
                    fclose(fid);
                end
           
           
            catch % for intstances with unstable VAR root
                
                continue                 
            end 
        end              
    end
end
ptoc('\nstopping\n')
fprintf('\nLongest f: %d \nShortest f: %d\n',len_long, len_short)
%% Calculate average
save('/home/shall/results/plots/f_all_channels_serial_0303.mat', 'f');
gc_one_average = gc_one(:,:,:) / division_mean;
gc_two_average = gc_two(:,:,:)  / division_mean;
%% Plot granger causality

f_all = [gc_one_average; gc_two_average];
copy_gc_plot(f_all, fs); 
saveas(gcf,'/home/shall/results/plots/v1v4_serial_0303.fig');