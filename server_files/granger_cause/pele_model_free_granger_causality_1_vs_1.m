% %% Set up paths local linux
% addpath('/home/jordan/neuro_thesis/cls_thesis/gc_hierarchies/')
% addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
% %addpath('D:/MATLAB/mvgc_v1.0')
% startup
% addpath('/home/jordan/common/matlab/fieldtrip-20210411')
% ft_defaults
% 
% format short;
% clear all;
% close all;clc;
% ptic('starting\n')
% 
% load('/home/jordan/neuro_thesis/data/monkey_pele_data/bipolar_post_data/pele_p_all_AttIn.mat')
% monkey = 'pele';

% %% Set up paths local windows
% addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\cls_thesis\gc_hierarchies\')
% addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\cls_thesis\helper_functions\')
% %addpath('D:/MATLAB/mvgc_v1.0')
% startup
% addpath('C:\Users\Jordan\Documents\cls_thesis\matlab\fieldtrip-20210411')
% ft_defaults
% 
% format short;
% clear all;
% close all;clc;
% ptic('starting\n')
% 
% load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_all_AttIn.mat')
% monkey = 'kurt';

%% Set up paths Server
addpath('/home/12297127/cls_thesis/gc_hierarchies/')
addpath('/home/12297127/cls_thesis/helper_functions/')
addpath('/home/12297127/matlab/MVGC');
addpath('/home/cbosman1/matlab/fieldtrip/');

startup
ft_defaults

format short;
clear all;
close all;clc;
ptic('starting\n')

monkey = 'pele';
load(sprintf('/home/12297127/data/no_bad_channels/%s_p_all_AttIn.mat', monkey))
load(sprintf('/home/12297127/data/no_bad_channels/%s_p_all_AttOut.mat', monkey))

cfg = [];
cfg.keepsampleinfo = 'yes';
all_attention = ft_appenddata(cfg, all_AttIn, all_AttOut);

clear all_AttIn all_AttOut

%% Parameters

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default) OG LWR

model_order    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
% morder    = 'AIC'; 
momax     = 30;     % maximum model order for model order estimation (<= 30)

acmaxlags = '';   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1000;    % sample rate (Hz)

% fres      = fs/nobs; % for FFT defined as fs/N samples
fres      = 1000;
fnq       = fs/2;

seed      = 0;      % random seed (0 for unseeded)
verb = false; 

%% Variables to help with counting and data locations.  Not sure if needed
sizes = [];
errors = [];
bad_combinations = [];

%% Storage for all GC

gc_forward = {};
gc_backward = {};

%% Identify the regions
fig_6_ROIS;

if strcmp(monkey, 'kurt')
    monkey_caps = 'Kurt';
else
    monkey_caps = 'Pele';
end

bad_channels = {};

%% Do GC
for r1 = 1 : length(regions)-1
    % Get the regions 
    roi_cfg = [];
    roi_cfg.channel = regions{r1};
    region1 = ft_selectdata(roi_cfg, all_attention);
    region1 = region1.trial;
    region1 = cat(3,region1{:});
    num_r1 = size(region1, 1);
    % Cycle over the channels for each of the two regions
    for r2 = r1+1 : length(regions)
        %Reset the gc storing arrays
        gc_one    = [];
        gc_two    = [];
        bad_chans = cell(0,2);
        % Get the regions 
        roi_cfg = [];
        roi_cfg.channel = regions{r2};
        region2 = ft_selectdata(roi_cfg, all_attention);
        region2 = region2.trial;
        region2 = cat(3,region2{:});
        num_r2 = size(region2, 1);
        
        % Check if the regions are the same, if they are, skip!
        for chan1 = 1 : num_r1
            for chan2 = 1 : num_r2
                try
                    fprintf('\nCurrent channel combination for regions %s %s is: %d, %d\n', region_names(r1), region_names(r2), chan1, chan2)
                    region_comp = cat(1, region1(chan1,:,:), region2(chan2,:,:));

                    f_res_bucket_width = fs/size(region_comp,2);
                    fres = floor(fnq/f_res_bucket_width);
                    %% Start the calcs of granger cause
                    % Find the cpsd
                    cpsd_trial = tsdata_to_cpsd(region_comp, fres,'MT', size(region_comp, 2), floor(size(region_comp, 2)/2),3, 4);

                    % find the autocov

                    [G,q] = cpsd_to_autocov(cpsd_trial);
                catch 
                    fprintf('\nCurrent channel combination is of regions %s %s: %d, %d, was not pos def for chol!\n', region_names(r1), region_names(r2), chan1, chan2)
                    bad_combinations(end+1,:) = [r1 r2 chan1 chan2];
                end

                try
    %                 var_info(info,true); % report results (and bail out on error)

                    f= autocov_to_spwcgc(G,fres);
                    % Check for failed spectral GC calculation
                    assert(~isbad(f,false),'spectral GC calculation failed');

                    if ~isreal(f)
                        fprintf('\nCurrent channel combination is complex %s, %s: %d, %d, was not pos def for chol!\n', region_names(r1), region_names(r2), chan1, chan2)
                    end

%                     %Old version
%                     gc_one(end+1,:) = squeeze(f(1,2,:)); 
%                     gc_two(end+1,:) = squeeze(f(2,1,:));
                    % Updated from source code first dim is to, second is
                    % from
                    
                    if isreal(sqrt(squeeze(f(1,2,:)))) && isreal(sqrt(squeeze(f(2,1,:))))

                        % Only collect f of size maz_len_f
                        gc_one(end+1,:) = squeeze(f(1,2,:)); 
                        gc_two(end+1,:) = squeeze(f(2,1,:));

                        sizes(end+1) = size(f,3);

                    else
                        bad_chans(end+1,:) = {chan1 chan2};
                    end

                catch % for intstances with unstable VAR root

                    fprintf('\nCurrent channel combination of regions %s, %s: %d, %d, didnt work!\n', region_names(r1), region_names(r2), chan1, chan2)
                    errors(end+1,:) = [r1 r2 chan1 chan2];
                    continue   
                end 
            end
        end
        gc_forward{end+1} = gc_one;
        gc_backward{end+1} = gc_two;
        bad_channels{end+1} = bad_chans;
    end
end    

clear all_AttIn;
save(sprintf('/home/12297127/cls_thesis/server_files/results/%s_all_attention_model_free_gc_one_v_one', monkey));