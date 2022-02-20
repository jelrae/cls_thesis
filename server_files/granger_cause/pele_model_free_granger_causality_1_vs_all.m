%% Set up paths
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

load('/home/12297127/data/no_bad_channels/pele_p_all_AttIn.mat')
monkey = 'pele';

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
    % Bad channels giving power in 100's
    a7A(6) = [];
    a7A(5) = [];
else
    monkey_caps = 'Pele';
end

%% Do GC
% Get all the regions
roi_cfg = [];
roi_cfg.channel = all_channels;
all_regions = ft_selectdata(roi_cfg, all_AttIn);
all_regions = all_regions.trial;
all_regions = cat(3,all_regions{:});
ntrials   = size(all_regions,3);     % no. trials
nobs      = size(all_regions,2);   % no. obs per trial

num_chan_all = size(all_regions, 1);

figure(1);
% regions(6) = [];
% region_names(6) = [];
% for each region
for region = 1 : length(regions)
    %Reset the gc storing arrays
    gc_one    = [];
    gc_two    = [];
    % Get the regions 
    roi_cfg = [];
    roi_cfg.channel = regions{region};
    current_region = ft_selectdata(roi_cfg, all_AttIn);
    current_region = current_region.trial;
    current_region = cat(3,current_region{:});
    num_chan_cur = size(current_region, 1);
    % Cycle over the channels for each of the two (all regions and current)
    % for testing purposes we are just going to do one of each and try to
    % plot
    for chan_all = 1 : num_chan_all
        for chan_cur = 1 : num_chan_cur
            try
                fprintf('\nCurrent channel combination for region %s is: %d, %d\n', region_names(region), chan_cur, chan_all)
                region_comp = cat(1, current_region(chan_cur,:,:), all_regions(chan_all,:,:));
                
                f_res_bucket_width = fs/size(region_comp,2);
                fres = floor(fnq/f_res_bucket_width);
                %% Start the calcs of granger cause
                % Find the cpsd
                cpsd_trial = tsdata_to_cpsd(region_comp, fres,'MT', size(region_comp, 2), floor(size(region_comp, 2)/2),3, 4);

                % find the autocov

                [G,q] = cpsd_to_autocov(cpsd_trial);

            catch 

                fprintf('\nCurrent channel combination is of region %d: %d, %d, was not pos def for chol!\n', region, chan_cur, chan_all)
                bad_combinations(end+1,:) = [region chan_cur chan_all];
            end
            
            try

                f= autocov_to_spwcgc(G,fres);

                % Check for failed spectral GC calculation
                assert(~isbad(f,false),'spectral GC calculation failed');

                if ~isreal(f)
                    fprintf('\nCurrent channel combination is complex: %d, %d, didnt work!\n', chan_cur, chan_all)
                end


%                     %Old version
%                     gc_one(end+1,:) = squeeze(f(1,2,:)); 
%                     gc_two(end+1,:) = squeeze(f(2,1,:));
                % Updated from source code first dim is to, second is
                % from
                gc_one(end+1,:) = squeeze(f(2,1,:)); 
                gc_two(end+1,:) = squeeze(f(1,2,:));

                sizes(end+1) = size(f,3);

            catch % for intstances with unstable VAR root

                fprintf('\nCurrent channel combination is of region %d: %d, %d, didnt work!\n', region, chan_cur, chan_all)
                errors(end+1,:) = [region chan_cur chan_all];
                continue   
            end 
        end
    end
    
    gc_forward{end+1} = gc_one;
    gc_backward{end+1} = gc_two;

end

clear all_AttIn;
save('/home/12297127/cls_thesis/server_files/results/pele_new_model free_gc_one_v_all');