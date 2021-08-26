%% Set up and load data
addpath('/home/12297127/matlab/MVGC');
addpath('/home/cbosman1/matlab/fieldtrip/');

ft_defaults
startup

format short;
clear all;
close all;clc;

default_home = '/home/12297127/';

%%Load data

load('/home/12297127/data/no_bad_channels/pele_p_v1_AttIn.mat')
load('/home/12297127/data/no_bad_channels/pele_p_v4_AttIn.mat')

v1_in = v1_AttIn.trial;
v4_in = v4_AttIn.trial;

v1_in = cat(3,v1_in{:});
v4_in = cat(3,v4_in{:});

clear v1_AttIn v4_AttIn

%% Parameters
ntrials   = size(v1_in,3);     % no. trials
nobs      = size(v1_in,2);   % no. obs per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

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

        if strcmpi(model_order,'AIC')
            morder = moaic;
            fprintf('\nusing AIC best model order = %d\n',morder);
        elseif strcmpi(model_order,'BIC')
            morder = mobic;
            fprintf('\nusing BIC best model order = %d\n',morder);
        else
            fprintf('\nusing specified model order = %d\n',morder);
            morder = model_order;
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
            
            if ~isreal(f)
                fprintf('\nCurrent channel combination is complex: %d, %d, didnt work!\n', chan1, chan4)
            end
            
            % Only collect f of size maz_len_f
            gc_one(end+1,:) = squeeze(f(1,2,:)); 
            gc_two(end+1,:) = squeeze(f(2,1,:));
            
            sizes(end+1) = size(f,3);

        catch % for intstances with unstable VAR root
            
            fprintf('\nCurrent channel combination is: %d, %d, didnt work!\n', chan1, chan4)
            errors(end+1,:) = [chan1 chan4];
            continue    
        end
    end
end

% average the results and plot
gc_1_ave = mean(gc_one, 1);
gc_2_ave = mean(gc_two, 1);

%% do some plotting
x_range = (1:1:length(gc_1_ave))./length(gc_1_ave);
x_range = x_range.*fnq;
blue=[.3 .8 .9];red=[1 .1 .1];

plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
set(gca,'FontSize',24,'LineWidth',5,'TickLength',[0.02 0.02])
set(gca,'box','off');legend('FF','FB');legend boxoff;
xlim([0 fnq]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
ylim([0 zgc]);title('V1-V4');set(gca,'Layer','top');
ylabel('Granger causality');xlabel('Frequency (Hz)');