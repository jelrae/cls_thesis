function out = InfoGathering_Neurons_Ramp_to_choice(dosave,doproj)

% out = InfoGathering_Neurons_Ramp_to_choice(dosave)
%
% averages neuron responses in different bins of reaction time, and with
% different reaction times
%
% output structure:
%     out.av_RespRT_LR_even (units*timebins*RTbins*2) - response-locked firing rate on
%       even trials, split by left choices (:,:,:,1) and right choices (:,:,:,2)
%     out.av_RespRT_LR_odd  - as above, for odd trials
%     out.av_RespRT_even (units*timebins*RTbins) - as above, but averaged across left/right choices
%     out.av_RespRT_odd - as above, for odd trials
%     out.timebins_ResponseLong - vector of timestamps for response-locked data
%
%     out.av_Cue1RT_LR_even - (units*timebins*RTbins*2) - cue1-locked firing rate on
%       even trials, split by left choices (:,:,:,1) and right choices (:,:,:,2)
%     out.av_Cue1RT_LR_odd - as above, for odd trials
%     out.av_Cue1L_even (units*timebins*RTbins) - as above, but averaged across left/right choices
%     out.av_Cue1L_odd - as above, for odd trials
%     out.timebins_Cue1Long - vector of timestamps for cue-locked data
%
%     out.Cue1_ON_list - list of times when Cue1 was presented, for all units/trials
%     out.Cue2_ON_list - list of times when Cue2 was presented, for all units/trials
%     out.Cue3_ON_list - list of times when Cue3 was presented, for all units/trials
%     out.Cue4_ON_list - list of times when Cue4 was presented, for all units/trials
%     out.Response_time_list - list of times when response was made, for all units/trials
%
% if dosave==1, saves out to neuronal_regression_results/info_gathering_ramp_to_choice.mat

%% save results after calculating averages?
if nargin<1 
    dosave = 0;
end

%% do random subsamples of trials requested by reviewer 3?

if nargin<2
    doproj = 0;
end



%% find pre-processed data files

bd = get_basedir;
addpath(genpath(bd));
rmpath(genpath(fullfile(bd,'nish_scripts')));
datadir = fullfile(bd,'neuronal_data','LTH_processed_data');
filenames=dir(fullfile(datadir,'Unit*.mat'));

if doproj %have to have run Info_gathering_Neurons_Cue123... first, to store which trials to sample
    %note that we load randproj_trials_inv because we want to use the half
    %of trials that WEREN'T originally used for regression.
    load(fullfile(bd,'neuronal_regression_results',...
        'info_gathering_cue123_alltrials_attentional_final.mat'),'randproj_trials','nProj');
end

%% these are going to become lists of *all* Cue1_ON, Cue2_ON, etc. times, so these latencies can be plotted in figure
Cue1_ON_list = [];
Cue2_ON_list = [];
Cue3_ON_list = [];
Cue4_ON_list = [];
Response_time_list = [];

%% loop over units
for u=1:length(filenames)
    fprintf('Unit %0.0f/%0.0f\n',u,length(filenames)); %display current neuron
    
    %% load in data
    fname = fullfile(datadir,sprintf('Unit%04.0f.mat',u));
    load(fname);

    %% build a 'long' response matrix - from 3s pre-response to 500ms post-response
    Response_Matrix_Long = neuron_matrix(3000,Response_time,500,SpikeTable);
    Response_Matrix_LongRaw = Response_Matrix_Long/1000; %raster, with 1 for spike and 0 for no spike
    Response_Matrix_Long=(Response_Matrix_Long-men)/sdev; %normalised raster
    
    %% build a 'long' Cue1 Matrix - from 500ms pre-Cue1 to 3s post-Cue1
    Cue1_Matrix_Long = neuron_matrix(500,Cue1_ON,3000,SpikeTable);
    Cue1_Matrix_LongRaw = Cue1_Matrix_Long/1000; %raster, with 1 for spike and 0 for no spike
    Cue1_Matrix_Long=(Cue1_Matrix_Long-men)/sdev; %normalised raster
   
    %% smooth these long matrices
    Response_Matrix_smoothed = conv2(1,ones(1,200)./200,Response_Matrix_Long,'valid');
    timebins_ResponseLong = -2900:401; %because original timebins went from -3000 to 500, 'valid' timebins are these
    Cue1_Matrix_smoothed = conv2(1,ones(1,200)./200,Cue1_Matrix_Long,'valid');
    timebins_Cue1Long = -400:2901; %because original timebins went from -500 to 3000, 'valid' timebins are these
    
    %% number of trials/timebins (note that this is the same for both 'long' matrices)
    nTr = size(Response_Matrix_smoothed,1);
    nTimebins = size(Response_Matrix_smoothed,2);
    
    %% subsample smoothed data
    slide_width=10;
    BINS=[1:slide_width:nTimebins];
    
    Response_Matrix_smoothed = Response_Matrix_smoothed(:,BINS);
    timebins_ResponseLong = timebins_ResponseLong(BINS);
    Cue1_Matrix_smoothed = Cue1_Matrix_smoothed(:,BINS);
    timebins_Cue1Long = timebins_Cue1Long(BINS);

    %% calculate 'overall trial' reaction time, and divide into 5 bins
    rt = Response_time-Cue1_ON;
    rtbins{1} = rt<1200;
    rtbins{2} = rt>=1200&rt<1400;
    rtbins{3} = rt>=1400&rt<1600;
    rtbins{4} = rt>=1600&rt<1800;
    rtbins{5} = rt>=1800;

    %% split into odd/even trials, to avoid double-dipping
    odd_trials = logical(zeros(1,nTr));odd_trials(1:2:end) = 1;
    even_trials = logical(zeros(1,nTr));even_trials(2:2:end) = 1;
    
    %% now average different subsets of conditions
    for i = 1:5 %loop over different RT bins
        av_RespRT_even(u,:,i) = mean(Response_Matrix_smoothed(even_trials&rtbins{i},:));
        av_RespRT_odd(u,:,i) = mean(Response_Matrix_smoothed(odd_trials&rtbins{i},:));
        av_RespRT_LR_even(u,:,i,1) = mean(Response_Matrix_smoothed(even_trials&Left_choice==1&rtbins{i},:));
        av_RespRT_LR_even(u,:,i,2) = mean(Response_Matrix_smoothed(even_trials&Left_choice==0&rtbins{i},:));
        av_RespRT_LR_odd(u,:,i,1) = mean(Response_Matrix_smoothed(odd_trials&Left_choice==1&rtbins{i},:));
        av_RespRT_LR_odd(u,:,i,2) = mean(Response_Matrix_smoothed(odd_trials&Left_choice==0&rtbins{i},:));
        
        av_Cue1L_even(u,:,i) = mean(Cue1_Matrix_smoothed(even_trials&rtbins{i},:));
        av_Cue1L_odd(u,:,i) = mean(Cue1_Matrix_smoothed(odd_trials&rtbins{i},:));
        av_Cue1RT_LR_even(u,:,i,1) = mean(Cue1_Matrix_smoothed(even_trials&Left_choice==1&rtbins{i},:));
        av_Cue1RT_LR_even(u,:,i,2) = mean(Cue1_Matrix_smoothed(even_trials&Left_choice==0&rtbins{i},:));
        av_Cue1RT_LR_odd(u,:,i,1) = mean(Cue1_Matrix_smoothed(odd_trials&Left_choice==1&rtbins{i},:));
        av_Cue1RT_LR_odd(u,:,i,2) = mean(Cue1_Matrix_smoothed(odd_trials&Left_choice==0&rtbins{i},:));
        
        if doproj % do random projections requested by reviewer 3?
            for pp = 1:nProj
                included_trials = logical(ones(1,nTr));
                
                %find trials where npv>2 (as these were used to define
                %randproj_trials, and then remove the trials that we used
                %in this projections for the regression
                npv_GT_2 = find(num_pics_viewed>2);
                exclude_these = npv_GT_2(randproj_trials{u}(:,pp));
                included_trials(exclude_these) = 0;
                
                av_RespRT_proj(u,:,i,pp) = mean(Response_Matrix_smoothed(included_trials&rtbins{i},:));
                av_RespRT_LR_proj(u,:,i,1,pp) = mean(Response_Matrix_smoothed(included_trials&Left_choice==1&rtbins{i},:));
                av_RespRT_LR_proj(u,:,i,2,pp) = mean(Response_Matrix_smoothed(included_trials&Left_choice==0&rtbins{i},:));
                av_Cue1L_proj(u,:,i,pp) = mean(Cue1_Matrix_smoothed(included_trials&rtbins{i},:));
                av_Cue1RT_LR_proj(u,:,i,1,pp) = mean(Cue1_Matrix_smoothed(included_trials&Left_choice==1&rtbins{i},:));
                av_Cue1RT_LR_proj(u,:,i,2,pp) = mean(Cue1_Matrix_smoothed(included_trials&Left_choice==0&rtbins{i},:));
            end
        end
    end
    
    %% build up list of event onset times, for figure
    Cue1_ON_list = [Cue1_ON_list Cue1_ON];
    Cue2_ON_list = [Cue2_ON_list Cue2_ON];
    Cue3_ON_list = [Cue3_ON_list Cue3_ON];
    Cue4_ON_list = [Cue4_ON_list Cue4_ON];
    Response_time_list = [Response_time_list Response_time];
    
end

%% put into output structure
out.av_RespRT_even = av_RespRT_even;
out.av_RespRT_odd = av_RespRT_odd;
out.av_RespRT_LR_even = av_RespRT_LR_even;
out.av_RespRT_LR_odd = av_RespRT_LR_odd;
out.av_Cue1L_even = av_Cue1L_even;
out.av_Cue1L_odd = av_Cue1L_odd;
out.av_Cue1RT_LR_even = av_Cue1RT_LR_even;
out.av_Cue1RT_LR_odd = av_Cue1RT_LR_odd;

out.timebins_ResponseLong = timebins_ResponseLong;
out.timebins_Cue1Long = timebins_Cue1Long;

out.Cue1_ON_list = Cue1_ON_list;
out.Cue2_ON_list = Cue2_ON_list;
out.Cue3_ON_list = Cue3_ON_list;
out.Cue4_ON_list = Cue4_ON_list;
out.Response_time_list = Response_time_list;

%% save, if requested
if dosave
    if doproj
        save(fullfile(bd,'neuronal_regression_results','info_gathering_ramp_to_choice.mat'),'out','av*proj');
    else
        save(fullfile(bd,'neuronal_regression_results','info_gathering_ramp_to_choice.mat'),'out');
    end
    
end