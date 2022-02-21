clear

addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
% set_paths_roi;

addpath('/home/jordan/common/matlab/fieldtrip-20210411')

clear all
clc
ft_defaults
%% set up base directory where data/scripts are stored
[bd] = get_basedir;
addpath(genpath(fullfile(bd)));

datadir = fullfile(bd,'neuronal_data');
filenames=dir(fullfile(datadir,'*.mat'));

load('frank_area');frank_area=area_index; clear area_index;
load('miles_area');

%% Testing things for making it in fieldtrip
data.label = {18;20};
data.fsample = 1000;

%1 row per channel, with following structure:
% column 1 = recording session number
% column 2 = channel ID
% column 3 = area index (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
% column 4 = hemisphere (1 = left hemisphere; 2 = right hemisphere)
% column 5 = anterior-posterior coordinate, stereotactic space (mm)
% column 6 = X-coordinate (medial-lateral position), chamber space (mm)
% column 7 = number of rotations of microdrive used to reach recording site (each rotation = 1/3 mm lowering of electrode)
% column 8 = ??
% column 9 = mm electrode travelled to reach recording site (mm, column 7 divided by 3)
% column 10= ??
% column 11= length of guide tube (mm)
% column 12= which session of a given stimulus set (0 = learning day
%            when stimulus values were learnt, not information sampling session;
%            1 = first day with this stimulus set in information sampling task;
%            2= second day; 3 = third day);

%% main loop over each unit from information gathering experiment

new_correct=100; %this is the minimum length of time for a picture to be considered 'viewed', used below

session_and_channels = [];
sessions_F = {};
sessions_M = {};
length(filenames)
% Variables for the storing of information
ACC_channels = {};
DLPFC_channels = {};
OFC_channels = {};
channel_names = {};
session_number = 0;
tmp_data_holder1 = [];
tmp_data_holder2 = [];
tmp_data_holder_pre_choice_300 = [];
tmp_data_holder_pre_choice_500 = [];
tmp_data_holder_pre_feedback = [];
tmp_data_holder_first_fix = [];
tmp_data_holder_finish_fix = [];


for u=1:length(filenames) %u = 'unit'

    %% get information about session/channel from filename
    sb=filenames(u).name(1); %subject name, F=frank, M=miles
    % change name(2:5) to 2:4
    snum=str2num(filenames(u).name(2:4)); %session number, 1-49 for Miles, 1-36 for Frank
    chnum=str2num(filenames(u).name(strfind(filenames(u).name,'C')+1:strfind(filenames(u).name,'U')-2)); %channel number

    %% load in the data
    load(filenames(u).name);
    u;
    filenames(u).name
    sstr = append('S', filenames(u).name(2:4));
    cstr = append('C', filenames(u).name(strfind(filenames(u).name,'C')+1:strfind(filenames(u).name,'U')-2));
    
    if session_number == 0
        session_number = snum;
        session_string = sstr;
        session_monkey = sb;
    end

%           disp('Epoch/pre-process');

   %% isolate information gathering trials where error code doesn't indicate
   %  timeout etc., and therefore where subject completed the trial
   inf_trials=find(BhvInfo.ConditionNumber==3|BhvInfo.ConditionNumber==4);
   trial_error=BhvInfo.TrialError; % 0 = 'Correct Response' (chose side with higher expected value),
   % 6 = 'Incorrect Response' (chose side with lower expected value)
   % 1 = No Response, 2 = Late Response, 3 = Break Fixation, 4 = No Fixation
   % 5 = Early Response, 7 = Lever Break, 8 = Ignored, 9 =Aborted
%    inf_responded_trials=inf_trials(trial_error(inf_trials)==0|trial_error(inf_trials)==6);
   inf_responded_trials=inf_trials(trial_error(inf_trials)==0);
   inf_cond=BhvInfo.ConditionNumber(inf_responded_trials);
   condition=BhvInfo.ConditionNumber(inf_responded_trials); %3 = option trial; 4 = attribute trial

    %% For future use
    %find which row of miles_area/frank_area this corresponds to, and then
    %find out which brain region was recorded from
    if strcmp(sb,'F')
        str_area=find(frank_area(:,1)==snum & frank_area(:,2)==chnum);
        brain_region=frank_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
        if brain_region == 1 || brain_region == 2 || brain_region == 3
            sessions_F(end+1,:) = {sstr};
        end 
    else
        str_area=find(miles_area(:,1)==snum&miles_area(:,2)==chnum);
        brain_region=miles_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
        if brain_region == 1 || brain_region == 2 || brain_region == 3
            sessions_M(end+1,:) = {sstr};
        end 
    end

    % Store which region the current combo is in

    %% To just figure out whats in this bitch
   session_and_channels(end+1,:) = [snum, chnum, brain_region];

%                disp(str_area);
%                disp(brain_region);

    good_lfp = BhvInfo.PlexonFieldPotentials(inf_responded_trials);

    trial_lfp_first_fix = [];
    trial_lfp_finish_fix = [];
    trial_lfp_img1 = [];
    trial_lfp_img2 = [];
    trial_lfp_img3 = [];
    trial_lfp_img4 = [];
    trial_lfp_end = [];
    trial_lfp_pre_choice_300 = [];
    trial_lfp_pre_choice_500 = [];
    trial_lfp_pre_feedback = [];

    for tr = 1:length(inf_responded_trials)
        ct=inf_responded_trials(tr);
%         % pull the codes for the time stamps and the times
        Codes{tr}(:,1)=BhvInfo.CodeNumbers{1,ct};
        Codes{tr}(:,2)=BhvInfo.CodeTimes{1,ct};

        %times when pictures were first fixated, position_times, and were removed, pic_end_times
        %(91=top left, 92 = bottom left, 93 = top right, 94 = bottom right):
%         picture_codes{tr}=BhvInfo.PlexonStrobes{ct}(find(BhvInfo.PlexonStrobes{ct}(:,1)>90&BhvInfo.PlexonStrobes{ct}(:,1)<95),:);
%         position_times{tr}=BhvInfo.CodeTimes{ct}(find(BhvInfo.CodeNumbers{ct}>90&BhvInfo.CodeNumbers{ct}<95));
        pic_end_times{tr}(:,1)=BhvInfo.CodeNumbers{ct}(find(BhvInfo.CodeNumbers{ct}==11|BhvInfo.CodeNumbers{ct}==12|BhvInfo.CodeNumbers{ct}==90|BhvInfo.CodeNumbers{ct}==72|BhvInfo.CodeNumbers{ct}==73));
        pic_end_times{tr}(:,2)=BhvInfo.CodeTimes{ct}(find(BhvInfo.CodeNumbers{ct}==11|BhvInfo.CodeNumbers{ct}==12|BhvInfo.CodeNumbers{ct}==90|BhvInfo.CodeNumbers{ct}==72|BhvInfo.CodeNumbers{ct}==73));
        
        for et = 1:(length(pic_end_times{tr}))
            end_time = pic_end_times{tr}(et,2);
            if pic_end_times{tr}(et,1) == 90
                eval(sprintf('trial_lfp_img%d(end+1,:) = good_lfp{tr}(end_time - 299 :end_time);',et-2));
%                     trial_lfp_img1(end+1,:) = good_lfp{tr}(first_end - 299 :first_end);
            elseif pic_end_times{tr}(et,1) == 11
                trial_lfp_first_fix(end+1,:) = good_lfp{tr}(end_time - 149 :end_time);
            elseif pic_end_times{tr}(et,1) == 12
                trial_lfp_finish_fix(end+1,:) = good_lfp{tr}(end_time - 499 :end_time);
            elseif pic_end_times{tr}(et,1) == 72
                trial_lfp_pre_choice_300(end+1,:) = good_lfp{tr}(end_time - 299 :end_time);
                trial_lfp_pre_choice_500(end+1,:) = good_lfp{tr}(end_time - 499 :end_time);
            elseif pic_end_times{tr}(et,1) == 73
                trial_lfp_pre_feedback(end+1,:) = good_lfp{tr}(end_time - 999 :end_time);
            end
        end
        
%         % Collect the last one before choice.
%         end_time = pic_end_times{tr}(end,2);
%         trial_lfp_pre_choice_300(end+1,:) = good_lfp{tr}(end_time - 299 :end_time);

    end


    %% formate data and save it in arrays for later use...  Maybe I should just put it into a cell now?

    % todo: I should also store what regions it has, so a list of which 
    % regions are non-empty
    % used_regions = [];
    % if ~isempty(region)
    % used_regions(end+1) = region;
    % ....

    %  if the two are different
    if snum ~= session_number

        var_name = append(session_monkey,'_',session_string, '_data');
        tmp_data1.label = channel_names;
        tmp_data1.fsample = 1000;
        eval(sprintf('%s_ACC_channels = ACC_channels;', var_name));
        eval(sprintf('%s_DLPFC_channels = DLPFC_channels;', var_name));
        eval(sprintf('%s_OFC_channels = OFC_channels;', var_name));
        tmp_data2 = tmp_data1;
        tmp_data_first_fix = tmp_data1;
        tmp_data_finish_fix = tmp_data1;
        tmp_data_choice_300 = tmp_data1;
        tmp_data_choice_500 = tmp_data1;
        tmp_data_pre_feedback = tmp_data1;

        time = 0:(300-1);
        time = time./data.fsample;
        time_150 = 0:(150-1);
        time_150 = time_150./data.fsample;
        time_500 = 0:(500-1);
        time_500 = time_500./data.fsample;
        time_1000 = 0:(1000-1);
        time_1000 = time_1000./data.fsample;
        
        full_time_1 = cell(1,size(tmp_data_holder1,1));
        full_time_2 = cell(1,size(tmp_data_holder2,1));
        full_time_choice_300 = cell(1,size(tmp_data_holder_pre_choice_300,1));
        full_time_choice_500 = cell(1,size(tmp_data_holder_pre_choice_500,1));
        full_time_pre_feedback = cell(1,size(tmp_data_holder_pre_feedback,1));
        full_time_first_fix = cell(1,size(tmp_data_holder_first_fix,1));
        full_time_finish_fix = cell(1,size(tmp_data_holder_finish_fix,1));
        
        cell_img_1_lfp = cell(1,size(tmp_data_holder1,1));
        cell_img_2_lfp = cell(1,size(tmp_data_holder2,1));
        cell_img_pre_choice_300 = cell(1,size(tmp_data_holder_pre_choice_300,1));
        cell_img_pre_choice_500 = cell(1,size(tmp_data_holder_pre_choice_500,1));
        cell_img_pre_feedback = cell(1,size(tmp_data_holder_pre_feedback,1));
        cell_img_first_fix = cell(1,size(tmp_data_holder_first_fix,1));
        cell_img_finish_fix = cell(1,size(tmp_data_holder_finish_fix,1));
        
        for irow = 1:size(tmp_data_holder1,1)
            cell_img_1_lfp{1,irow} = squeeze(tmp_data_holder1(irow,:,:));
            cell_img_2_lfp{1,irow} = squeeze(tmp_data_holder2(irow,:,:));
            cell_img_pre_choice_300{1,irow} = squeeze(tmp_data_holder_pre_choice_300(irow,:,:));
            cell_img_pre_choice_500{1,irow} = squeeze(tmp_data_holder_pre_choice_500(irow,:,:));
            cell_img_pre_feedback{1,irow} = squeeze(tmp_data_holder_pre_feedback(irow,:,:));
            cell_img_first_fix{1,irow} = squeeze(tmp_data_holder_first_fix(irow,:,:));
            cell_img_finish_fix{1,irow} = squeeze(tmp_data_holder_finish_fix(irow,:,:));
            
            full_time_1{1,irow} = time;
            full_time_2{1,irow} = time;
            full_time_choice_300{1,irow} = time;
            full_time_choice_500{1,irow} = time_500;
            full_time_pre_feedback{1,irow} = time_1000;
            full_time_first_fix{1,irow} = time_150;
            full_time_finish_fix{1,irow} = time_500;
        end
        
        tmp_data1.trial = cell_img_1_lfp;
        tmp_data1.time = full_time_1;
        tmp_data2.trial = cell_img_2_lfp;
        tmp_data2.time = full_time_2;
        tmp_data_choice_300.trial = cell_img_pre_choice_300;
        tmp_data_choice_300.time = full_time_choice_300;
        tmp_data_choice_500.trial = cell_img_pre_choice_500;
        tmp_data_choice_500.time = full_time_choice_500;
        tmp_data_pre_feedback.trial = cell_img_pre_feedback;
        tmp_data_pre_feedback.time = full_time_pre_feedback;
        tmp_data_first_fix.trial = cell_img_first_fix;
        tmp_data_first_fix.time = full_time_first_fix;
        tmp_data_finish_fix.trial = cell_img_finish_fix;
        tmp_data_finish_fix.time = full_time_finish_fix;

        % do the preprossesing now

%         p_cfg = [];
%         p_cfg.demean = 'yes';
%         p_cfg.detrend = 'yes';
%         tmp_data1 = ft_preprocessing(p_cfg, tmp_data1);
%         tmp_data2 = ft_preprocessing(p_cfg, tmp_data2);

        eval(sprintf('%s_1 = tmp_data1;', var_name));
        eval(sprintf('%s_2 = tmp_data2;', var_name));
        eval(sprintf('%s_pre_choice_300 = tmp_data_choice_300;', var_name));
        eval(sprintf('%s_pre_choice_500 = tmp_data_choice_500;', var_name));
        eval(sprintf('%s_pre_feedback = tmp_data_pre_feedback;', var_name));
        eval(sprintf('%s_first_fix = tmp_data_first_fix;', var_name));
        eval(sprintf('%s_finish_fix = tmp_data_finish_fix;', var_name));
        
        clear tmp_data1 tmp_data2 tmp_data_choice_300 tmp_data_choice_500 tmp_data_pre_feedback tmp_data_first_fix tmp_data_finish_fix;

        % save things off
        save(sprintf('%s.mat', var_name), sprintf('%s_1', var_name), ...
            sprintf('%s_2', var_name), sprintf('%s_pre_choice_300', var_name), ...
            sprintf('%s_pre_choice_500', var_name), ...
            sprintf('%s_pre_feedback', var_name), ...
            sprintf('%s_first_fix', var_name), ...
            sprintf('%s_finish_fix', var_name), ...
            sprintf('%s_ACC_channels', var_name),  ...
            sprintf('%s_DLPFC_channels', var_name),  ...
            sprintf('%s_OFC_channels', var_name));
        
        clear tmp_data1 tmp_data2 tmp_data_choice_300 tmp_data_choice_500 ...
            tmp_data_pre_feedback tmp_data_holder1 tmp_data_holder2 ...
            tmp_data_holder_pre_choice_300 tmp_data_holder_pre_choice_500 ...
            tmp_data_holder_pre_feedback;

        if brain_region == 1
            channel_names = {append(sstr,cstr)};
            ACC_channels = {append(sstr,cstr)};
            DLPFC_channels = {};
            OFC_channels = {};
        elseif brain_region == 2
            channel_names = {append(sstr,cstr)};
            ACC_channels = {};
            DLPFC_channels = {append(sstr,cstr)};
            OFC_channels = {};
        elseif brain_region == 3
            channel_names = {append(sstr,cstr)};
            ACC_channels = {};
            DLPFC_channels = {};
            OFC_channels = {append(sstr,cstr)};
        else
            channel_names = {};
            ACC_channels = {};
            DLPFC_channels = {};
            OFC_channels = {};
        end

        session_number = snum;
        session_string = sstr;
        session_monkey = sb;
        tmp_data_holder1 = [];
        tmp_data_holder2 = [];
        tmp_data_holder_pre_choice_300 = [];
        tmp_data_holder_pre_choice_500 = [];
        tmp_data_holder_pre_feedback = [];
        tmp_data_holder_first_fix = [];
        tmp_data_holder_finish_fix = [];
        

        % only store if the brain region is 1, 2 or 3
        if brain_region == 1 || brain_region == 2 || brain_region == 3
            tmp_data_holder1(:,end+1,:) = reshape(trial_lfp_img1, [size(trial_lfp_img1,1), 1, size(trial_lfp_img1,2)]);
            tmp_data_holder2(:,end+1,:) = reshape(trial_lfp_img2, [size(trial_lfp_img2,1), 1, size(trial_lfp_img2,2)]);
            tmp_data_holder_pre_choice_300(:,end+1,:) = reshape(trial_lfp_pre_choice_300, [size(trial_lfp_pre_choice_300,1), 1, size(trial_lfp_pre_choice_300,2)]);
            tmp_data_holder_pre_choice_500(:,end+1,:) = reshape(trial_lfp_pre_choice_500, [size(trial_lfp_pre_choice_500,1), 1, size(trial_lfp_pre_choice_500,2)]);
            tmp_data_holder_pre_feedback(:,end+1,:) = reshape(trial_lfp_pre_feedback, [size(trial_lfp_pre_feedback,1), 1, size(trial_lfp_pre_feedback,2)]);
            tmp_data_holder_first_fix(:,end+1,:) = reshape(trial_lfp_first_fix, [size(trial_lfp_first_fix,1), 1, size(trial_lfp_first_fix,2)]);
            tmp_data_holder_finish_fix(:,end+1,:) = reshape(trial_lfp_finish_fix, [size(trial_lfp_finish_fix,1), 1, size(trial_lfp_finish_fix,2)]);
        end 

    % if we hit the end of the filenames
    elseif u == length(filenames)

        % First store the last information
        if brain_region == 1
            channel_names(end+1,1) = {append(sstr,cstr)};
            ACC_channels(end+1) = {append(sstr,cstr)};
        elseif brain_region == 2
            channel_names(end+1,1) = {append(sstr,cstr)};
            DLPFC_channels(end+1) = {append(sstr,cstr)};
        elseif brain_region == 3
            channel_names(end+1,1) = {append(sstr,cstr)};
            OFC_channels(end+1) = {append(sstr,cstr)};
        end

        % only store if the brain region is 1, 2 or 3
        if brain_region == 1 || brain_region == 2 || brain_region == 3
%             disp('hit here')
            tmp_data_holder1(:,end+1,:) = reshape(trial_lfp_img1, [size(trial_lfp_img1,1), 1, size(trial_lfp_img1,2)]);
            tmp_data_holder2(:,end+1,:) = reshape(trial_lfp_img2, [size(trial_lfp_img2,1), 1, size(trial_lfp_img2,2)]);
            tmp_data_holder_pre_choice_300(:,end+1,:) = reshape(trial_lfp_pre_choice_300, [size(trial_lfp_pre_choice_300,1), 1, size(trial_lfp_pre_choice_300,2)]);
            tmp_data_holder_pre_choice_500(:,end+1,:) = reshape(trial_lfp_pre_choice_500, [size(trial_lfp_pre_choice_500,1), 1, size(trial_lfp_pre_choice_500,2)]);
            tmp_data_holder_pre_feedback(:,end+1,:) = reshape(trial_lfp_pre_feedback, [size(trial_lfp_pre_feedback,1), 1, size(trial_lfp_pre_feedback,2)]);
            tmp_data_holder_first_fix(:,end+1,:) = reshape(trial_lfp_first_fix, [size(trial_lfp_first_fix,1), 1, size(trial_lfp_first_fix,2)]);
            tmp_data_holder_finish_fix(:,end+1,:) = reshape(trial_lfp_finish_fix, [size(trial_lfp_finish_fix,1), 1, size(trial_lfp_finish_fix,2)]);
        end

        % todo : Change the ACC cchannels bit to not save them into the tmp
        % data.  
        var_name = append(session_monkey,'_', session_string, '_data');
        tmp_data1.label = channel_names;
        tmp_data1.fsample = 1000;
        eval(sprintf('%s_ACC_channels = ACC_channels;', var_name));
        eval(sprintf('%s_DLPFC_channels = DLPFC_channels;', var_name));
        eval(sprintf('%s_OFC_channels = OFC_channels;', var_name));
        tmp_data2 = tmp_data1;
        tmp_data_choice_300 = tmp_data1;

        time = 0:(300-1);
        time = time./data.fsample;
        time_150 = 0:(150-1);
        time_150 = time_150./data.fsample;
        time_500 = 0:(500-1);
        time_500 = time_500./data.fsample;
        time_1000 = 0:(1000-1);
        time_1000 = time_1000./data.fsample;
        
        full_time_1 = cell(1,size(tmp_data_holder1,1));
        full_time_2 = cell(1,size(tmp_data_holder2,1));
        full_time_choice_300 = cell(1,size(tmp_data_holder_pre_choice_300,1));
        full_time_choice_500 = cell(1,size(tmp_data_holder_pre_choice_500,1));
        full_time_pre_feedback = cell(1,size(tmp_data_holder_pre_feedback,1));
        full_time_first_fix = cell(1,size(tmp_data_holder_first_fix,1));
        full_time_finish_fix = cell(1,size(tmp_data_holder_finish_fix,1));
        
        cell_img_1_lfp = cell(1,size(tmp_data_holder1,1));
        cell_img_2_lfp = cell(1,size(tmp_data_holder2,1));
        cell_img_pre_choice_300 = cell(1,size(tmp_data_holder_pre_choice_300,1));
        cell_img_pre_choice_500 = cell(1,size(tmp_data_holder_pre_choice_500,1));
        cell_img_pre_feedback = cell(1,size(tmp_data_holder_pre_feedback,1));
        cell_img_first_fix = cell(1,size(tmp_data_holder_first_fix,1));
        cell_img_finish_fix = cell(1,size(tmp_data_holder_finish_fix,1));
        
        for irow = 1:size(tmp_data_holder1,1)
            cell_img_1_lfp{1,irow} = squeeze(tmp_data_holder1(irow,:,:));
            cell_img_2_lfp{1,irow} = squeeze(tmp_data_holder2(irow,:,:));
            cell_img_pre_choice_300{1,irow} = squeeze(tmp_data_holder_pre_choice_300(irow,:,:));
            cell_img_pre_choice_500{1,irow} = squeeze(tmp_data_holder_pre_choice_500(irow,:,:));
            cell_img_pre_feedback{1,irow} = squeeze(tmp_data_holder_pre_feedback(irow,:,:));
            cell_img_first_fix{1,irow} = squeeze(tmp_data_holder_first_fix(irow,:,:));
            cell_img_finish_fix{1,irow} = squeeze(tmp_data_holder_finish_fix(irow,:,:));
            
            full_time_1{1,irow} = time;
            full_time_2{1,irow} = time;
            full_time_choice_300{1,irow} = time;
            full_time_choice_500{1,irow} = time_500;
            full_time_pre_feedback{1,irow} = time_1000;
            full_time_first_fix{1,irow} = time_150;
            full_time_finish_fix{1,irow} = time_500;
        end
        
        tmp_data1.trial = cell_img_1_lfp;
        tmp_data1.time = full_time_1;
        tmp_data2.trial = cell_img_2_lfp;
        tmp_data2.time = full_time_2;
        tmp_data_choice_300.trial = cell_img_pre_choice_300;
        tmp_data_choice_300.time = full_time_choice_300;
        tmp_data_choice_500.trial = cell_img_pre_choice_500;
        tmp_data_choice_500.time = full_time_choice_500;
        tmp_data_pre_feedback.trial = cell_img_pre_feedback;
        tmp_data_pre_feedback.time = full_time_pre_feedback;
        tmp_data_first_fix.trial = cell_img_first_fix;
        tmp_data_first_fix.time = full_time_first_fix;
        tmp_data_finish_fix.trial = cell_img_finish_fix;
        tmp_data_finish_fix.time = full_time_finish_fix;

        % do the preprossesing now

%         p_cfg = [];
%         p_cfg.demean = 'yes';
%         p_cfg.detrend = 'yes';
%         tmp_data1 = ft_preprocessing(p_cfg, tmp_data1);
%         tmp_data2 = ft_preprocessing(p_cfg, tmp_data2);

        eval(sprintf('%s_1 = tmp_data1;', var_name));
        eval(sprintf('%s_2 = tmp_data2;', var_name));
        eval(sprintf('%s_pre_choice_300 = tmp_data_choice_300;', var_name));
        eval(sprintf('%s_pre_choice_500 = tmp_data_choice_500;', var_name));
        eval(sprintf('%s_pre_feedback = tmp_data_pre_feedback;', var_name));
        eval(sprintf('%s_first_fix = tmp_data_first_fix;', var_name));
        eval(sprintf('%s_finish_fix = tmp_data_finish_fix;', var_name));
        
        clear tmp_data1 tmp_data2 tmp_data_choice_300 tmp_data_choice_500 tmp_data_pre_feedback tmp_data_first_fix tmp_data_finish_fix;

        % save things off
        save(sprintf('%s.mat', var_name), sprintf('%s_1', var_name), ...
            sprintf('%s_2', var_name), sprintf('%s_pre_choice_300', var_name), ...
            sprintf('%s_pre_choice_500', var_name), ...
            sprintf('%s_pre_feedback', var_name), ...
            sprintf('%s_first_fix', var_name), ...
            sprintf('%s_finish_fix', var_name), ...
            sprintf('%s_ACC_channels', var_name),  ...
            sprintf('%s_DLPFC_channels', var_name),  ...
            sprintf('%s_OFC_channels', var_name));
        
        clear tmp_data1 tmp_data2 tmp_data_choice_300 tmp_data_choice_500 ...
            tmp_data_pre_feedback tmp_data_holder1 tmp_data_holder2 ...
            tmp_data_holder_pre_choice_300 tmp_data_holder_pre_choice_500 ...
            tmp_data_holder_pre_feedback;

    % if the 2 are the same
    else
%         disp(brain_region)
        if brain_region == 1
            channel_names(end+1,1) = {append(sstr,cstr)};
            ACC_channels(end+1) = {append(sstr,cstr)};
        elseif brain_region == 2
            channel_names(end+1,1) = {append(sstr,cstr)};
            DLPFC_channels(end+1) = {append(sstr,cstr)};
        elseif brain_region == 3
            channel_names(end+1,1) = {append(sstr,cstr)};
            OFC_channels(end+1) = {append(sstr,cstr)};
        end

        % only store if the brain region is 1, 2 or 3
        if brain_region == 1 || brain_region == 2 || brain_region == 3
%             disp('hit here')
            tmp_data_holder1(:,end+1,:) = reshape(trial_lfp_img1, [size(trial_lfp_img1,1), 1, size(trial_lfp_img1,2)]);
            tmp_data_holder2(:,end+1,:) = reshape(trial_lfp_img2, [size(trial_lfp_img2,1), 1, size(trial_lfp_img2,2)]);
            tmp_data_holder_pre_choice_300(:,end+1,:) = reshape(trial_lfp_pre_choice_300, [size(trial_lfp_pre_choice_300,1), 1, size(trial_lfp_pre_choice_300,2)]);
            tmp_data_holder_pre_choice_500(:,end+1,:) = reshape(trial_lfp_pre_choice_500, [size(trial_lfp_pre_choice_500,1), 1, size(trial_lfp_pre_choice_500,2)]);
            tmp_data_holder_pre_feedback(:,end+1,:) = reshape(trial_lfp_pre_feedback, [size(trial_lfp_pre_feedback,1), 1, size(trial_lfp_pre_feedback,2)]);
            tmp_data_holder_first_fix(:,end+1,:) = reshape(trial_lfp_first_fix, [size(trial_lfp_first_fix,1), 1, size(trial_lfp_first_fix,2)]);
            tmp_data_holder_finish_fix(:,end+1,:) = reshape(trial_lfp_finish_fix, [size(trial_lfp_finish_fix,1), 1, size(trial_lfp_finish_fix,2)]);
        end
    end


%         %% Testing things for making it in fieldtrip
%         data.label = {'15'};
%         data.fsample = 1000;
%         time = 0:(300-1);
%         time = time./data.fsample;
%         full_time = repmat(time,length(trial_lfp_img1),1);
%         data.time = num2cell(full_time,2)';
%         data.trial = num2cell(trial_lfp_img1, 2)';


    clear Codes picture_codes position_times pic_end_times end_time
end

sessions_M = unique(sessions_M);
sessions_F = unique(sessions_F);

save('monkey_sessions_M.mat', 'sessions_M');
save('monkey_sessions_F.mat', 'sessions_F');

% save(filename,variables)

%         %% Power analysis
%         tapsmofrq = [8 4 2]; %last with 12 instead 0f 9
%         foilim = {[20 140],[8 52],[2 15]}; %this is old{[20 140] [5 30] [2 10]};
%         pow_cfg = [];
%         pow_cfg.average = 'yes';
%         pow_cfg.method = 'mtmfft';
%         pow_cfg.keeptapers = 'yes';
%         pow_cfg.output = 'fourier';
%         pow_cfg.channel = S003_data_1.label;
% 
%         pow_cfg.keeptrials = 'yes';%'no';
%         pow_cfg.pad='maxperlen'; % Padding: not adding zeros
%         pow_cfg.flag = 0;
% 
%         pow_cfg.tapsmofrq = 1;
%         pow_cfg.foilim = [20 140];
%         pow_cfg.taper = 'hanning';
% 
%         tmp = ft_freqanalysis(pow_cfg,S003_data_1);
% 
%         tmp.freq;
% 
%         pow2_cfg = [];
%         pow2_cfg.psi = 'no'; % Phase slope Index
%         pow2_cfg.channel = S003_data_1.label;
%         pow2_cfg.jackknife = 'yes';
%         pow2_cfg.avgChOI = 'yes';
%         powtap8 = ft_freqdescriptives(pow2_cfg, tmp);
% %         loglog(powtap8.freq,(powtap8.powspctrm), ...
% %             'DisplayName', "Test")
%         loglog(powtap8.freq,(mean(powtap8.powspctrm(:,:))), ...
%             'DisplayName', "Test")
%         
%         
        