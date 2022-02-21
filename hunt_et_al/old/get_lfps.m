%% Extracting LFP data

% This script is based off the scripts provided by Hunt et al. 

%% Add paths and load relevant information
addpath('D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/get_basedir.m');
addpath('D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etall/Data/unprocessed/meta_data');
load('D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/meta_data/all_units_info.mat', 'all_units_info.mat');
load('D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/meta_data/frank_area.mat', 'frank_area.mat');
load('D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/meta_data/miles_area.mat', 'miles_area.mat');

%set up base directory where data/scripts are stored
[bd] = get_basedir; 
addpath(genpath(fullfile(bd)));

datadir = bd;
filenames=dir(fullfile(datadir,'*.mat'));

% load in information about where each electrode (channel) is recorded from, in each animal

load('frank_area');frank_area=area_index; clear area_index;
load('miles_area');

%% Get LFP data + information
%TODO: this needs to loop over every file in directory and extract the
%information - i.e. name needs to change with each loop

% Get the LFP data
for u=1:length(filenames) %u = 'unit'
    tic
    %% load in the data
    load(filenames(u).name);
    u;
    
     disp('Epoch/pre-process');
    
    % isolate information gathering trials where error code doesn't indicate
    %  timeout etc., and therefore where subject completed the trial
    inf_trials=find(BhvInfo.ConditionNumber==3|BhvInfo.ConditionNumber==4);
    trial_error=BhvInfo.TrialError; % 0 = 'Correct Response' (chose side with higher expected value),
    % 6 = 'Incorrect Response' (chose side with lower expected value)
    % 1 = No Response, 2 = Late Response, 3 = Break Fixation, 4 = No Fixation
    % 5 = Early Response, 7 = Lever Break, 8 = Ignored, 9 =Aborted
    inf_responded_trials=inf_trials(trial_error(inf_trials)==0|trial_error(inf_trials)==6);
    inf_cond=BhvInfo.ConditionNumber(inf_responded_trials);
    condition=BhvInfo.ConditionNumber(inf_responded_trials); %3 = option trial; 4 = attribute trial
    
    % get information about session/channel from filename
    sb=filenames(u).name(1); %subject name, F=frank, M=miles

    snum=str2num(filenames(u).name(2:4)); %session number, 1-49 for Miles, 1-36 for Frank
    chnum=str2num(filenames(u).name(strfind(filenames(u).name,'C')+1:strfind(filenames(u).name,'U')-2)); %channel number
    
    %find which row of miles_area/frank_area this corresponds to, and then
    %find out which brain region was recorded from
    if strcmp(sb,'F'),
        str_area=find(frank_area(:,1)==snum&frank_area(:,2)==chnum);
        brain_region=frank_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
    else
        continue
%         str_area=find(miles_area(:,1)==snum&miles_area(:,2)==chnum);
%         brain_region=miles_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
    end
    
    InfoLFP = getfield(BhvInfo,'PlexonFieldPotentials');
    TimingStrobe = getfield(BhvInfo,'PlexonStrobes');
    nInfoGatheringTrials = length(inf_responded_trials); %number of InfoGathering Trials
    
    % now loop over trials, calculating key variables about stimuli presented and event onsets
    for tr=1:nInfoGatheringTrials    %tr = 'trial'
        ct=inf_responded_trials(tr);
        Plex_codes{tr}=BhvInfo.PlexonStrobes{ct};  %all Plexon codes and their timings for this trial
        Realign_plex_code(tr)=Plex_codes{tr}(1,2); %time of first code, to which timings may be realigned - normally 0
        Plex_codes{tr}(:,2)=round((Plex_codes{tr}(:,2)-Realign_plex_code(tr))); %realign timings

        %timing of response, end of fixation, reward onset,
        Response_time(tr)=Plex_codes{tr}(find(Plex_codes{tr}(:,1)==72),2);
        FixationEND(tr)=Plex_codes{tr}(find(Plex_codes{tr}(:,1)==12),2);
        Reward_ON(tr)=Plex_codes{tr}(find(Plex_codes{tr}(:,1)==1),2);

    end
%     Cue1_Matrix = neuron_matrix(pre_cue,Cue1_ON,post_cue);
    
    %% Save files and clear workdspace
    fname = fullfile(datadir,sprintf('LTH_processed_data_Unit%04.0f.mat',u));
    % Save all variables except those listed
    save(fname,'-regexp',['^(?!(CM1|CM2|CM3|RM|Joy_time|smoothed_std_neurons|' ...
        'std_neurons|frank_area|miles_area|outM|datadir|' ...
        'filenames|bd|outF|u|fname)$).']); 
    
    clearvars -except frank_area miles_area  outM ...
        datadir filenames bd outF u new_correct;
%     toc
end
 


% end
