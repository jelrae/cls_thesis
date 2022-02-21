%% InfoGathering_Preprocess_Neurons
% preprocesses neurons, creates some useful variables for analysis, strips
% out analogue data (LFP, eye-tracking etc.) and stores results in
% neuronal_data/LTH_processed_data
clear

%% set up base directory where data/scripts are stored
[bd] = get_basedir;
addpath(genpath(fullfile(bd)));

datadir = fullfile(bd,'neuronal_data');
filenames=dir(fullfile(datadir,'*.mat'));

%% load in information about where each electrode (channel) is recorded from, in each animal

load('frank_area');frank_area=area_index; clear area_index;
load('miles_area');

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

for u=1:length(filenames) %u = 'unit'
    tic
    %% load in the data
    load(filenames(u).name);
    u
    filenames(u)
    
    disp('Epoch/pre-process');
    
    %% isolate information gathering trials where error code doesn't indicate
    %  timeout etc., and therefore where subject completed the trial
    inf_trials=find(BhvInfo.ConditionNumber==3|BhvInfo.ConditionNumber==4);
    trial_error=BhvInfo.TrialError; % 0 = 'Correct Response' (chose side with higher expected value),
    % 6 = 'Incorrect Response' (chose side with lower expected value)
    % 1 = No Response, 2 = Late Response, 3 = Break Fixation, 4 = No Fixation
    % 5 = Early Response, 7 = Lever Break, 8 = Ignored, 9 =Aborted
    inf_responded_trials=inf_trials(trial_error(inf_trials)==0|trial_error(inf_trials)==6);
    inf_cond=BhvInfo.ConditionNumber(inf_responded_trials);
    condition=BhvInfo.ConditionNumber(inf_responded_trials); %3 = option trial; 4 = attribute trial
    
    %% get information about session/channel from filename
    sb=filenames(u).name(1); %subject name, F=frank, M=miles
    % change name(2:5) to 2:4
    snum=str2num(filenames(u).name(2:4)); %session number, 1-49 for Miles, 1-36 for Frank
    chnum=str2num(filenames(u).name(strfind(filenames(u).name,'C')+1:strfind(filenames(u).name,'U')-2)); %channel number
    
    %find which row of miles_area/frank_area this corresponds to, and then
    %find out which brain region was recorded from
    if strcmp(sb,'F'),
        str_area=find(frank_area(:,1)==snum & frank_area(:,2)==chnum);
        brain_region=frank_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
    else
        str_area=find(miles_area(:,1)==snum&miles_area(:,2)==chnum);
        brain_region=miles_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
    end
    
    %% obtain the smoothed firing rate of the neuron across the entire session, and then calculate std/mean of this
    std_neurons(1:max(BhvInfo.OverallSpikeCodes))=0;
    std_neurons(BhvInfo.OverallSpikeCodes(2:end))=1000;
    
    smoothed_std_neurons=smooth_wa(std_neurons,200,'boxcar');
    
    sdev=nanstd(smoothed_std_neurons);
    men=nanmean(smoothed_std_neurons);
    
    %% now loop over trials, calculating key variables about stimuli presented and event onsets
    
    nInfoGatheringTrials = length(inf_responded_trials); %number of InfoGathering Trials
    
    for tr=1:nInfoGatheringTrials    %tr = 'trial'
        ct=inf_responded_trials(tr); %ct = current trial, note that this
        %     is different from tr because
        %     we are only interested in
        %     *completed info-gathering
        %     trials*
        
        %% Key variables for each trial:
        %Payoff/probability/expected value of each side:
        Left_pay(tr,1)=BhvInfo.Uservars{1,ct}.PayoffTracker(1);
        Left_prob(tr,1)=BhvInfo.Uservars{1,ct}.ProbabilityTracker(1);
        Right_pay(tr,1)=BhvInfo.Uservars{1,ct}.PayoffTracker(2);
        Right_prob(tr,1)=BhvInfo.Uservars{1,ct}.ProbabilityTracker(2);
        Left_EV(tr,1)=BhvInfo.Uservars{1,ct}.EVTracker(1);
        Right_EV(tr,1)=BhvInfo.Uservars{1,ct}.EVTracker(2);
        Pay_diff(tr)=Left_pay(tr)-Right_pay(tr);
        Prob_diff(tr)=Left_prob(tr)-Right_pay(tr);
        
        %Which option was eventually chosen:
        Chosen_target(tr,1)=BhvInfo.Uservars{1,ct}.ChosenTarget; %1 if left, 2 if right
        
        %Whether probability was on top/bottom row
        Prob_top(tr)=BhvInfo.Uservars{1,ct}.ProbTop; %1 if top, 2 if bottom
        
        %Code numbers/times - will use these to work out the sequence of
        %information sampling, amongst other things
        Codes{tr}(:,1)=BhvInfo.CodeNumbers{1,ct};
        Codes{tr}(:,2)=BhvInfo.CodeTimes{1,ct};
        
        %% Now calculate the timing of key events in the trial
        
        Fix_end(tr)=Codes{tr}(find(Codes{tr}(:,1)==12),2); %end of fixation period
        Response(tr)=Codes{tr}(find(Codes{tr}(:,1)==72),2); %time of joystick response
        
        %time of first sizeable joystick movement, Move_time (should precede 'Response'):
        Joy_time{tr}=BhvInfo.AnalogData{ct}.Joystick(Fix_end(tr):Response(tr),1); %analog data for X-position of joystick
        temp=min(find(abs(Joy_time{tr})>2)); %find first time that this exceeds 2
        if isempty(temp)
            Move_time(tr)=NaN;
        else
            Move_time(tr)=temp+Fix_end(tr);
        end
        
        %times when pictures were first fixated, position_times, and were removed, pic_end_times
        %(91=top left, 92 = bottom left, 93 = top right, 94 = bottom right):
        position_code{tr}=BhvInfo.CodeNumbers{ct}(find(BhvInfo.CodeNumbers{ct}>90&BhvInfo.CodeNumbers{ct}<95));
        position_times{tr}=BhvInfo.CodeTimes{ct}(find(BhvInfo.CodeNumbers{ct}>90&BhvInfo.CodeNumbers{ct}<95));
        pic_end_times{tr}(:,1)=BhvInfo.CodeNumbers{ct}(find(BhvInfo.CodeNumbers{ct}==90|BhvInfo.CodeNumbers{ct}==72));
        pic_end_times{tr}(:,2)=BhvInfo.CodeTimes{ct}(find(BhvInfo.CodeNumbers{ct}==90|BhvInfo.CodeNumbers{ct}==72));
        
        %length of time final picture was viewed for (note that this can be >300ms
        % if subject took long time to respond after cue 4)
        last_pic_view_time(tr)=pic_end_times{tr}(end,2)- position_times{tr}(end);
        
        % an alternative way of calculating this, using the time of first making
        % the joystick movement
        if ~isnan(Move_time(tr));
            new_last_time(tr)=Move_time(tr)-position_times{tr}(end);
        else
            new_last_time(tr)=300;
        end
        
        
        %% Now calculate information about the choice made, what influenced it
        
        %was this a 'brainer' trial (where probability told you to go one
        %way, and magnitude the opposite?)
        if sign(Pay_diff(tr))~=sign(Prob_diff(tr))&& ... %checks whether differences go in opposite directions
                sum(abs([sign(Prob_diff(tr)) sign(Pay_diff(tr))]))>1 %checks that neither difference is 0
            brainer(tr)=1; %brainer trial
        else
            brainer(tr)=0; %nobrainer trial
        end
        
        %was current trial a left choice? was current choice driven by
        %Prob_diff rather than Pay_diff?
        if Chosen_target(tr,1)==1,
            Left_choice(tr)=1;
            if Prob_diff(tr)==0&&Pay_diff(tr)==0
                P_choice(tr)=0;
            elseif  Prob_diff(tr)>=0&&Pay_diff(tr)<=0;
                P_choice(tr)=1;
            else
                P_choice(tr)=0;
            end
        else
            Left_choice(tr)=0;
            if Prob_diff(tr)==0&&Pay_diff(tr)==0;
                P_choice(tr)=0;
            elseif Prob_diff(tr)<=0&&Pay_diff(tr)>=0;
                P_choice(tr)=1;
            else
                P_choice(tr)=0;
            end
        end
        
        %% Now calculate information about information sampled during trial, spatial position, etc.
        
        %position of first and second spatial cues (1=top left, 2 = bottom
        %left, 3 = top right, 4 = bottom right)
        First_pos(tr,1)=BhvInfo.Uservars{1,ct}.Remove_index;
        Second_pos(tr,1)=BhvInfo.Uservars{1,ct}.Sind;
        
        %this is similar to position_code/position_times, but it is in Plexon's
        %timing, rather than MonkeyLogic's:
        picture_codes{tr}=BhvInfo.PlexonStrobes{ct}(find(BhvInfo.PlexonStrobes{ct}(:,1)>90&BhvInfo.PlexonStrobes{ct}(:,1)<95),:);
        picture_codes{tr}(:,2)=picture_codes{tr}(:,2)-BhvInfo.PlexonStrobes{ct}(1,2);
        
        %this information is similar to information cotained in
        %Pic_end_times, but it is in Plexon's timing, rather than
        %MonkeyLogic's:
        Pic_ends{tr}=BhvInfo.PlexonStrobes{ct}(find(BhvInfo.PlexonStrobes{ct}(:,1)==90),:);
        Pic_ends{tr}(:,2)=Pic_ends{tr}(:,2)-BhvInfo.PlexonStrobes{ct}(1,2);
        
        %was first side left (1) or right (0)?
        if picture_codes{tr}(1,1)==91||picture_codes{tr}(1,1)==92,
            First_side(tr)=1;
        else
            First_side(tr)=0;
        end
        
        %was second side left (1) or right (0)?
        if picture_codes{tr}(2,1)==91||picture_codes{tr}(2,1)==92,
            Second_side(tr)=1;
        else
            Second_side(tr)=0;
        end
        
        %was third side left (1) or right (0)?
        if length(picture_codes{tr}(:,1))>2,
            Third_pos(tr)=picture_codes{tr}(3,1);
            if picture_codes{tr}(3,1)==91||picture_codes{tr}(3,1)==92,
                Third_side(tr)=1;
            else
                Third_side(tr)=0;
            end
        else
            Third_side(tr)=NaN;
            Third_pos(tr)=NaN;
        end
        
        %was fourth side left (1) or right (0)?
        if length(picture_codes{tr}(:,1))>3,
            if picture_codes{tr}(4,1)==91||picture_codes{tr}(3,1)==92,
                Fourth_side(tr)=1;
            else
                Fourth_side(tr)=0;
            end
        else
            Fourth_side(tr)=NaN;
        end
        
        %how many pictures were revealed (but see below for correction by
        %how long final picture was viewed for)
        pictures_views(tr)=length(picture_codes{tr}(:,1));
        
        %set up 'spatial_values':
        % spatial_values{tr}(1,1) = top left value
        % spatial_values{tr}(1,2) = top right value
        % spatial_values{tr}(2,1) = bottom left value
        % spatial_values{tr}(2,2) = bottom right value
        if Prob_top(tr)==1,
            spatial_values{tr}=[Left_prob(tr), Right_prob(tr);
                Left_pay(tr), Right_pay(tr)];
        elseif Prob_top(tr)==2;
            spatial_values{tr}=[Left_pay(tr), Right_pay(tr);
                Left_prob(tr), Right_prob(tr)];
        end
        
        %% Calculate further timings of key events
        
        %timing of Cue 1,2,3,4 onsets:
        Cue1_ON(tr)=picture_codes{tr}(1,2);
        Cue2_ON(tr)=picture_codes{tr}(2,2);
        if length(picture_codes{tr}(:,1))>2,
            Cue3_ON(tr)=picture_codes{tr}(3,2);
        else
            Cue3_ON(tr)=NaN;
        end
        if length(picture_codes{tr}(:,1))>3,
            Cue4_ON(tr)=picture_codes{tr}(4,2);
        else
            Cue4_ON(tr)=NaN;
        end
        
        %timing of last cue, 2nd to last cue, 3rd to last cue, etc. stored in Last_Cue_Mat:
        rev_picture_codes{tr}=flipud(picture_codes{tr});
        Last_Cue_Mat(tr,1)=rev_picture_codes{tr}(1,2);
        Last_Cue_Mat(tr,2)=rev_picture_codes{tr}(2,2);
        if length(rev_picture_codes{tr}(:,1))>2,
            Last_Cue_Mat(tr,3)=rev_picture_codes{tr}(3,2);
        else
            Last_Cue_Mat(tr,3)=NaN;
        end
        if length(rev_picture_codes{tr}(:,1))>3,
            Last_Cue_Mat(tr,4)=rev_picture_codes{tr}(4,2);
        else
            Last_Cue_Mat(tr,4)=NaN;
        end
        
        Plex_codes{tr}=BhvInfo.PlexonStrobes{ct};  %all Plexon codes and their timings for this trial
        Realign_plex_code(tr)=Plex_codes{tr}(1,2); %time of first code, to which timings may be realigned - normally 0
        Plex_codes{tr}(:,2)=round((Plex_codes{tr}(:,2)-Realign_plex_code(tr))); %realign timings
        
        %timing of response, end of fixation, reward onset,
        Response_time(tr)=Plex_codes{tr}(find(Plex_codes{tr}(:,1)==72),2);
        FixationEND(tr)=Plex_codes{tr}(find(Plex_codes{tr}(:,1)==12),2);
        Reward_ON(tr)=Plex_codes{tr}(find(Plex_codes{tr}(:,1)==1),2);
        GO_cue(tr)=Pic_ends{tr}(1,2); %we no longer have explicit go cue, but red 'no-go' cue is removed here
        
        
        %% Now get the cell's firing on this trial
        SpikeTable{tr}=BhvInfo.SpikeCodes{ct}; %timing of action potentials in the current trial
        
        %calculate number of spikes in various events:
        %**N.B. Most analyses don't use these - use Cue1_Matrix etc.,
        %       created below. **
        %FSpikes is spikes during fixation epoch
        %C1Spikes, C2Spikes, C3Spikes, C4Spikes are spikes during Cues 1-4
        %LSpikes are spikes during final cue presentation
        
        FSpikes(tr)=(sum(SpikeTable{tr}>FixationEND(tr)-299&SpikeTable{tr}<=FixationEND(tr))/300)*1000;
        C1Spikes(tr)=(sum(SpikeTable{tr}>Cue1_ON(tr)&SpikeTable{tr}<=Cue1_ON(tr)+299)/300)*1000;
        C2Spikes(tr)=(sum(SpikeTable{tr}>Cue2_ON(tr)&SpikeTable{tr}<=Cue2_ON(tr)+299)/300)*1000;
        if ~isnan(Cue3_ON(tr)),
            C3Spikes(tr)=(sum(SpikeTable{tr}>Cue3_ON(tr)&SpikeTable{tr}<=Cue3_ON(tr)+299)/300)*1000;
        else
            C3Spikes(tr)=NaN;
        end
        if ~isnan(Cue4_ON(tr)),
            C4Spikes(tr)=(sum(SpikeTable{tr}>Cue4_ON(tr)&SpikeTable{tr}<=Cue4_ON(tr)+299)/300)*1000;
        else
            C4Spikes(tr)=NaN;
        end
        LSpikes(tr,1)=sum((SpikeTable{tr}>Last_Cue_Mat(tr,1)&SpikeTable{tr}<=Last_Cue_Mat(tr,1)+299)/300)*1000;
        LSpikes(tr,2)=sum((SpikeTable{tr}>Last_Cue_Mat(tr,2)&SpikeTable{tr}<=Last_Cue_Mat(tr,2)+299)/300)*1000;
        if ~isnan(Last_Cue_Mat(tr,3)),
            LSpikes(tr,3)=sum((SpikeTable{tr}>Last_Cue_Mat(tr,3)&SpikeTable{tr}<=Last_Cue_Mat(tr,3)+299)/300)*1000;
        else
            LSpikes(tr,3)=NaN;
        end
        if ~isnan(Last_Cue_Mat(tr,4)),
            LSpikes(tr,4)=sum((SpikeTable{tr}>Last_Cue_Mat(tr,4)&SpikeTable{tr}<=Last_Cue_Mat(tr,4)+299)/300)*1000;
        else
            LSpikes(tr,4)=NaN;
        end
        
    end %of loop over trials
    
    %% did subject choose side on which first,second,third,fourth cue was presented?
    First_side_chosen=First_side==Left_choice;
    Second_side_chosen=Second_side==Left_choice;
    Third_side_chosen=Third_side==Left_choice;
    Fourth_side_chosen=Fourth_side==Left_choice;
    
    %% now epoch all trials to make rasters
    pre_cue=199; %time pre-event (in ms)
    post_cue=600; %time post-event (in ms)
    
    %make rasters:
    Cue1_Matrix = neuron_matrix(pre_cue,Cue1_ON,post_cue,SpikeTable);
    Cue2_Matrix = neuron_matrix(pre_cue,Cue2_ON,post_cue,SpikeTable);
    Cue3_Matrix = neuron_matrix(pre_cue,Cue3_ON,post_cue,SpikeTable);
    Response_Matrix = neuron_matrix(post_cue,Response_time,pre_cue,SpikeTable);
    % Cue4_Matrix = neuron_matrix(pre_cue,Cue4_ON,post_cue,SpikeTable);
    FixMatrix = neuron_matrix(250,FixationEND-250,250,SpikeTable);
    
    %store 'original' rasters, that simply contain a 1 for each spike, 0
    %elsewhere, as Cue1_MatrixRaw etc.
    Cue1_MatrixRaw = Cue1_Matrix/1000; %divide by 1000, as neuron_matrix places 1000 in each spike location
    Cue2_MatrixRaw = Cue2_Matrix/1000;
    Cue3_MatrixRaw = Cue3_Matrix/1000;
    %Cue4_MatrixRaw = Cue4_Matrix/1000;
    Response_MatrixRaw = Response_Matrix/1000;
    FixMatrixRaw = FixMatrix/1000;
    
    %normalise rasters by baseline average firing rate/standard deviation
    Cue1_Matrix=(Cue1_Matrix-men)/sdev;
    Cue2_Matrix=(Cue2_Matrix-men)/sdev;
    Cue3_Matrix=(Cue3_Matrix-men)/sdev;
    Response_Matrix=(Response_Matrix-men)/sdev;
    % Cue4_Matrix=(Cue4_Matrix-men)/sdev;
    
    %% now do some further processing on cues viewed by subject on each trial
    
    % Lreg will be a matrix of nTrials*4
    % If a cue appeared on the left-hand side on Cue 1,2,3 or 4, it will
    % contain the *picture rank* of this picture - i.e. lowest
    % probability/magnitude picture = 1, highest probability/magnitude
    % picture = 5
    % If this cue was presented on the right-hand cue, or if cue was not
    % sampled (can only apply to cues 3 and 4), LReg will be valued 0.
    Lreg=zeros(nInfoGatheringTrials,4);
    
    % Rreg will be similar, but when cues were presented on RHS of screen
    % rather than LHS
    Rreg=zeros(nInfoGatheringTrials,4);
    
    %lastLreg/lastRreg will be similar, but column 1 is final cue presented, column 2 is
    %penultimate cue presented, column 3 is antepenultimate (if >2 cues
    %sampled), column 4 is preantepenultimate cue (if >3 cues sampled)
    lastLreg=zeros(nInfoGatheringTrials,4);
    lastRreg=zeros(nInfoGatheringTrials,4);
    
    for t=1:nInfoGatheringTrials %loop over Info Gathering trials
        
        %% calculate number of pictures viewed
        num_pics_viewed(t)=length(position_code{t});
        % check whether final cue was viewed sufficently long to be
        % considered 'viewed' (new_correct is in ms, defined above)
        if (new_last_time(t)<=new_correct) & new_correct~=0 & num_pics_viewed(t)>2,
            num_pics_viewed(t)=num_pics_viewed(t)-1;
        end
        
        last_pos_code{t}=fliplr(position_code{t}); %cue locations sorted from end to beginning
        
        temp_pos_code{t}=position_code{t}(1:num_pics_viewed(t)); %only consider pictures viewed
        
        %% calculate the probability/magnitude seen on each side (set to 0 if not seen)
        [LP_seen(t),LM_seen(t),RP_seen(t),RM_seen(t)]=calculate_atts_seen...
            (temp_pos_code{t},Prob_top(t),Left_prob(t),Left_pay(t),...
            Right_prob(t),Right_pay(t));
        
        %% now fill up Lreg, Rreg, lastLreg, lastRreg for this trial (see above for description)
        for y=1:length(position_code{t})
            if position_code{t}(y)==91
                Lreg(t,y)=calc_rank(spatial_values{t}(1,1));
            elseif position_code{t}(y)==92
                Lreg(t,y)=calc_rank(spatial_values{t}(2,1));
            elseif position_code{t}(y)==93
                Rreg(t,y)=calc_rank(spatial_values{t}(1,2));
            elseif position_code{t}(y)==94
                Rreg(t,y)=calc_rank(spatial_values{t}(2,2));
            end
            
            if last_pos_code{t}(y)==91
                lastLreg(t,y)=calc_rank(spatial_values{t}(1,1));
            elseif last_pos_code{t}(y)==92
                lastLreg(t,y)=calc_rank(spatial_values{t}(2,1));
            elseif last_pos_code{t}(y)==93
                lastRreg(t,y)=calc_rank(spatial_values{t}(1,2));
            elseif last_pos_code{t}(y)==94
                lastRreg(t,y)=calc_rank(spatial_values{t}(2,2));
            end
        end %loop over positions seen on this trial
        
        %% calculate value of first picture
        if position_code{t}(1)==91
            first_picture(t)=spatial_values{t}(1,1);
        elseif position_code{t}(1)==92
            first_picture(t)=spatial_values{t}(2,1);
        elseif position_code{t}(1)==93
            first_picture(t)=spatial_values{t}(1,2);
        elseif position_code{t}(1)==94
            first_picture(t)=spatial_values{t}(2,2);
        end
        
        %because probabilities are 0.1:0.2:0.9; magnitudes are
        %0.15:0.2:0.95...
        if rem(first_picture(t),0.1)==0 %if this is true
            P_first(t)=1; %then proability was first cue on this trial
        else
            P_first(t)=0; %magnitude was first cue
        end
        
        %% now do same as above for second cue:
        if position_code{t}(2)==91
            second_picture(t)=spatial_values{t}(1,1);
        elseif position_code{t}(2)==92
            second_picture(t)=spatial_values{t}(2,1);
        elseif position_code{t}(2)==93
            second_picture(t)=spatial_values{t}(1,2);
        elseif position_code{t}(2)==94
            second_picture(t)=spatial_values{t}(2,2);
        end
        if rem(second_picture(t),0.1)==0
            P_second(t)=1;
        else
            P_second(t)=0;
        end
        
        %% now do same as above for third cue:
        if num_pics_viewed(t)>2
            if position_code{t}(3)==91
                third_picture(t)=spatial_values{t}(1,1);
            elseif position_code{t}(3)==92
                third_picture(t)=spatial_values{t}(2,1);
            elseif position_code{t}(3)==93
                third_picture(t)=spatial_values{t}(1,2);
            elseif position_code{t}(3)==94
                third_picture(t)=spatial_values{t}(2,2);
            end
            
            if rem(third_picture(t),0.1)==0
                P_third(t)=1;
            else
                P_third(t)=0;
            end
        end
        if num_pics_viewed(t)>2
            if rem(spatial_values{t}(position_code{t}(3)-90),0.1)~=0
                thirdattP(t)=0;
            else
                thirdattP(t)=1;
            end
        else
            thirdattP(t)=NaN;
        end
        
        %% now calcalate whether third saccade was vertical on attribute trials:
        if num_pics_viewed(t)>2&&inf_cond(t)==4 %attribute trial with >2 cues seen
            if position_code{t}(2)==91||position_code{t}(2)==92 %second cue on left option
                if position_code{t}(3)==92||position_code{t}(3)==91 %third cue also on left option
                    vert(t)=1; %hence vertical saccade
                else
                    vert(t)=0; %diagonal saccade back to first option
                end
            elseif position_code{t}(2)==93||position_code{t}(2)==94 %second cue right
                if position_code{t}(3)==94||position_code{t}(3)==93 %third cue also right
                    vert(t)=1;
                else
                    vert(t)=0;
                end
            else
                vert(t)=NaN;
            end
        else
            vert(t)=NaN; %option trial, or attribute trial with only 2 cues
        end
        
        %% now calculate whether third saccade was horizontal on option trials:
        if num_pics_viewed(t)>2&&inf_cond(t)==3 %option trial with >2 cues seen
            if position_code{t}(2)==91||position_code{t}(2)==93 %second cue top, etc.
                if position_code{t}(3)==93||position_code{t}(3)==91
                    horz(t)=1;
                else
                    horz(t)=0;
                end
            elseif position_code{t}(2)==92||position_code{t}(2)==94
                if position_code{t}(3)==94||position_code{t}(3)==92
                    horz(t)=1;
                else
                    horz(t)=0;
                end
            else
                horz(t)=NaN;
            end
        else
            horz(t)=NaN;
        end
        
    end %of loop over InfoGathering trials

    %% a few more variables to create before we save...
    
    %reg_mat - this is all cues in sequence, and mid-ranked
    %picture (3) is set to have average value of 0 (note that for columns 3
    %and 4, this means that non-viewed pictures by default have average
    %value, as they were already set for zero
    reg_mat=Lreg+Rreg;
    reg_mat(reg_mat~=0)=reg_mat(reg_mat~=0)-3;
    
    %thirdattP is adjusted to be -1 on option trials where third attribute is
    %magniutde, and zero on attribute trials/trials where third cue was not seen - 
    %this means it will be useful as a regressor for *which attribute* the subject
    %decided to look at on the third saccade of option trials
    thirdattP(thirdattP==0)=-1;
    thirdattP(isnan(thirdattP))=0;
    thirdattP(inf_cond==4)=0; %attribute trials
    
    %similarly for thirdtop,  which is set to be 1 on option trials
    %where the third saccade was to the top cue, and -1 on option
    %trials where it was to the bottom cue
    Third_pos=str2num(num2str(Third_pos));
    thirdtop=ones(length(Third_pos),1);
    thirdtop(Third_pos==92|Third_pos==94)=-1;
    thirdtop(inf_cond==4)=0; %attribute trials
    
    %similarly for thirdsacdir, which is set to be 1 on attribute trials
    %where the third saccade was to the left cue, and -1 on attribute
    %trials where it was to the right cue
    thirdsacdir=(Third_side==1)-(Third_side==0);
    thirdsacdir(inf_cond==3)=0; %option trials


    %similarly, attendsac1 and attendsac2 tell us on attribute trials
    %whether the third saccade was directed back towards the first viewed
    %option (1), or towards the second viewed option (-1)
    attendsac1=(First_side==Third_side)-(First_side~=Third_side);
    attendsac1(inf_cond==3)=0; %option trials
    attendsac2=(Second_side==Third_side)-(Second_side~=Third_side);
    attendsac2(inf_cond==3)=0; %option trials

    %on attribute trials, optsacc tries to ask: do I saccade to the 
    %side with current higher value (1) or not (-1)?
    optdir4=(reg_mat(:,2)>reg_mat(:,1))&inf_cond'==4;
    optsacc=double((vert'==optdir4)&inf_cond'==4);
    optsacc((vert'~=optdir4)&inf_cond'==4)=-1;
        
    npv=num_pics_viewed;

    
    %% these variables are used in subsequent analyses as smoothing parameters
    window_size=200; %ms
    slide_width=10; %ms
    epoch_duration=length(Cue3_Matrix(1,:)); %ms
    pre_cue=pre_cue; %ms
    num_perms=10; %ms
    
    variables(1)=window_size;
    variables(2)=slide_width;
    variables(3)=epoch_duration;
    variables(4)=pre_cue;
    variables(5)=num_perms;
    
    %% remove analog data from BhvInfo, as this is what takes up the space
    BhvInfo = rmfield(BhvInfo,'AnalogData');
    BhvInfo = rmfield(BhvInfo,'PlexonEyeData1');
    BhvInfo = rmfield(BhvInfo,'PlexonEyeData2');
    BhvInfo = rmfield(BhvInfo,'PlexonEyeData3');
    BhvInfo = rmfield(BhvInfo,'PlexonFieldPotentials');
    if isfield(BhvInfo,'RewardRecord')
        BhvInfo = rmfield(BhvInfo,'RewardRecord');
    end
    
%     %% define filename
%     fname = fullfile(datadir,'LTH_processed_data',sprintf('Unit%04.0f.mat',u));
% 
%     %% save file
%     save(fname,'-regexp',['^(?!(CM1|CM2|CM3|RM|Joy_time|smoothed_std_neurons|' ...
%         'std_neurons|frank_area|miles_area|outM|datadir|' ...
%         'filenames|bd|outF|u|fname)$).']);
    
    %% clear variables for next iteration of loop
    clearvars -except frank_area miles_area  outM ...
        datadir filenames bd outF u new_correct;
    toc
end %of loop over units