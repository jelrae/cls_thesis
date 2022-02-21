function out = InfoGathering_Behaviour(zInfoGatheringInfo)

% InfoGathering_Behaviour(zInfoGatheringInfo)
%
% Behavioural analyses of Informtation Gathering task
%
% zInfoGatheringInfo is a choice behaviour structure,
% loaded in from Miles/Frank InfoGathering All Data

%% pull out session numbers and trial number etc.

for s=1:length(zInfoGatheringInfo)
    if s==1
        MTrialNum=zInfoGatheringInfo{s}.TrialNumber;
        Sess_l(s)=length(MTrialNum);
        Sess_ind(1:Sess_l(s)) = 1;
    else
        MTrialNum(end+1:(end+length(zInfoGatheringInfo{s}.TrialNumber)))=zInfoGatheringInfo{s}.TrialNumber;
        Sess_l(s)=length(MTrialNum);
        Sess_ind(Sess_l(s-1)+1:Sess_l(s)) = s;
    end
end

%% grab all uservars, analogdata and codes

for s=1:length(zInfoGatheringInfo)
    if s==1 %first session
        MCoditionNum=zInfoGatheringInfo{s}.ConditionNumber;
        MTrialError=zInfoGatheringInfo{s}.TrialError;
        MNumCodes=zInfoGatheringInfo{s}.NumCodes;
        for t=1:length(zInfoGatheringInfo{s}.ConditionNumber)
            if t==1 %first trial
                MCodeNumbers{1}=zInfoGatheringInfo{s}.CodeNumbers(t);
                MCodeTimes{1}=zInfoGatheringInfo{s}.CodeTimes(t);
                MAnalog{1}=zInfoGatheringInfo{s}.AnalogData(t);
                MUservars{1}=zInfoGatheringInfo{s}.Uservars(t);
            else
                MCodeNumbers{t}=zInfoGatheringInfo{s}.CodeNumbers(t);
                MCodeTimes{t}=zInfoGatheringInfo{s}.CodeTimes(t);
                MAnalog{t}=zInfoGatheringInfo{s}.AnalogData(t);
                MUservars{t}=zInfoGatheringInfo{s}.Uservars(t);
            end
        end
    else
        for t=1:length(zInfoGatheringInfo{s}.ConditionNumber),
            MCoditionNum(end+1)=zInfoGatheringInfo{s}.ConditionNumber(t);
            MTrialError(end+1)=zInfoGatheringInfo{s}.TrialError(t);
            MNumCodes(end+1)=zInfoGatheringInfo{s}.NumCodes(t);
            MCodeNumbers{end+1}=zInfoGatheringInfo{s}.CodeNumbers(t);
            MCodeTimes{end+1}=zInfoGatheringInfo{s}.CodeTimes(t);
            MAnalog{end+1}=zInfoGatheringInfo{s}.AnalogData(t);
            MUservars{end+1}=zInfoGatheringInfo{s}.Uservars(t);
        end
    end
end

%% extract information about values, probabilities, magnitudes, etc.
% also, preprocesses joystick/saccade info

info_trials=find((MTrialError==0|MTrialError==6));
new_correct=0; %threshold (in ms) for excluding 'unseen' stimuli when they are still fixated

for v=1:length(info_trials) %loop across trials
    Left_pay(v,1)=MUservars{v}{1}{1}.PayoffTracker(1); %0.15 to 0.95
    Left_prob(v,1)=MUservars{v}{1}{1}.ProbabilityTracker(1); %0.1 to 0.9
    Right_pay(v,1)=MUservars{v}{1}{1}.PayoffTracker(2); %0.15 to 0.95
    Right_prob(v,1)=MUservars{v}{1}{1}.ProbabilityTracker(2); %0.1 to 0.9
    Left_EV(v)=MUservars{v}{1}{1}.EVTracker(1);
    Right_EV(v)=MUservars{v}{1}{1}.EVTracker(2);
    chosen_target(v)=MUservars{v}{1}{1}.ChosenTarget; % 1 = left chosen, 2 = right chosen
    Prob_top(v)=MUservars{v}{1}{1}.ProbTop; % defines whether probability is on top or bottom - 1 means probability on top, 2 means bottom
    Prob_diff(v)=Left_prob(v)-Right_prob(v);
    Pay_diff(v)=Left_pay(v)-Right_pay(v);
    Rewarded(v)=MUservars{v}{1}{1}.Rewarded(1);
    
    % timestamps of key events
    Codes{v}(:,1)=MCodeNumbers{v};
    Codes{v}(:,2)=MCodeTimes{v};
    Fix_start(v)=Codes{v}{2}(find(Codes{v}{1}==11)); %timestamp of first entering fixation point at beginning of trial
    Fix_end(v)=Codes{v}{2}(find(Codes{v}{1}==12)); %timestamp of when fixation finishes
    Response(v)=Codes{v}{2}(find(Codes{v}{1}==72));     %timestamp of when joystick enters target area
    
    % calculate when joystick movement exceeded a threshold:
    
    Fix_Joy(v)=mean(MAnalog{1,v}{1}.Joystick(Fix_end(v)-100:Fix_end(v))); %Joystick during fixation period
    Joy_time{v}=smooth_wa(MAnalog{1,v}{1}.Joystick(Fix_end(v):Response(v),1),5,'boxcar');
    if abs(Fix_Joy(v))<0.8 %calibrated for Miles
        temp=find(abs(Joy_time{v})>1);
    else %calibrated for Frank
        temp=find(abs(Joy_time{v})>2.9);
    end
    
    Move_time(v)=temp(1)+Fix_end(v); %Move_time is when the joystick movement is considered to be 'committed':
    %any pictures viewed beyond this time
    %can be disregarded (see below)
    
    position_code{v}=Codes{v}{1}(find(Codes{v}{1}>90&Codes{v}{1}<95)); %code numbers (91-94) for four spatial locations, when they first fixate it
    position_times{v}=Codes{v}{2}(find(Codes{v}{1}>90&Codes{v}{1}<95)); %times corresponding to codes
    pic_end_times{v}(:,1)=Codes{v}{1}(find(Codes{v}{1}==90|Codes{v}{1}==72)); %when stimulus disappears - 90 when fixated for 300ms, 72 if choice has been made
    pic_end_times{v}(:,2)=Codes{v}{2}(find(Codes{v}{1}==90|Codes{v}{1}==72));
    
    %last_pic_view_time(v)=pic_end_times{v}(end,2)- position_times{v}(end);
    
    new_last_time(v)=Move_time(v)-position_times{v}(end); %this replaces last_pic_view_time,
    %better metric of how many pictures
    %seen before joystick movement
    
    
    if Prob_top(v)==1 % probability on top row, magnitude on bottom
        spatial_values{v}=[Left_prob(v) Right_prob(v); Left_pay(v) Right_pay(v)];
    elseif Prob_top(v)==2 % magnitude on top row, probability on bottom
        spatial_values{v}=[Left_pay(v) Right_pay(v); Left_prob(v) Right_prob(v)];
    end
    
    num_pics_viewed(v)=length(position_code{v}); %total number of pictures revealed IGNORING time of joystick movement
    
    Eye_data{v}=MAnalog{1,v}{1}.EyeSignal; %first column is X-position, second column is Y-position
    %     freq_cut=35;
    %     sampling_rate=1000;
    %     butterworth_order=2;
    %     frw=freq_cut/(sampling_rate/2);
    %     [Bco Aco]=butter(butterworth_order,frw);
    %     Eye_data{v}(:,1)=filtfilt(Bco, Aco,Eye_data{v}(:,1));
    %     Eye_data{v}(:,2)=filtfilt(Bco, Aco,Eye_data{v}(:,2));
end

Left_choice=chosen_target==1;

%Lreg and Rreg will contain value of stimulus if subject is fixating that
%picture, and 0 if fixating the other side
Lreg(1:length(Left_EV),1:4)=zeros;
Rreg(1:length(Left_EV),1:4)=zeros;
lastLreg(1:length(Left_EV),1:4)=zeros;
lastRreg(1:length(Left_EV),1:4)=zeros;

lastLreg250(1:length(Left_EV),1:4)=zeros;
lastRreg250(1:length(Left_EV),1:4)=zeros;

Third_Cue_Position = nan(1,length(Left_EV));

%%
for t=1:length(Left_EV) %loop again across trials
    %position_code{t}=error_correct(position_code{t}); %LH commented out
    
    last_pos_code{t}=flipud(position_code{t});
    if (new_last_time(t)<=new_correct) & new_correct~=0 & num_pics_viewed(t)>2,
        num_pics_viewed(t)=num_pics_viewed(t)-1;    
    end
    
    %% calculate number of pictures viewed based upon time of joystick movement, and length of final fixation:
    
    %joystick movement before last picture was fixated?
    if new_last_time(t)<0 & num_pics_viewed(t)>2
        num_pics_viewedprezero(t)=num_pics_viewed(t)-1;
    else
        num_pics_viewedprezero(t)=num_pics_viewed(t);
    end
    
    %joystick movement before last picture was fixated for 100ms?
    if new_last_time(t)<=100 & num_pics_viewed(t)>2,
        num_pics_viewed100(t)=num_pics_viewed(t)-1;
    else
        num_pics_viewed100(t)=num_pics_viewed(t);
    end
    
    %joystick movement before last picture was fixated for 200ms?
    if new_last_time(t)<=200 & num_pics_viewed(t)>2,
        num_pics_viewed200(t)=num_pics_viewed(t)-1;
    else
        num_pics_viewed200(t)=num_pics_viewed(t);
    end
    
    %joystick movement before last picture was fixated for 250ms?
    if new_last_time(t)<=250 & num_pics_viewed(t)>2,
        num_pics_viewed250(t)=num_pics_viewed(t)-1;
    else
        num_pics_viewed250(t)=num_pics_viewed(t);
    end
    
    %joystick movement before last picture was fixated for 250ms?
    if new_last_time(t)<=300 & num_pics_viewed(t)>2,
        num_pics_viewed300(t)=num_pics_viewed(t)-1;
    else
        num_pics_viewed300(t)=num_pics_viewed(t);
    end
    
    temp_pos_code{t}=position_code{t}(1:num_pics_viewed(t));
    temp_pos_codeprezero{t}=position_code{t}(1:num_pics_viewedprezero(t));
    temp_pos_code100{t}=position_code{t}(1:num_pics_viewed100(t));
    temp_pos_code200{t}=position_code{t}(1:num_pics_viewed200(t));
    temp_pos_code250{t}=position_code{t}(1:num_pics_viewed250(t));
    temp_pos_code300{t}=position_code{t}(1:num_pics_viewed300(t));
    
    [LP_seen(t),LM_seen(t),RP_seen(t),RM_seen(t)]=calculate_atts_seen(temp_pos_code{t},Prob_top(t),Left_prob(t),Left_pay(t),Right_prob(t),Right_pay(t));
    [LP_seenprezero(t),LM_seenprezero(t),RP_seenprezero(t),RM_seenprezero(t)]=calculate_atts_seen(temp_pos_codeprezero{t},Prob_top(t),Left_prob(t),Left_pay(t),Right_prob(t),Right_pay(t));
    [LP_seen100(t),LM_seen100(t),RP_seen100(t),RM_seen100(t)]=calculate_atts_seen(temp_pos_code100{t},Prob_top(t),Left_prob(t),Left_pay(t),Right_prob(t),Right_pay(t));
    [LP_seen200(t),LM_seen200(t),RP_seen200(t),RM_seen200(t)]=calculate_atts_seen(temp_pos_code200{t},Prob_top(t),Left_prob(t),Left_pay(t),Right_prob(t),Right_pay(t));
    [LP_seen250(t),LM_seen250(t),RP_seen250(t),RM_seen250(t)]=calculate_atts_seen(temp_pos_code250{t},Prob_top(t),Left_prob(t),Left_pay(t),Right_prob(t),Right_pay(t));
    [LP_seen300(t),LM_seen300(t),RP_seen300(t),RM_seen300(t)]=calculate_atts_seen(temp_pos_code300{t},Prob_top(t),Left_prob(t),Left_pay(t),Right_prob(t),Right_pay(t));
    
    for y=1:length(temp_pos_code{t})  % loop through pictures seen
        if position_code{t}(y)==91
            Lreg(t,y)=calc_rank(spatial_values{t}(1,1)); %calcrank transforms 0.1 to 0.9 into 1 to 5, and 0.15 to 0.95 into 1 to 5
        elseif position_code{t}(y)==92
            Lreg(t,y)=calc_rank(spatial_values{t}(2,1));
        elseif position_code{t}(y)==93
            Rreg(t,y)=calc_rank(spatial_values{t}(1,2));
        elseif position_code{t}(y)==94
            Rreg(t,y)=calc_rank(spatial_values{t}(2,2));
        end
        
        if last_pos_code{t}(y)==91 %for coding backwards from last fixated picture
            lastLreg(t,y)=calc_rank(spatial_values{t}(1,1));
        elseif last_pos_code{t}(y)==92
            lastLreg(t,y)=calc_rank(spatial_values{t}(2,1));
        elseif last_pos_code{t}(y)==93
            lastRreg(t,y)=calc_rank(spatial_values{t}(1,2));
        elseif last_pos_code{t}(y)==94
            lastRreg(t,y)=calc_rank(spatial_values{t}(2,2));
        end
        
        if y==3
            Third_Cue_Position(t)=temp_pos_code{t}(y)-90;
        end
    end
    
    if Lreg(t,2)>Rreg(t,2)
        second_sideL(t)=1;
    else
        second_sideL(t)=0;
    end
    
    for y=1:length(temp_pos_code250{t})
        if last_pos_code{t}(y)==91
            lastLreg250(t,y)=calc_rank(spatial_values{t}(1,1));
        elseif last_pos_code{t}(y)==92
            lastLreg250(t,y)=calc_rank(spatial_values{t}(2,1));
        elseif last_pos_code{t}(y)==93
            lastRreg250(t,y)=calc_rank(spatial_values{t}(1,2));
        elseif last_pos_code{t}(y)==94
            lastRreg250(t,y)=calc_rank(spatial_values{t}(2,2));
        end
    end
    
    if position_code{t}(1)==91|position_code{t}(1)==92
        first_side(t)=1; %first side is left
    else
        first_side(t)=0; %first side is right
    end

    %first picture value
    if position_code{t}(1)==91
        first_picture(t)=spatial_values{t}(1,1);
    elseif position_code{t}(1)==92
        first_picture(t)=spatial_values{t}(1,2);
    elseif position_code{t}(1)==93
        first_picture(t)=spatial_values{t}(2,1);
    elseif position_code{t}(1)==94
        first_picture(t)=spatial_values{t}(2,2);
    end
    
    if rem(first_picture(t),0.1)~=0 %easy check whether it is a prob/mag pic
        M_first(t)=1; %Magnitude first - 0.15 to 0.95
    else
        M_first(t)=0; %prob first - 0.1 to 0.9
    end

    %Next section works out if I saccade horizontally/diagonally on
    %option trials, and vertically/diagonally on attribute trials.
    if MCoditionNum(t)==3 %option trials
        vert(t)=NaN; %do I make a vertical saccade - nan in option trial
        if num_pics_viewed(t)>2
            if position_code{t}(2)==91||position_code{t}(2)==93
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
            end
            
        else
            horz(t)=NaN;
        end

    else %attribute trials
        horz(t)=NaN;
        if num_pics_viewed(t)>2
            if Lreg(t,3)>Rreg(t,3)
                third_sideL(t)=1;
            else
                third_sideL(t)=0;
            end

            if position_code{t}(2)==91||position_code{t}(2)==92
                if position_code{t}(3)==92||position_code{t}(3)==91
                    vert(t)=1;
                else
                    vert(t)=0;
                end
                
            elseif position_code{t}(2)==93||position_code{t}(2)==94
                if position_code{t}(3)==94||position_code{t}(3)==93
                    vert(t)=1;
                else
                    vert(t)=0;
                end
            end
            
        else
            vert(t)=NaN;
            third_sideL(t)=NaN;
        end
    end

    first_picture_rank=Lreg(:,1)+Rreg(:,1);
    second_picture_rank=Lreg(:,2)+Rreg(:,2);
    third_picture_rank=Lreg(:,3)+Rreg(:,3);
    
    %LH:
    first_position_code(t) = position_code{t}(1);
    second_position_code(t) = position_code{t}(2);
    
    top_first(t) = first_position_code(t)==91|first_position_code(t)==93; %was top presented first?
    prob_first(t) = xor(top_first(t),(Prob_top(t)==2)); %was probability presented first?

    if MCoditionNum(t)==3 %option trials
        if (first_position_code(t)==91)||(first_position_code(t)==93)
            top_picture_rank_opt(t) = first_picture_rank(t);
            bottom_picture_rank_opt(t) = second_picture_rank(t);
        elseif (first_position_code(t)==92)||(first_position_code(t)==94)
            top_picture_rank_opt(t) = second_picture_rank(t);
            bottom_picture_rank_opt(t) = first_picture_rank(t);
        else
            error;
        end
        
        if prob_first(t)
            prob_picture_rank_opt(t) = first_picture_rank(t);
            mag_picture_rank_opt(t) = second_picture_rank(t);
        else
            mag_picture_rank_opt(t) = first_picture_rank(t);
            prob_picture_rank_opt(t) = second_picture_rank(t);            
        end
        
        
        left_picture_rank_att(t) = nan;
        right_picture_rank_att(t) = nan;
        
        if Third_Cue_Position(t)==1||Third_Cue_Position(t)==3
            saccaded_up_opt(t) = 1;
            saccaded_prob_opt(t) = ~xor(top_first(t),prob_first(t));
        elseif Third_Cue_Position(t)==2||Third_Cue_Position(t)==4
            saccaded_up_opt(t) = 0;
            saccaded_prob_opt(t) = xor(top_first(t),prob_first(t));
        else
            saccaded_up_opt(t) = nan;
            saccaded_prob_opt(t) = nan;
        end
        saccaded_left_att(t) = nan;
     
    elseif MCoditionNum(t)==4 %option trials
        if (first_position_code(t)==91)|(first_position_code(t)==92)
            left_picture_rank_att(t) = first_picture_rank(t);
            right_picture_rank_att(t) = second_picture_rank(t);
        elseif (first_position_code(t)==93)|(first_position_code(t)==94)
            left_picture_rank_att(t) = second_picture_rank(t);
            right_picture_rank_att(t) = first_picture_rank(t);
        else
            error;
        end
        top_picture_rank_opt(t) = nan;
        bottom_picture_rank_opt(t) = nan;
                
        if Third_Cue_Position(t)==1|Third_Cue_Position(t)==2
            saccaded_left_att(t) = 1;
        elseif Third_Cue_Position(t)==3|Third_Cue_Position(t)==4
            saccaded_left_att(t) = 0;
        else
            saccaded_left_att(t) = nan;
        end
        saccaded_up_opt(t) = nan;
        saccaded_prob_opt(t) = nan;
    else
        left_picture_rank_att(t) = nan;
        right_picture_rank_att(t) = nan;
        top_picture_rank_opt(t) = nan;
        bottom_picture_rank_opt(t) = nan;
        saccaded_left_att(t) = nan;
        saccaded_up_opt(t) = nan;
        saccaded_prob_opt(t) = nan;
    end
end

%some useful variables w.r.t. choice, whether a side was looked at, etc.
first_side_chosen=first_side==Left_choice;
second_side_chosen=second_sideL==Left_choice;
third_side_chosen=third_sideL==Left_choice;

Left_side_seen=sum(Lreg'~=0)==2;
Right_side_seen=sum(Rreg'~=0)==2;

third_side(1:length(Left_EV))=NaN;
third_side(Lreg(:,3)~=0)=1;
third_side(Rreg(:,3)~=0)=0;


%% Choice logistic Regression: do subjects assign similar weight to probability and magnitude?

clear cond tmep beh_reg_mat_LR
LM_seen250(rem(LM_seen250,0.1)~=0)=LM_seen250(rem(LM_seen250,0.1)~=0)-0.05;
RM_seen250(rem(RM_seen250,0.1)~=0)=RM_seen250(rem(RM_seen250,0.1)~=0)-0.05;

LM_seen250(LM_seen250==0)=0.5;
RM_seen250(RM_seen250==0)=0.5;
LP_seen250(LP_seen250==0)=0.5;
RP_seen250(RP_seen250==0)=0.5;

tmep=[ones(length(Left_EV),1),LP_seen250'-RP_seen250',LM_seen250'-RM_seen250',first_side'];

cond(1:length(Left_EV),8)=0;
cond(MCoditionNum==3,[1:4])=1;
cond(MCoditionNum==4,[5:8])=1;

beh_reg_mat_LR=[tmep,tmep].*cond;

[out.choice_b,~,out.choice_b_stats] = glmfit(beh_reg_mat_LR,Left_choice','binomial','link','logit','constant','off');

%% Choice logistic Regression: do subjects assign similar weight to different cues?
rescale_vec = [0 0.1:0.2:0.9];
lastLreg_rescale = rescale_vec(lastLreg+1);
lastRreg_rescale = rescale_vec(lastRreg+1);
tmep = [ones(length(Left_EV),1),lastLreg_rescale-lastRreg_rescale];
[out.seq_choice_b,~,out.seq_choice_b_stats] = glmfit(tmep,Left_choice','binomial','link','logit','constant','off');

%% average number of pictures viewed on option/attribute trials, sorted by Pic1/Pic2

for i = 1:5
    for j = 1:5
        nPicView12Opt{i,j} = num_pics_viewed250(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==3);
        nPicView12Opt_chose1{i,j} = num_pics_viewed250(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==3&chosen_target'==1);
        nPicView12Opt_chose2{i,j} = num_pics_viewed250(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==3&chosen_target'==2);
        nPicView12Att{i,j} = num_pics_viewed250(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==4);
        nPicViewLRAtt{i,j} = num_pics_viewed250(left_picture_rank_att==i&...
            right_picture_rank_att==j&MCoditionNum==4);
    end
end
nPicView12Opt = cellfun(@mean,nPicView12Opt);
nPicView12Opt_chose1 = cellfun(@mean,nPicView12Opt_chose1);
nPicView12Opt_chose2 = cellfun(@mean,nPicView12Opt_chose2);
nPicView12Att = cellfun(@mean,nPicView12Att);
nPicViewLRAtt = cellfun(@mean,nPicViewLRAtt);
out.nPicView12Opt = nPicView12Opt;
out.nPicView12Opt_chose1 = nPicView12Opt_chose1;
out.nPicView12Opt_chose2 = nPicView12Opt_chose2;
out.nPicView12Att = nPicView12Att;
out.nPicViewLRAtt = nPicViewLRAtt;

%% average choice on option/attribute trials, sorted by Pic1/Pic2

for i = 1:5
    for j = 1:5
        AvChoice12Opt{i,j} = first_side_chosen(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==3);
        AvChoice12Att{i,j} = first_side_chosen(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==4);
    end
end
AvChoice12Opt = cellfun(@mean,AvChoice12Opt);
AvChoice12Att = cellfun(@mean,AvChoice12Att);
out.AvChoice12Opt = AvChoice12Opt;
out.AvChoice12Att = AvChoice12Att;

%% average saccade direction on option/attribute trials, sorted by Pic1/Pic2

for i = 1:5
    for j = 1:5
        AvSaccDir12Opt{i,j} = horz(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==3&num_pics_viewed'>2);
        AvSaccDir12Att{i,j} = vert(first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==4&num_pics_viewed'>2);
    end
end
AvSaccDir12Opt = cellfun(@mean,AvSaccDir12Opt);
AvSaccDir12Att = cellfun(@mean,AvSaccDir12Att);
out.AvSaccDir12Opt = AvSaccDir12Opt;
out.AvSaccDir12Att = AvSaccDir12Att;

%% average saccade direction on option/attribute trials, sorted by Left/Right or Up/Down

for i = 1:5
    for j = 1:5
        AvSaccDirUDOpt{i,j} = saccaded_up_opt(top_picture_rank_opt==i&...
            bottom_picture_rank_opt==j&MCoditionNum==3&num_pics_viewed>2);
        AvSaccDirLRAtt{i,j} = saccaded_left_att(left_picture_rank_att==i&...
            right_picture_rank_att==j&MCoditionNum==4&num_pics_viewed>2);
    end
end
AvSaccDirUDOpt = cellfun(@mean,AvSaccDirUDOpt);
AvSaccDirLRAtt = cellfun(@mean,AvSaccDirLRAtt);
out.AvSaccDirUDOpt = AvSaccDirUDOpt;
out.AvSaccDirLRAtt = AvSaccDirLRAtt;

%% average saccade direction on option trials, sorted by Prob/Mag

for i = 1:5
    for j = 1:5
        AvSaccDirMPOpt{i,j} = [horz(M_first'&first_picture_rank==i&...
            second_picture_rank==j&MCoditionNum'==3&num_pics_viewed'>2) ...
            (1-horz(~M_first'&first_picture_rank==j&...
            second_picture_rank==i&MCoditionNum'==3&num_pics_viewed'>2))];
    end
end

AvSaccDirMPOpt = cellfun(@mean,AvSaccDirMPOpt);
out.AvSaccDirMPOpt = AvSaccDirMPOpt;

%% accuracy and early choices - reported in main text
out.accuracy = mean(((Left_EV>=Right_EV)&Left_choice)|((Right_EV>=Left_EV)&~Left_choice));
out.earlychoices = 1 - mean(num_pics_viewedprezero==4); %average frequency of moving the joystick prior to seeing all four pictures

%% what if we redefine this according to the pictures seen by the subject? - asked for by reviewer
Left_EV_seen = LM_seen250.*LP_seen250;
Right_EV_seen = RM_seen250.*RP_seen250;
out.accuracy_seen = mean(((Left_EV_seen>=Right_EV_seen)&Left_choice)|...
                        ((Right_EV_seen>=Left_EV_seen)&~Left_choice));
                    
for i = 1:max(Sess_ind)
    %what if we look at option vs. attribute trials?
    at = MCoditionNum(:)==4; %attribute
    ot = MCoditionNum(:)==3; %option
    
    cc = at'&Sess_ind==i; %trials of interest
    out.accuracy_seen_att(i) = mean(((Left_EV_seen(cc)>=Right_EV_seen(cc))&Left_choice(cc))|...
        ((Right_EV_seen(cc)>=Left_EV_seen(cc))&~Left_choice(cc)));
    cc = ot'&Sess_ind==i; %trials of interest
    out.accuracy_seen_opt(i) = mean(((Left_EV_seen(cc)>=Right_EV_seen(cc))&Left_choice(cc))|...
        ((Right_EV_seen(cc)>=Left_EV_seen(cc))&~Left_choice(cc)));
    
    %what about the number of pictures viewed?
    npv = num_pics_viewed250;
    cc = npv==2 &Sess_ind==i;
    out.accuracyNPV(i,1) = mean(((Left_EV_seen(cc)>=Right_EV_seen(cc))&Left_choice(cc))|...
        ((Right_EV_seen(cc)>=Left_EV_seen(cc))&~Left_choice(cc)));
    out.mvd(i,1) = mean(abs(Left_EV(cc)-Right_EV(cc))); %mean value diff.
    
    cc = npv==3 &Sess_ind==i;
    out.accuracyNPV(i,2) = mean(((Left_EV_seen(cc)>=Right_EV_seen(cc))&Left_choice(cc))|...
        ((Right_EV_seen(cc)>=Left_EV_seen(cc))&~Left_choice(cc)));
    out.mvd(i,2) = mean(abs(Left_EV(cc)-Right_EV(cc))); %mean value diff.
    cc = npv==4 &Sess_ind==i;
    out.accuracyNPV(i,3) = mean(((Left_EV_seen(cc)>=Right_EV_seen(cc))&Left_choice(cc))|...
        ((Right_EV_seen(cc)>=Left_EV_seen(cc))&~Left_choice(cc)));
    out.mvd(i,3) = mean(abs(Left_EV(cc)-Right_EV(cc))); %mean value diff.
end


%% what if we look at trial t-1 (reward, choice) - asked for by reviewer

rlt = find([0 Rewarded(1:end-1)]); %rewarded last trial
out.accuracy_reward_Tminus1 = mean(((Left_EV_seen(rlt)>=Right_EV_seen(rlt))&Left_choice(rlt))|((Right_EV_seen(rlt)>=Left_EV_seen(rlt))&~Left_choice(rlt)));

nrlt = find([0 1-Rewarded(1:end-1)]); %not rewarded last trial
out.accuracy_noreward_Tminus1 = mean(((Left_EV_seen(nrlt)>=Right_EV_seen(nrlt))&Left_choice(nrlt))|((Right_EV_seen(nrlt)>=Left_EV_seen(nrlt))&~Left_choice(nrlt)));

cllt_stay =   find([0 Left_choice(1:end-1)] & Left_EV_seen>=Right_EV_seen); %chose left last trial - should choose left
cllt_switch = find([0 Left_choice(1:end-1)] & Right_EV_seen>=Left_EV_seen); %chose left last trial - should choose right
crlt_stay =   find([0 ~Left_choice(1:end-1)] & Right_EV_seen>=Left_EV_seen); %chose right last trial - should choose right
crlt_switch = find([0 ~Left_choice(1:end-1)] & Left_EV_seen>=Right_EV_seen); %chose right last trial - should choose left

out.accuracy_cllt_stay   = mean(Left_choice(cllt_stay)); 
out.accuracy_cllt_switch = mean(~Left_choice(cllt_switch)); 
out.accuracy_crlt_stay   = mean(~Left_choice(crlt_stay)); 
out.accuracy_crlt_switch = mean(Left_choice(crlt_switch)); 

%accuracy as function of chose left/right last trial
out.accuracy_cllt = (out.accuracy_cllt_stay +out.accuracy_cllt_switch)/2;
out.accuracy_crlt = (out.accuracy_crlt_stay +out.accuracy_crlt_switch)/2;

%accuracy as function of stay/switch from previous action
out.accuracy_stay = (out.accuracy_cllt_stay +out.accuracy_crlt_stay)/2;
out.accuracy_switch = (out.accuracy_cllt_switch +out.accuracy_crlt_switch)/2;


%% stats for left/right choices, first/second choices, for main fig 1

trind = MCoditionNum==4;
dm = [left_picture_rank_att-right_picture_rank_att];
[attchoice_b,~,attchoice_b_stats] = glmfit(dm(:,trind)',Left_choice(trind)','binomial','link','logit');

out.attchoice_b = attchoice_b;
out.attchoice_b_stats = attchoice_b_stats;

%
trind = MCoditionNum==3;
dm = [(first_picture_rank-3)+(second_picture_rank-3)]';
[optchoice_b,~,optchoice_b_stats] = glmfit(dm(:,trind)',first_side_chosen(trind)','binomial','link','logit');

out.optchoice_b = optchoice_b;
out.optchoice_b_stats = optchoice_b_stats;

%
trind = MCoditionNum==3;
dm = [(first_picture_rank-3)+(second_picture_rank-3)]';
dm = [ones(size(dm)); dm; dm.^2];
[optinfo_b,~,optinfo_t_stats] = ols(num_pics_viewed250(trind)',dm(:,trind)',eye(3));
out.optinfo_b = optinfo_b;
out.optinfo_t_stats = optinfo_t_stats;

%
trind = MCoditionNum==4;
dm = [left_picture_rank_att-right_picture_rank_att];
dm = [ones(size(dm)); dm; dm.^2];
[attinfo_b,~,attinfo_t_stats] = ols(num_pics_viewed250(trind)',dm(:,trind)',eye(3));
out.attinfo_b = attinfo_b;
out.attinfo_t_stats = attinfo_t_stats;

%%

out.totalNtrials = length(Left_choice);

