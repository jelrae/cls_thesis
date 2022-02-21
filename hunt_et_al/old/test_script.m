% Trying to rip whats useful to get the specific information we want from
% the dataset

tr = 1;
inf_trials=find(BhvInfo.ConditionNumber==3|BhvInfo.ConditionNumber==4);
trial_error=BhvInfo.TrialError; % 0 = 'Correct Response' (chose side with higher expected value),
inf_responded_trials=inf_trials(trial_error(inf_trials)==0|trial_error(inf_trials)==6);
ct=inf_responded_trials(tr);
nInfoGatheringTrials = length(inf_responded_trials); %number of InfoGathering Trials
picture_codes{tr}=BhvInfo.PlexonStrobes{ct}(find(BhvInfo.PlexonStrobes{ct}(:,1)>90&BhvInfo.PlexonStrobes{ct}(:,1)<95),:);
picture_codes{tr}(:,2)=picture_codes{tr}(:,2)-BhvInfo.PlexonStrobes{ct}(1,2);

%was first side left (1) or right (0)?
if picture_codes{tr}(1,1)==91||picture_codes{tr}(1,1)==92,
    First_side(tr)=1;
else
    First_side(tr)=0;
end
