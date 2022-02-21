function out = InfoGathering_Neurons_Cue123_alltrials_attentional_value_final(dosave,doproj)

% out = InfoGathering_Neurons_Cue123_alltrials_attentional_value_final(dosave,doproj)
%
% Runs regression of Cue 1,2,3,4 value, belief confirmation, on firing
%   rates timelocked to cue 1, 2, 3 presentation and response 
%
% if dosave == 1 (default 0), stores results in
% neuronal_regression_results/info_gathering_cue123_alltrials_attentional_final.mat
%
% Key output variables include:
%   significant_neuron (nUnits*nContrasts) - whether a unit significantly
%       encodes this contrast of parameter estimates, at p<0.05 corrected for
%       multiple comparisons across time
%   cpd1,2,3,R (nUnits*nRegressors*nTimeBins) - the coefficient of partial
%       determination for 18 different regressors of interest, time-locked 
%       to cue 1, 2, 3, response
%   t_stats1,2,3,R - cell array of units, each nTimeBins*nContrasts - sliding
%       T-statistic for each unit for each contrast of parameter estimates
%       across time, time-locked to cue 1, 2, 3, response
%
% see comments inside this function for meaning of each of the contrasts of
%   parameter estimates (lines 123 onwards) and CPD variables (lines 50
%   onwards)

% save results after running regression?
if nargin<1 
    dosave = 0;
end

% do regression on lots of different subselections of trials? (reviewer 3 comment)
if nargin<2 
    doproj = 0;
end

%% find pre-processed data files

bd = get_basedir;
addpath(genpath(bd));

datadir = fullfile(bd,'neuronal_data','LTH_processed_data');
filenames=dir(fullfile(datadir,'Unit*.mat'));

%% loop over units
for u=1:length(filenames)

    fprintf('Unit %0.0f/%0.0f\n',u,length(filenames)); %display current neuron
    fname = fullfile(datadir,sprintf('Unit%04.0f.mat',u));
    load(fname,'reg_mat','condition','attendsac1','Left_choice','Cue1_MatrixRaw',...
        'Cue2_MatrixRaw','Cue3_MatrixRaw','Response_MatrixRaw','num_pics_viewed',...
        'pre_cue','post_cue');


    %% build design matrix
    
    %reminder of variable meanings: 
    % reg_mat    = demeaned regressor of (1st, 2nd, 3rd, 4th) picture cue
    %             values (ranked -2 for worst to +2 for best)
    % condition = whether current trial is option (3) or attribute (4) trial
    % attendsac1= on attribute trials, is third saccade directed back
    %             towards first viewed option (1) or second viewed option (-1)
    
    %EVs 1-4: cue values 1-4 on option trials:
    regression_matrix(:,1:4)=reg_mat.*repmat(condition'==3,1,4); 
    
    %EVs 5-6: cue values 1-2 on attribute trials
    regression_matrix(:,5:6)=reg_mat(:,1:2).*repmat(condition'==4,1,2);
        
    %EVs 7-8: cue values 3-4 on attribute trials where subject attended option 1 at cue 3
    regression_matrix(:,7:8) = reg_mat(:,3:4).*repmat(attendsac1'==1,1,2);
       
    %EVs 9-10: cue values 3-4 on attribute trials where subject attended option 2 at cue 3
    regression_matrix(:,9:10) = reg_mat(:,3:4).*repmat(attendsac1'==-1,1,2);

    %EVs 11-12: indicator variable for option trials and attribute trials
    regression_matrix(:,11:12) = [condition'==3 condition'==4];
    
    %EV 13: belief confirmation at cue 2, option trials
    regression_matrix(:,13)= ((regression_matrix(:,1)>0)-(regression_matrix(:,1)<0)).*regression_matrix(:,2);
    
    %EV 14: belief confirmation at cue 2, attribute trials
    regression_matrix(:,14)= ((regression_matrix(:,5)<0)-(regression_matrix(:,5)>0)).*regression_matrix(:,6);
    
    %EV 15: belief confirmation at cue 3, option trials
    regression_matrix(:,15) = (((regression_matrix(:,1)+regression_matrix(:,2))<0)-...
        ((regression_matrix(:,1)+regression_matrix(:,2))>0)).*regression_matrix(:,3);
    
    %EV 16: belief confirmation at cue 3, attribute trials
    regression_matrix(:,16) = (((regression_matrix(:,5)-regression_matrix(:,6))<0)-...
        ((regression_matrix(:,5)-regression_matrix(:,6))>0)).*regression_matrix(:,9)+... %when option 2 was attended
        (((regression_matrix(:,5)-regression_matrix(:,6))>0)-...
        ((regression_matrix(:,5)-regression_matrix(:,6))<0)).*regression_matrix(:,7); %when option 1 was attended
    
    %Left_choice is currently valued 1 for left and 0 for right, but we want it to contrast left and right choices:
    Left_choice(Left_choice==0) = -1;

    %EV 17-18: left minus right choice on option and attribute trials
    regression_matrix(:,17:18) = [Left_choice'.*(condition'==3) Left_choice'.*(condition'==4)];
    
    DM_crosscorr(:,:,u) = corrcoef(regression_matrix);

    %% now pick out odd trials and even trials, to allow for cross-validation in subsequent analyses
    odd_trials = logical(zeros(size(regression_matrix,1),1));
    odd_trials(1:2:end) = 1;
    even_trials = logical(zeros(size(regression_matrix,1),1));
    even_trials(2:2:end) = 1;
    
    if doproj
        %repeated with random subsets of the data, rahter than just odd and
        %even, as requested by reviewer 3
        nProj = 100; %number of projections
        randproj_trials{u} = logical(zeros(size(regression_matrix,1),nProj));
        randproj_trials{u}(rand(size(randproj_trials{u}))>0.5) = 1;
        randproj_trials_inv{u} = ~randproj_trials{u};
    end
    
    %% only analyse trials where >=3 pictures have been viewed
    regression_matrix(num_pics_viewed<3,:)=[];
    odd_trials(num_pics_viewed<3) = [];
    even_trials(num_pics_viewed<3) = [];
    Cue1_MatrixRaw(num_pics_viewed<3,:)=[];
    Cue2_MatrixRaw(num_pics_viewed<3,:)=[];
    Cue3_MatrixRaw(num_pics_viewed<3,:)=[];
    Response_MatrixRaw(num_pics_viewed<3,:)=[];
    if doproj
        randproj_trials{u}(num_pics_viewed<3,:) = [];
        randproj_trials_inv{u}(num_pics_viewed<3,:) = [];
    end
    
    %% build contrast matrix
    
    contrasts1=eye(18);                                               % Initial regressors
    
    contrasts1(19:20,:) = contrasts1(1:2,:) + contrasts1(5:6,:);      % 19/20 - sum of VCue1/VCue2 on option/attribute
    contrasts1(21:22,:) = contrasts1(7:8,:) + contrasts1(9:10,:);     % 21/22 - sum of VCue3/Vue4 on attribute
    contrasts1(23,:) = contrasts1(13,:) + contrasts1(14,:);           % 23    - belief confirmation at cue 2
    contrasts1(24:25,:) = contrasts1(7:8,:) - contrasts1(9:10,:);     % 24/25 - saccade-contingent preference for 
                                                                      %         VCue3/VCue4 on attribute trials
    contrasts1(26:27,:) = contrasts1(1:2,:) -(contrasts1(5:6,:));     % 26/27 - preference for VCue1/2 on option>attribute
    contrasts1(28,:) = contrasts1(13,:) - contrasts1(14,:);           % 28    - preference for belief confirmation 
                                                                      %         on option>attribute
    contrasts1(29,:) = contrasts1(17,:) + contrasts1(18,:);           % 29    - left choice preference
    contrasts1(30,:) = contrasts1(15,:) + contrasts1(16,:);           % 30    - belief confirmation at cue 2
    
    %looking for signatures of value comparison at single neuron level
    contrasts1(31,:) = contrasts1(6,:) - contrasts1(5,:);               % 31    - VCue2 - VCue1 (attribute)
    contrasts1(32,:) = contrasts1(3,:) - contrasts1(1,:);               % 32    - VCue3 - VCue1 (option)
    contrasts1(33,:) = contrasts1(3,:) - contrasts1(2,:);               % 33    - VCue3 - VCue2 (option)
    contrasts1(34,:) = contrasts1(1,:) + contrasts1(2,:);               % 34    - VCue1 + VCue2 (option)
    contrasts1(35,:) = contrasts1(5,:) + contrasts1(6,:);               % 35    - VCue1 + VCue2 (attribute)

    %% variables to pass into regression function
    
    window_size=200;
    slide_width=10;
    epoch_duration=length(Cue1_MatrixRaw(1,:));
    pre_cue=pre_cue;
    num_perms=1000;
    
    variables(1)=window_size;
    variables(2)=slide_width;
    variables(3)=epoch_duration;
    variables(4)=pre_cue;
    variables(5)=num_perms;
    
    %note that we want to permute the design matrix the *same way* for
    %regression at different cues - we therefore construct permuted design
    %matrices to pass into the regression funciton below
    for p = 1:num_perms; perm_mat(:,p) = randperm(size(regression_matrix,1)); end
    
    %% Cue 1 regression

    [betas1{u},t_stats1{u},min_contrast_thresholds1{u},max_contrast_thresholds1{u},c_first_bin]=...
        sliding_regression_ols_perm(regression_matrix,Cue1_MatrixRaw,contrasts1,variables,perm_mat);
    [betas1_even{u},t_stats1_even{u}]=...
        sliding_regression_ols_perm(regression_matrix(even_trials,:),...
        Cue1_MatrixRaw(even_trials,:),contrasts1,variables);
    [betas1_odd{u},t_stats1_odd{u}]=...
        sliding_regression_ols_perm(regression_matrix(odd_trials,:),...
        Cue1_MatrixRaw(odd_trials,:),contrasts1,variables);
    if doproj
        for k = 1:nProj
            [~,t_stats1_proj(:,:,k)]=...
                sliding_regression_ols_perm(regression_matrix(randproj_trials{u}(:,k),:),...
                Cue1_MatrixRaw(randproj_trials{u}(:,k),:),contrasts1,variables);
        end
    end
    
    timebins = (-pre_cue-1+window_size/2):slide_width:post_cue-window_size/2;
    
    %calculate fraction of significant neurons
    for i=1:length(t_stats1{u}(1,:))
        log_mat(:,i)=t_stats1{u}(c_first_bin:end,i)>=max_contrast_thresholds1{u}(i)|t_stats1{u}(c_first_bin:end,i)<=...
            min_contrast_thresholds1{u}(i);
        time_encoding_matrix1(u,i)=sum(log_mat(:,i)); % total number of bins where neuron was 
                                                      % significantly encoding this contrast
        tmp=find(log_mat(:,i)>0);
        if ~isempty(tmp)
            latency_encoding_matrix1(u,i)=tmp(1);     % latency of bin when neuron was first encoding this contrast
        else
            latency_encoding_matrix1(u,i)=0;
        end
    end
    
    significant_neuron1(u,:)=time_encoding_matrix1(u,:)>0;
    
    %calculate coefficient of partial determination
    cpd1(u,:,:) =sliding_regression_CPD(regression_matrix,Cue1_MatrixRaw,variables,perm_mat);

    %% Cue 2 regression
    
    [betas2{u},t_stats2{u},min_contrast_thresholds2{u},max_contrast_thresholds2{u},c_first_bin]=...
        sliding_regression_ols_perm(regression_matrix,Cue2_MatrixRaw,contrasts1,variables,perm_mat);
    [betas2_even{u},t_stats2_even{u}]=...
        sliding_regression_ols_perm(regression_matrix(even_trials,:),...
        Cue2_MatrixRaw(even_trials,:),contrasts1,variables);
    [betas2_odd{u},t_stats2_odd{u}]=...
        sliding_regression_ols_perm(regression_matrix(odd_trials,:),...
        Cue2_MatrixRaw(odd_trials,:),contrasts1,variables);
    if doproj
        for k = 1:nProj
            [~,t_stats2_proj(:,:,k)]=...
                sliding_regression_ols_perm(regression_matrix(randproj_trials{u}(:,k),:),...
                Cue2_MatrixRaw(randproj_trials{u}(:,k),:),contrasts1,variables);
        end
    end
    timebins = (-pre_cue-1+window_size/2):slide_width:post_cue-window_size/2;
    
    %calculate fraction of significant neurons
    for i=1:length(t_stats2{u}(1,:))
        log_mat(:,i)=t_stats2{u}(c_first_bin:end,i)>=max_contrast_thresholds2{u}(i)|t_stats2{u}(c_first_bin:end,i)<=...
            min_contrast_thresholds2{u}(i);
        time_encoding_matrix2(u,i)=sum(log_mat(:,i)); % total number of bins where neuron was 
                                                      % significantly encoding this contrast
        tmp=find(log_mat(:,i)>0);
        if ~isempty(tmp)
            latency_encoding_matrix2(u,i)=tmp(1);     % latency of bin when neuron was first encoding this contrast
        else
            latency_encoding_matrix2(u,i)=0;
        end
    end
    
    significant_neuron2(u,:)=time_encoding_matrix2(u,:)>0;
    
    %calculate coefficient of partial determination
    cpd2(u,:,:) =sliding_regression_CPD(regression_matrix,Cue2_MatrixRaw,variables,perm_mat);
    
    %% Cue 3 regression
    [betas3{u},t_stats3{u},min_contrast_thresholds3{u},max_contrast_thresholds3{u},c_first_bin]=...
        sliding_regression_ols_perm(regression_matrix,Cue3_MatrixRaw,contrasts1,variables,perm_mat);
    [betas3_even{u},t_stats3_even{u}]=...
        sliding_regression_ols_perm(regression_matrix(even_trials,:),...
        Cue3_MatrixRaw(even_trials,:),contrasts1,variables);
    [betas3_odd{u},t_stats3_odd{u}]=...
        sliding_regression_ols_perm(regression_matrix(odd_trials,:),...
        Cue3_MatrixRaw(odd_trials,:),contrasts1,variables);
    if doproj
        for k = 1:nProj
            [~,t_stats3_proj(:,:,k)]=...
                sliding_regression_ols_perm(regression_matrix(randproj_trials{u}(:,k),:),...
                Cue3_MatrixRaw(randproj_trials{u}(:,k),:),contrasts1,variables);
        end
    end
    timebins = (-pre_cue-1+window_size/2):slide_width:post_cue-window_size/2;
    
    %calculate fraction of significant neurons
    for i=1:length(t_stats3{u}(1,:)),
        log_mat(:,i)=t_stats3{u}(c_first_bin:end,i)>=max_contrast_thresholds3{u}(i)|t_stats3{u}(c_first_bin:end,i)<=...
            min_contrast_thresholds3{u}(i);
        time_encoding_matrix3(u,i)=sum(log_mat(:,i)); % total number of bins where neuron was 
                                                      % significantly encoding this contrast
        tmp=find(log_mat(:,i)>0);
        if ~isempty(tmp)
            latency_encoding_matrix3(u,i)=tmp(1);     % latency of bin when neuron was first encoding this contrast
        else
            latency_encoding_matrix3(u,i)=0;
        end
    end
    significant_neuron3(u,:)=time_encoding_matrix3(u,:)>0;
    
    %calculate coefficient of partial determination
    cpd3(u,:,:) =sliding_regression_CPD(regression_matrix,Cue3_MatrixRaw,variables,perm_mat);
    
    %% Response regression
    [betasR{u},t_statsR{u},min_contrast_thresholdsR{u},max_contrast_thresholdsR{u},c_first_bin]=...
        sliding_regression_ols_perm(regression_matrix,Response_MatrixRaw,contrasts1,variables,perm_mat);
    [betasR_even{u},t_statsR_even{u}]=...
        sliding_regression_ols_perm(regression_matrix(even_trials,:),...
        Response_MatrixRaw(even_trials,:),contrasts1,variables);
    [betasR_odd{u},t_statsR_odd{u}]=...
        sliding_regression_ols_perm(regression_matrix(odd_trials,:),...
        Response_MatrixRaw(odd_trials,:),contrasts1,variables);
    if doproj
        for k = 1:nProj
            [~,t_statsR_proj(:,:,k)]=...
                sliding_regression_ols_perm(regression_matrix(randproj_trials{u}(:,k),:),...
                Response_MatrixRaw(randproj_trials{u}(:,k),:),contrasts1,variables);
        end
    end
    timebins = (-pre_cue-1+window_size/2):slide_width:post_cue-window_size/2;
    
    if doproj
        subspace_belief_proj(u,:) = squeeze((t_stats3_proj(41,15,:)+t_stats3_proj(41,16,:)+...
            t_stats2_proj(41,14,:)+t_stats2_proj(41,13,:))/4);
        subspace_RL_proj(u,:) = (squeeze(t_statsR_proj(30,17,:))+squeeze(t_statsR_proj(30,18,:)))/2;
        subspace_O_proj(u,:) = squeeze(t_statsR_proj(30,11,:));
        subspace_A_proj(u,:) = squeeze(t_statsR_proj(30,12,:));
    end
    
    for i=1:length(t_statsR{u}(1,:)),
        log_mat(:,i)=t_statsR{u}(c_first_bin:end,i)>=max_contrast_thresholdsR{u}(i)|t_statsR{u}(c_first_bin:end,i)<=...
            min_contrast_thresholdsR{u}(i);
        time_encoding_matrixR(u,i)=sum(log_mat(:,i)); % total number of bins where neuron was 
                                                      % significantly encoding this contrast
        tmp=find(log_mat(:,i)>0);
        if ~isempty(tmp)
            latency_encoding_matrixR(u,i)=tmp(1);     % latency of bin when neuron was first encoding this contrast
        else
            latency_encoding_matrixR(u,i)=0;
        end
    end
    significant_neuronR(u,:)=time_encoding_matrixR(u,:)>0;

    %calculate coefficient of partial determination
    cpdR(u,:,:)=sliding_regression_CPD(regression_matrix,Response_MatrixRaw,variables,perm_mat);

    %% tidy up
    clear regression_matrix contrasts1 perm_mat

end %of loop over units
clear t_stats*proj 

if dosave
    save(fullfile(bd,'neuronal_regression_results','info_gathering_cue123_alltrials_attentional_final.mat'))
end

%% tidy up output
out = whos;
C = cell(2,length(out));
C(1,:) = {out.name};
out = struct(C{:});
tmp = structvars(3,out,0);
for i = 1:size(tmp,1)
    eval(tmp(i,:));
end
