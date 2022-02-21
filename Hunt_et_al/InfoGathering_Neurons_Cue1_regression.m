function out = InfoGathering_Neurons_Cue1_regression(dosave)

% out = InfoGathering_Neruons_Cue1_regression(dosave)
%
% Runs regression of Attended Value, Left/Right Value, Prob/Mag. Value,
%   Top/Bottom Value on firing rates timelocked to Cue 1 presentation
%
% if dosave == 1 (default 0), stores results in
% neuronal_regression_results/info_gathering_cue1_both.mat
%
% Key output variables include:
%   significant_neuron (nUnits*nContrasts) - whether a unit significantly
%       encodes this contrast of parameter estimates, at p<0.05 corrected for
%       multiple comparisons across time
%   varCPD (nUnits*nRegressors*nTimeBins) - the coefficient of partial
%       determination for 8 different regressors of interest
%   t_stats1 - cell array of units, each nTimeBins*nContrasts - sliding
%       T-statistic for each unit for each contrast of parameter estimates
%       across time
%
% see comments inside this function for meaning of each of the contrasts of
%   parameter estimates (lines 186 onwards) and CPD variables (lines 296
%   onwards)

%%
% save results after running regression?
if nargin<1 
    dosave = 0;
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
    load(fname,'inf_responded_trials','Left_prob','Left_pay','Right_prob','Right_pay',...
        'First_pos','Second_pos','BhvInfo','Left_EV','Right_EV','Left_choice',...
        'First_side','pictures_views','Cue1_MatrixRaw','pre_cue','post_cue');

    
    %% compute some further trial variables for the regression
    for tr=1:length(inf_responded_trials) %loop over info-gathering trials
        ct=inf_responded_trials(tr);
        Prob_Top(tr)=BhvInfo.Uservars{ct}.ProbTop;

        if Prob_Top(tr)==1
            spatial_vals(tr,1)=Left_prob(tr);       % left-upper value
            spatial_vals(tr,2)=Left_pay(tr);        % left-lower value
            spatial_vals(tr,3)=Right_prob(tr);      % right-upper value
            spatial_vals(tr,4)=Right_pay(tr);       % right-lower value
            
            %now calculate 'P_first_val' - value of probability if P was
            %first, 0 otherwise - and likewise, 'M_first_val'
            if First_pos(tr)==1||First_pos(tr)==3
                P_first_val(tr)=spatial_vals(tr,First_pos(tr));
                M_first_val(tr)=0;
            else
                M_first_val(tr)=spatial_vals(tr,First_pos(tr));
                P_first_val(tr)=0;
            end
            
            %now calculate 'P_second_val' - value of probability if P was
            %second, 0 otherwise - and likewise, 'M_second_val'
            if Second_pos(tr)==2||Second_pos(tr)==4
                M_second_val(tr)=spatial_vals(tr,Second_pos(tr));
                P_second_val(tr)=0;
            else
                M_second_val(tr)=0;
                P_second_val(tr)=spatial_vals(tr,Second_pos(tr));
            end
        else
            spatial_vals(tr,1)=Left_pay(tr);        % left-upper value
            spatial_vals(tr,2)=Left_prob(tr);       % left-lower value
            spatial_vals(tr,3)=Right_pay(tr);       % right-upper value
            spatial_vals(tr,4)=Right_prob(tr);      % right-lower value
            
            %now calculate 'P_first_val' - value of probability if P was
            %first, 0 otherwise - and likewise, 'M_first_val'
            if First_pos(tr)==2||First_pos(tr)==4
                P_first_val(tr)=spatial_vals(tr,First_pos(tr));
                M_first_val(tr)=0;
            else
                P_first_val(tr)=0;
                M_first_val(tr)=spatial_vals(tr,First_pos(tr));
            end
            
            %now calculate 'P_second_val' - value of probability if P was
            %second, 0 otherwise - and likewise, 'M_second_val'
            if Second_pos(tr)==1||Second_pos(tr)==3
                M_second_val(tr)=spatial_vals(tr,Second_pos(tr));
                P_second_val(tr)=0;
            else
                M_second_val(tr)=0;
                P_second_val(tr)=spatial_vals(tr,Second_pos(tr));
            end
        end
        
    end %of loop over info gathering trials
    
    %% store key behaviour for each unit in a structure, behaviour
    behaviour{u}=struct;
    behaviour{u}.LP=Left_prob;
    behaviour{u}.RP=Right_prob;
    behaviour{u}.LM=Left_pay;
    behaviour{u}.RM=Right_pay;
    behaviour{u}.LEV=Left_EV;
    behaviour{u}.REV=Right_EV;
    behaviour{u}.choice=Left_choice;
    behaviour{u}.First_pos=First_pos;
    behaviour{u}.Second_pos=Second_pos;
    behaviour{u}.Prob_top=Prob_Top;
    behaviour{u}.First_picture_value=M_first_val+P_first_val;
    behaviour{u}.Second_picture_value=M_second_val+P_second_val;

    %% construct regressors for design matrix
    
    %M1LT will be value of Magnitude on trials where 1st cue was
    %Magnitude/Left/Top. Same for Left/Right, Top/Bottom,
    %Probability/Magnitude
    M1LT(1:length(First_side),1)=zeros;
    M1RT(1:length(First_side),1)=zeros;
    M1LB(1:length(First_side),1)=zeros;
    M1RB(1:length(First_side),1)=zeros;
    M1LT(First_pos==1)=M_first_val(First_pos==1);
    M1RT(First_pos==3)=M_first_val(First_pos==3);
    M1LB(First_pos==2)=M_first_val(First_pos==2);
    M1RB(First_pos==4)=M_first_val(First_pos==4);
    
    P1LT(1:length(First_side),1)=zeros;
    P1RT(1:length(First_side),1)=zeros;
    P1LB(1:length(First_side),1)=zeros;
    P1RB(1:length(First_side),1)=zeros;
    P1LT(First_pos==1)=P_first_val(First_pos==1);
    P1RT(First_pos==3)=P_first_val(First_pos==3);
    P1LB(First_pos==2)=P_first_val(First_pos==2);
    P1RB(First_pos==4)=P_first_val(First_pos==4);
    
    %value of cue 1 when on left, right, top, bottom:
    LVal=sum([M1LT, M1LB, P1LT, P1LB]');   
    RVal=sum([M1RT, M1RB, P1RT, P1RB]');
    TVal=sum([M1LT, M1RT, P1LT, P1RT]');  
    BVal=sum([M1LB, M1RB, P1LB, P1RB]');
    
    LRVal=LVal-RVal; %regressor for cells that prefer left value to right value
    TBVal=TVal-BVal; %regressor for cells that prefer top value to bottom value
    PMVal=P_first_val-M_first_val; %regressor for cells that prefer probability value to magnitude value
    Val=LVal+RVal; %regressor for cells that respond to currently attended value
    
    % Demean mag regressors
    M1LT(M1LT~=0)=M1LT(M1LT~=0)-0.55;
    M1RT(M1RT~=0)=M1RT(M1RT~=0)-0.55;
    M1LB(M1LB~=0)=M1LB(M1LB~=0)-0.55;
    M1RB(M1RB~=0)=M1RB(M1RB~=0)-0.55;
      
    % Demean prob regressors
    P1LT(P1LT~=0)=P1LT(P1LT~=0)-0.5;
    P1RT(P1RT~=0)=P1RT(P1RT~=0)-0.5;
    P1LB(P1LB~=0)=P1LB(P1LB~=0)-0.5;
    P1RB(P1RB~=0)=P1RB(P1RB~=0)-0.5;
    
    % Construct constant terms for each time a cue is presented in each
    % spatial position
    POS1(1:length(First_side),1)=zeros;
    POS2(1:length(First_side),1)=zeros;
    POS3(1:length(First_side),1)=zeros;
    POS4(1:length(First_side),1)=zeros;
    POS1(First_pos==1)=1;
    POS2(First_pos==2)=1;
    POS3(First_pos==3)=1;
    POS4(First_pos==4)=1;
    
    %lastly, include a regressor for the numbers of pictures viewed
    reg_target=demean(pictures_views'/4);
    
    %% build design matrix
    regression_matrix=[M1LT, M1RT, M1LB, M1RB, P1LT, P1RT, P1LB, P1RB, POS1, POS2, POS3, POS4, reg_target];
    DM_crosscorr(:,:,u) = corrcoef(regression_matrix); %correlation between different regressors, for this unit

    %% build contrast matrix
    contrasts1=eye(13);                                                 % Initial regressors
    
    contrasts1(14,:)=sum(contrasts1(1:2,:));                            % 14 Mag top
    contrasts1(15,:)=sum(contrasts1(3:4,:));                            % 15 Mag bottom
    contrasts1(16,:)=sum(contrasts1(5:6,:)) ;                           % 16 Prob top
    contrasts1(17,:)=sum(contrasts1(7:8,:));                            % 17 Prob bottom
    
    contrasts1(18,:)=sum(contrasts1([1,3],:)) ;                         % 18 Mag L
    contrasts1(19,:)=sum(contrasts1([2,4],:)) ;                         % 19 Mag R
    contrasts1(20,:)=sum(contrasts1([5,7],:));                          % 20 Prob L
    contrasts1(21,:)=sum(contrasts1([6,8],:));                          % 21 Prob R
    
    contrasts1(22,:)=sum(contrasts1([1,2,5,6],:));                      % 22 Top value
    contrasts1(23,:)=sum(contrasts1([3,4,7,8],:));                      % 23 bottom value
    contrasts1(24,:)=sum(contrasts1([1,3,5,7],:));                      % 24 Left value
    contrasts1(25,:)=sum(contrasts1([2,4,6,8],:));                      % 25 Right value
    contrasts1(26,:)=sum(contrasts1([1:8],:));                          % 26 Value
    
    contrasts1(27,:)=sum(contrasts1([9,10],:));                         % 27 Left pos
    contrasts1(28,:)=sum(contrasts1([11,12],:));                        % 28 Right pos
    contrasts1(29,:)=sum(contrasts1([9,11],:));                         % 29 Top pos
    contrasts1(30,:)=sum(contrasts1([10,12],:));                        % 30 Bottom pos
    
    contrasts1(31,:)=(contrasts1([1],:))-...                            % 31 Mag top L-R
        (contrasts1([2],:));
    contrasts1(32,:)=(contrasts1([3],:))-...                            % 32 Mag bottom L-R
        (contrasts1([4],:));
    contrasts1(33,:)=sum(contrasts1([1,3],:))-...                       % 33 Mag L-R
        sum(contrasts1([2,4],:));
    contrasts1(34,:)=sum(contrasts1([1:2],:))-...                       % 34 Mag T-B
        sum(contrasts1([3:4],:));
    contrasts1(35,:)=(contrasts1([5],:))-...                            % 35 Prob top L-R
        (contrasts1([6],:));
    contrasts1(36,:)=(contrasts1([7],:))-...                            % 36 prob Bottom L-R
        (contrasts1([8],:));
    contrasts1(37,:)=sum(contrasts1([5,7],:))-...                       % 37 prob L-R
        sum(contrasts1([6,8],:));
    contrasts1(38,:)=sum(contrasts1([5,6],:))-...                       % 38 prob T-B
        sum(contrasts1([7,8],:));
    contrasts1(39,:)=sum(contrasts1([1,5],:))-...                       % 39 value top L-R
        sum(contrasts1([2,6],:));
    contrasts1(40,:)=sum(contrasts1([3,7],:))-...                       % 40 value bottom L-R
        sum(contrasts1([4,8],:));
    contrasts1(41,:)=sum(contrasts1([1,3,5,7],:))-...                   % 41 Value L-R
        sum(contrasts1([2,4,6,8],:));
    contrasts1(42,:)=sum(contrasts1([1:2,5:6],:))-...                   % 42 Value T-B
        sum(contrasts1([3:4,7:8],:));
    contrasts1(43,:)=(contrasts1(9,:))-...                              % 43 Top Left - Bottom left
        (contrasts1(10,:));
    contrasts1(44,:)=(contrasts1(11,:))-...                             % 44 Top Right - Bottom Right
        (contrasts1(12,:));
    contrasts1(45,:)=(contrasts1(9,:))-...                              % 45 Top Left - top right
        (contrasts1(11,:));
    contrasts1(46,:)=(contrasts1(10,:))-...                             % 46 Bottom left - Bottom right
        (contrasts1(12,:));
    contrasts1(47,:)=sum(contrasts1([9:10],:))-...                      % 47 Left - Right
        sum(contrasts1([11:12],:));
    contrasts1(48,:)=sum(contrasts1([9,11],:))-...                      % 48 Top - Bottom
        sum(contrasts1([10,12],:));
    
    contrasts1(49,:)=sum(contrasts1(1:4,:));                            % 49 Mag
    contrasts1(50,:)=sum(contrasts1(5:8,:));                            % 50 Prob
    contrasts1(51,:)=sum(contrasts1(5:8,:))-sum(contrasts1(1:4,:));     % 51 Mag-Prob
    contrasts1(52,:)=sum(contrasts1([1,3],:))-sum(contrasts1([5,7],:)); % 52 LM-LP
    contrasts1(53,:)=sum(contrasts1([2,4],:))-sum(contrasts1([6,8],:)); % 53 RM-RP
    
    contrasts1(54,:)=(contrasts1([9],:))-sum(contrasts1([10,11,12],:)); % 54 POS1 over others
    contrasts1(55,:)=(contrasts1([10],:))-sum(contrasts1([9,11,12],:)); % 55 POS2 over others
    contrasts1(56,:)=(contrasts1([11],:))-sum(contrasts1([9,10,12],:)); % 56 POS3 over others
    contrasts1(57,:)=(contrasts1([12],:))-sum(contrasts1([9,10,11],:)); % 57 POS4 over others
    
    %% Run sliding ordinary least squares regression, including permutation testing
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

    [betas1{u},t_stats1{u},min_contrast_thresholds1{u},max_contrast_thresholds1{u},c_first_bin]=...
        sliding_regression_ols_perm(regression_matrix,Cue1_MatrixRaw,contrasts1,variables);
    
    timebins = (-pre_cue-1+window_size/2):slide_width:post_cue-window_size/2;
    
    %% work out when neurons are significant: 
    %  compare T-stats against null distribution from permutation test in
    %  sliding_regression_ols_perm.m
    
    for i=1:size(contrasts1,1)

        log_mat(:,i)=t_stats1{u}(c_first_bin:end,i)>=max_contrast_thresholds1{u}(i)|t_stats1{u}(c_first_bin:end,i)<=min_contrast_thresholds1{u}(i);
        time_encoding_matrix(u,i)=sum(log_mat(:,i)); % total number of bins where neuron was significantly encoding this contrast

        tmp=find(log_mat(:,i)>0);
        
        if ~isempty(tmp)
            latency_encoding_matrix(u,i)=tmp(1);     % latency of bin when neuron was first encoding this contrast
        else            
            latency_encoding_matrix(u,i)=0;            
        end        
    end
    
    significant_neuron(u,:)=time_encoding_matrix(u,:)>0;
    
    %% CPD
    window_size=variables(1);
    slide_width=variables(2);
    epoch_duration=variables(3);
    pre_cue=variables(4);
    num_perms=variables(5);
    totbins2=(epoch_duration/slide_width)-((window_size/slide_width)-1);
    BINS=[1:slide_width:(((epoch_duration)-slide_width)+1)];
    
    for thisbin=1:floor(totbins2)
        smoothed_data(:,thisbin)=mean(Cue1_MatrixRaw(:,BINS(thisbin):(BINS(thisbin)+(window_size-1))),2);
    end

    [varCPD(u,:,:)] = cpd(smoothed_data,[Val',LRVal',PMVal', TBVal', POS1, POS2, POS3, POS4]);
    % CPD is as follows:
    % 1 - value
    % 2 - left minus right value
    % 3 - probability minus magnitude
    % 4 - top minus bottom value
    % 5-8 constant terms for different positions on screen
        
    %% clean up
    clear M1LT M1RT M1LB M1RB P1LT P1RT P1LB P1RB POS1 POS2 POS3 POS4 reg_target ...
        Val LRVal PMVal TBVal LVal RVal TVal BVal P_first_val P_second_val ...
        M_first_val M_second_val spatial_vals smoothed_data Cue1_MatrixRaw First_pos ...
        First_side inf_responded_trials Left_choice Left_EV Left_pay ...
        Left_prob Prob_Top regression_matrix Right_EV Right_pay Right_prob ...
        Second_pos tr ct fname i pictures_views tmp

end %of loop over units

if dosave
    save(fullfile(bd,'neuronal_regression_results','info_gathering_cue1_both.mat'))
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
