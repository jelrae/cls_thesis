clear all

des_data = 'post';
monkey = 'kurt';
data_used = 'all';

set_paths(des_data);

ft_defaults

%% Run script to get ROIs specific to the monkey
fig_6_ROIS;
% % Define region of interest 1
% roi = V1;
% %Define region of interest 2
% roi_2 = V4;
% roi_combine = [roi roi_2];
roi = all_channels;

%% Configuration parameters and bipolar merging
cfg = config_setup(roi, monkey, des_data, data_used, false);

%%  Parameters
equal = cfg.equal;
path_file = cfg.path;
monkey = cfg.monkey;
v_session = cfg.v_session;
numses = length(v_session);
baseline = cfg.bsln;

% clc
if strcmp(equal,'no')
    disp('No Trial uniformity')
else
    disp('Trial Uniformity Activated')
end
fprintf('loading %d sessions\n',numses);

%%
for ns=1:numses

    fprintf('session %d of %d\n',ns,numses);
    if baseline
        name_ses = sprintf('%s%s_%0.3d_baseline.mat',path_file,monkey,v_session(ns));
    else
        name_ses = sprintf('%s%s_%0.3d_rawdata.mat',path_file,monkey,v_session(ns));
    end
    
    coh_cfg.avgoverchancmb = 'yes';eval(sprintf('load %s;',name_ses));
    if ~isfield(data,'sampleinfo')   %otherwise get an error
        data.sampleinfo = [];
    end
    %update trl (add a new column to identify the session number)
   try
    tmp_trl = data.trialinfo;%findcfg(data.cfg,'trl');
   catch
       tmp_trl = data.cfg.trl(:,4:end);
   end
    tmp_trl(:,end+1) = ns;
    data.trialinfo = tmp_trl;
    numtr = size(data.trial,2);  

%     %% Remove mean
%     tmp_cfg = [];
%     tmp_cfg.demean = 'yes';
%     data = ft_preprocessing(tmp_cfg,data);


%     %%  Bipolar Derivation
%     
%     disp('Rereferencing the dataset')
%     if strcmp(cfg.bipolarder,'average')
%         tmp_cfg = [];
%         tmp_cfg.reref = 'yes';
%         tmp_cfg.refchannel = cfg.chan_list;%'all';
%     elseif strcmp(cfg.bipolarder,'ownmontage')
%         tmp_cfg = [];
%         tmp_cfg.reref = 'yes';
%         tmp_cfg.montage = cls_createmontage(data.label,cfg.bipolarder_list);
%         tmp_cfg.channel = data.label;
%     else
%         tmp_cfg = [];
%         tmp_cfg=load(cfg.bipolarder);
% %         tmp_cfg.reref = 'yes';
%     end
% 
%     data_bipolar = ft_preprocessing(tmp_cfg,data);
% 
%     if strcmp(cfg.bipolarder,'average') && ~any(strcmp(cfg.chan_list,'all'))
%         sel1 = match_str(data_bipolar.label,cfg.chan_list);
%         data_bipolar.label = data_bipolar.label(sel1);
%         numtr_ses = size(data_bipolar.trial,2);
%         for i=1:numtr_ses
%             fprintf('Trial %d\n',i)
%             data_bipolar.trial{i} = data_bipolar.trial{i}(sel1,:);
%         end
%     end
    
    %% normalization
%     if strcmp(cfg.norm,'yes')
%         s_cfg = [];
%         s_cfg.channelnormalise = 'yes';
%         data = ft_preprocessing(s_cfg,data_bipolar);
%         %data = session_norm2(data_bipolar);
%     else
%         data = data_bipolar;
%     end
% 
    eval(sprintf('data%0.2d = data;',ns));  
    clear data

end

%%  Merge all datasets
numtrperses = [];

if length(v_session)>1
    disp('concatenating trials')
    msg = 'data = ft_appenddata([]';
    for i=1:numses
        msg = sprintf('%s,data%0.2d',msg,i);
        eval(sprintf('numtrperses(i) = length(data%0.2d.trial);',i));
    end
    msg = sprintf('%s);',msg);
    
    
    
    eval(msg)
    for i=1:numses
        eval(sprintf('clear data%0.2d',i));
    end
else
    data = data01;
    clear data01;
   % numtrperses(ns) = length(data.trial);
end



%% Split into conditions
if ~baseline
    cnd = {'AttIn' 'AttOut'};
    trl = data.trialinfo;%findcfg(data.cfg,'trl');
    col= size(trl,2);
    
    if col==6
        col_in = 5;
    elseif col==7
        col_in = 5+1;
    end
    
    for a=1:length(cnd)
        index = find(trl(:,col_in)==a);
        newtrl = trl(index,:);
        tmp = data;
        tmp.trial = tmp.trial(index);
        tmp.time = tmp.time(index);
        tmp.trialinfo = newtrl;
        tmp.dimord = 'unknown';
        %tmp.EyeChan = tmp.EyeChan(index);
        eval(sprintf('dataN%s = tmp;',cnd{a}));
    end
    
    
    %% Equal Numbers of trials
    if strcmp(equal,'yes')
        disp('Calculating equal numbers of trials to each condition')
        
        trlAttIn  = dataNAttIn.trialinfo; %findcfg(dataAttIn.cfg,'trl');
        trlAttOut = dataNAttOut.trialinfo; %findcfg(dataAttOut.cfg,'trl');
        trlAttIn(:,end+1) = 0;
        trlAttOut(:,end+1) = 0;
        
        numses = trlAttIn(end,col_in+1);
        %dataNAttIn = dataAttIn;
        %dataNAttOut = dataAttOut;
        
        for i=1:numses
            
            numtrAttIn = trlAttIn(:,col_in+1) == i;
            numtrAttOut =trlAttOut(:,col_in+1) == i;
            tmp = sum(numtrAttIn)-sum(numtrAttOut);
            
            if tmp==0
                %do nothing
            elseif tmp>0
                
                tmptrl = [];
                tmptrl = trlAttIn(numtrAttIn,:);
                tmptrl(end+1-abs(tmp):end,end) = 1;
                trlAttIn(numtrAttIn,:) = tmptrl;
            elseif tmp<0
                
                tmptrl = [];
                tmptrl = trlAttOut(numtrAttOut,:);
                tmptrl(end+1-abs(tmp):end,end) = 1;
                trlAttOut(numtrAttOut,:) = tmptrl;
            end
        end
        
        
        
        
        dataNAttIn.trial(trlAttIn(:,end)==1) = [];
        dataNAttIn.time(trlAttIn(:,end)==1) = [];
        trlAttIn(trlAttIn(:,end)==1,:) = [];
        dataNAttIn.trialinfo = trlAttIn;
        %dataNAttIn.EyeChan(trlAttIn(:,10)==1) = [];
        
        
        dataNAttOut.trial(trlAttOut(:,end)==1) = [];
        dataNAttOut.time(trlAttOut(:,end)==1) = [];
        trlAttOut(trlAttOut(:,end)==1,:) = [];
        dataNAttOut.trialinfo = trlAttOut;
        %dataNAttOut.EyeChan(trlAttOut(:,10)==1) = [];
    end
    %%
    dataAttIn = dataNAttIn;
    dataAttOut = dataNAttOut;
    %%
    clear dataNAttIn dataNAttOut
else
    %% baseline calculation
    %we need to known how many trials per session; so we use mcfg.numtr to
    %that purpose
    
    dataAttIn = data;
    dataAttOut = [];
end

%Trim the data
COI = load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\kurt_montage_layout\kurt_COI_alongfiner2.mat');
COI = struct2cell(COI);
COI = COI{1};
all_cfg = [];
all_cfg.channel = permute(COI,[2,1]);
dataAttIn =  ft_selectdata(all_cfg, dataAttIn);
dataAttOut =  ft_selectdata(all_cfg, dataAttOut);

%plot the data
in_trials = dataAttIn.trial;
in_trials = cat(3,in_trials{:});
histogram(in_trials)
title('AttIn histogram');

%get some stats
mean_channel = mean(in_trials, [2,3]);
std_channel = std(in_trials, 0, [2,3]);

labels = dataAttIn.label;

%check whats missing
for i = 1:length(COI)
    if ~ismember(COI(i), labels)
        disp(COI(i));   
    end
end

% More stats

max_channel = max(in_trials, [],[2,3]);
min_channel = min(in_trials, [],[2,3]);

trial_mean = mean(in_trials, 2);
trial_std = std(in_trials,0, 2);

% A8L_regions = {'A08','A07','A09','A08','O04','O03','O06','O05','O07','O06','O08','O07','O02','O01','O03','O02'};
A8L_regions = {'A08','A07','O06','O05','O07','O02','O01','O03'};


power_check(dataAttIn, A8L_regions);

% %Data viewing in FT
% cfg = [];
% cfg.viewmode = 'vertical';
% artfct       = ft_databrowser(cfg, dataAttIn)