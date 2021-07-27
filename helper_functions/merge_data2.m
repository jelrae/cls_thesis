function [dataAttIn, dataAttOut] = merge_data2(cfg, bipolar)

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

    %% Remove Eye Channel to a separate field
    if 0

        eye_ch(1) = match_str(data.label,'Eye_Hor');
        eye_ch(2) = match_str(data.label,'Eye_Ver');
        for i=1:numtr
            data.EyeChan{i} = data.trial{i}(eye_ch,:);
            data.trial{i}(eye_ch,:) = [];
        end

        data.label(eye_ch) = [];
        data.EyeChan_label = {'Eye_Hor' 'Eye_Ver'};

        %% Sanity Check
        for i=1:numtr
            data.trial{i}(isnan(data.trial{i}(:,:)))= 0;
        end
    end

    %% Remove mean
    tmp_cfg = [];
    tmp_cfg.demean = 'yes';
    data = ft_preprocessing(tmp_cfg,data);


    %%  Bipolar Derivation
    
    if bipolar

        disp('Rereferencing the dataset')
        if strcmp(cfg.bipolarder,'average')
            tmp_cfg = [];
            tmp_cfg.reref = 'yes';
            tmp_cfg.refchannel = cfg.chan_list;%'all';
        elseif strcmp(cfg.bipolarder,'ownmontage')
            tmp_cfg = [];
            tmp_cfg.reref = 'yes';
            tmp_cfg.montage = cls_createmontage(data.label,cfg.bipolarder_list);
            tmp_cfg.channel = data.label;
        else
            tmp_cfg = [];
            tmp_cfg=load(cfg.bipolarder);
    %         tmp_cfg.reref = 'yes';
        end

        data_bipolar = ft_preprocessing(tmp_cfg,data);

        if strcmp(cfg.bipolarder,'average') && ~any(strcmp(cfg.chan_list,'all'))
            sel1 = match_str(data_bipolar.label,cfg.chan_list);
            data_bipolar.label = data_bipolar.label(sel1);
            numtr_ses = size(data_bipolar.trial,2);
            for i=1:numtr_ses
                fprintf('Trial %d\n',i)
                data_bipolar.trial{i} = data_bipolar.trial{i}(sel1,:);
            end
        end
        % numtr = size(data_bipolar,2);
        % 
        % for i=1:numtr
        %     tmp = nansum(data_bipolar.trial{i},2);
        %     index = find(tmp==0);
        %     data_bipolar.trial{i}(index,:) = [];
        % end
        % 
        % data_bipolar.label(index) = [];
        %data_bipolar.EyeChan = data.EyeChan;
        %data_bipolar.EyeChan_label = data.EyeChan_label;

        %% Remove channels with the sum of all trials is 0
        % a_indx = [];
        % for tt=1:length(data_bipolar.trial)
        %     tmpval= sum(data_bipolar.trial{tt},2);
        %     if any(tmpval==0)
        %         a_indx{tt}= find(tmpval==0);
        %     end
        %     
        % end
        % 
        % if 0
        % if ~isempty(a_indx)
        %     for tr=1:length(data_bipolar.trial)
        %         try
        %         if ~isempty(a_indx{tr})
        %             sprintf('Removing channel numbers %d',a_indx{tr}(:));
        %             data_bipolar.trial{tr}(a_indx{tr},:) = [];
        %         end
        %         data_bipolar.label{tr} = [];
        %         catch
        %         end
        %     end
        % end
        % end
    else
        % do the channel selection here
        tmp_cfg = [];
        load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\kurt_montage_layout\kurt_COI_alongfiner2.mat');
        tmp_cfg.channel = COI;
        data_bipolar = ft_preprocessing(tmp_cfg,data);
    end
    
    %%
    if strcmp(cfg.norm,'yes')
        s_cfg = [];
        s_cfg.channelnormalise = 'yes';
        data = ft_preprocessing(s_cfg,data_bipolar);
        %data = session_norm2(data_bipolar);
    else
        data = data_bipolar;
    end

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
else
    %% baseline calculation
    %we need to known how many trials per session; so we use mcfg.numtr to
    %that purpose
    
    dataAttIn = data;
    dataAttOut = [];
    
end
disp('Ready')
