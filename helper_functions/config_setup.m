function cfg = config_setup(roi, monkey, des_data, data_used)
    %% Configuration parameters - the following are consistent
    cfg = [];
    cfg.equal = 'yes';
    
    lpf = 1;
    
    if lpf
        cfg.lpfilter = 1;
    end
    
    if strcmp(des_data, 'baseline')
        cfg.bsln = 1;
    else
        cfg.bsln = 0;
    end
    cfg.norm = 'yes';
    cfg.montage = 'yes'; %here you need to add the name of the montage structure you want to use... alongfinger_v2
%     cfg.channel = roi; % This is one difference

    %% Assigning parameter 'monkey' above sets up the cfg with monkey specific info
    if strcmp(monkey,'kurt')
        % Monkey specific information
        monkey_code = 'K';
        cfg.monkey = 'ku';
        if strcmp(des_data, 'baseline')
            cfg.path = 'C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\raw_k_baseline_data\';
            if strcmp(data_used, 'all')
                cfg.v_session = [67 68 69];
            else
                cfg.v_session = [67 68];
            end
        else
            cfg.path = 'C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\raw_k_postbaseline_data\';
            if strcmp(data_used, 'all')
                cfg.v_session = [40 41 42 43 44 45 46 48]; % Vector referring to number of sessions
            else
                cfg.v_session = [40 41];
            end
        end
        cfg.bipolarder = 'C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\kurt_montage_layout\kurt_montage_bipolar_alongfinger_v2.mat';
    elseif strcmp(monkey,'pele')    
        monkey_code = 'P';
        cfg.monkey = 'pe';
        if strcmp(des_data, 'baseline')
            cfg.path = 'C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\raw_p_baseline_data\';
            if strcmp(data_used, 'all')
                cfg.v_session = [39 40 42 87 88 90]; % Vector referring to number of sessions
            else
                cfg.v_session = [39 40]; % Vector referring to number of sessions old ones [39 40]
            end
        else
            cfg.path = 'C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\raw_p_postbaseline_data\';
            if strcmp(data_used, 'all')
                cfg.v_session = [28 29 30 31 32 33 34 35 36 37 38 39 40 42]; % Vector referring to number of sessions
            else 
                cfg.v_session = [28 29]; % Vector referring to number of sessions [28 29]
            end
        end
        cfg.bipolarder = 'C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\pele_montage_layout\pele_montage_bipolar_alongfinger_v2.mat';
    end
    disp('---------------Config made---------------')