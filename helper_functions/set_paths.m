function set_paths(des_data)
    addpath('C:\Users\Jordan\Documents\cls_thesis\matlab\fieldtrip-20210411')
    addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\cls_thesis\helper_functions\')
    addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\kurt_montage_layout\')
    addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\pele_montage_layout\')

    if strcmp(des_data, 'baseline')
        addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\raw_k_baseline_data\')
        addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\raw_p_baseline_data\')
    else
        addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\raw_k_postbaseline_data\')
        addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\raw_p_postbaseline_data\')
    end