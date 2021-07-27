clear all
clc

monkeys = {'pele' 'kurt'};
% des_datas = {'baseline', 'post'}; 
des_datas = {'post'};
data_used = 'all';

for i=1:length(monkeys)
    for j=1:length(des_datas)
        disp(monkeys{i})
        disp(des_datas{j})
        monkey = monkeys{i};
        des_data = des_datas{j};
        exp_desc = 'Coherence';
        roi_string = 'V1 and V4'; % Set for plotting


        %% Add paths
        set_paths(des_data);

        ft_defaults

        %% Run script to get ROIs specific to the monkey
        vbm_rois;
        % Define region of interest 1
        roi = V1;
        %Define region of interest 2
        roi_2 = V4;
        roi_combine = [roi roi_2];

        %% Configuration parameters and bipolar merging
        cfg = config_setup(roi, monkey, des_data, data_used);
        mcfg = cfg;
    end
end