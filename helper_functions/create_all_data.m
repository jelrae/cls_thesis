clear all
clc

% monkeys = {'pele' 'kurt'};
% des_datas = {'baseline', 'post'}; 
monkeys = {'pele'};
des_datas = {'post'}; 
% des_datas = {'post'};
data_used = 'all';
bipolar = true;

for i=1:length(monkeys)
    for j=1:length(des_datas)
        disp(monkeys{i})
        disp(des_datas{j})
        if bipolar
            create_bipolarderivatives_function(monkeys{i}, des_datas{j}, data_used);
        else
            create_no_bipolar_data_function(monkeys{i}, des_datas{j}, data_used, bipolar);
        end
    end
end