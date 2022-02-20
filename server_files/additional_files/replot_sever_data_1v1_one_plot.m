%% Set up paths
addpath('/home/jordan/neuro_thesis/neuro_thesis/cls_thesis/gc_hierarchies/')
addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
%addpath('D:/MATLAB/mvgc_v1.0')
startup
addpath('/home/jordan/common/matlab/fieldtrip-20210411')
ft_defaults

format short;
clear all;
close all;clc;
disp('starting')

monkey = 'kurt';
model_free = 1;

if model_free
    if strcmp(monkey, 'kurt')
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt/kurt_model_free_gc_one_v_one.mat')
    else
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele/pele_model_free_gc_one_v_one.mat')
    end
else
    if strcmp(monkey, 'kurt')
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt/kurt_gc_one_v_one.mat')
    else
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele/pele_gc_one_v_one.mat')
    end
end

disp("done loading")

figure(1);

compare_counter = 1;
plot_counter = 1;

for r1 = 1 : length(regions)-1
    for r2 = 1 : length(regions)
        if r1 == r2
            compare_counter = compare_counter+1;
            continue
        else
%             if strcmp(region_names(r1), 'TEO') && strcmp(region_names(r2), 'V4')
                gc_1_ave = mean(gc_forward{compare_counter}, 1);
                gc_2_ave = mean(gc_backward{compare_counter}, 1);

                %% do some plotting
                x_range = (1:1:length(gc_1_ave))./length(gc_1_ave);
                x_range = x_range.*fnq;
                blue=[.3 .8 .9];red=[1 .1 .1];

                plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
                plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on
                legend('FF','FB');legend boxoff;
                xlim([0 140]);
%                 ylim([0 real(zgc)]);
                title(sprintf("%s v.s. %s", region_names(r1), region_names(r2)));
                set(gca,'Layer','top');
                saveas(gcf,sprintf('individual_plots/%s_vs_%s.png',region_names(r1),region_names(r2)));
                saveas(gcf,sprintf('individual_plots/%s_vs_%s.fig',region_names(r1),region_names(r2)));
                close(gcf)
                compare_counter = compare_counter+1;
%             end
        end 
    end
end

han = axes(figure(1), 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han,'Granger Causality');xlabel(han, 'Frequency (Hz)');
title(han,sprintf('%s Granger Causeality of Regions Against all Other \n', monkey));
