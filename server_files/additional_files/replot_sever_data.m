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

monkey = 'pele';
model_free = 1;

if model_free
    if strcmp(monkey, 'kurt')
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt/kurt_model_free_gc_one_v_all.mat')
    else
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele/pele_model_free_gc_one_v_all.mat')
    end
else
    if strcmp(monkey, 'kurt')
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt/kurt_gc_one_v_all.mat')
    else
        load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele/pele_gc_one_v_all.mat')
    end
end

disp("done loading")

figure(1);

for region = 1 : length(gc_forward)
    gc_1_ave = mean(gc_forward{region}, 1);
    gc_2_ave = mean(gc_backward{region}, 1);
    %% do some plotting
    x_range = (1:1:length(gc_1_ave))./length(gc_1_ave);
    x_range = x_range.*fnq;
    blue=[.3 .8 .9];red=[1 .1 .1];

    subplot(length(regions), 4, (4*region) - 3)
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
    legend('FF','FB');legend boxoff;
    xlim([0 10]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);
    title(sprintf("%s v.s. all Theta", region_names(region)));
    set(gca,'Layer','top');

    subplot(length(regions), 4, (4*region) - 2)
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
    legend('FF','FB');legend boxoff;
    xlim([10 30]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);
    title(sprintf("%s v.s. all Beta", region_names(region)));
    set(gca,'Layer','top');

    subplot(length(regions), 4, (4*region) - 1)
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
    legend('FF','FB');legend boxoff;
    xlim([30 100]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);
    title(sprintf("%s v.s. all Gamma", region_names(region)));
    set(gca,'Layer','top');

    subplot(length(regions), 4, (region*4))
    plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
    plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on
    legend('FF','FB');legend boxoff;
    xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
    ylim([0 zgc]);
    title(sprintf("%s v.s. all Full", region_names(region)));
    set(gca,'Layer','top');
    
end

han = axes(figure(1), 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han,'Granger Causality');xlabel(han, 'Frequency (Hz)');
title(han,sprintf('%s Granger Causeality of Regions Against all Other \n', monkey));
