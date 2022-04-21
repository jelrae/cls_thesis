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
model_free = 0;

% if model_free
%     if strcmp(monkey, 'kurt')
%         load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt/kurt_model_free_gc_one_v_all.mat')
%     else
%         load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele/pele_model_free_gc_one_v_all.mat')
%     end
% else
%     if strcmp(monkey, 'kurt')
%         load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt/kurt_gc_one_v_all.mat')
%     else
%         load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele/pele_gc_one_v_all.mat')
%     end
% end

load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/kurt_one_v_all_mf_short_data')

load('/home/jordan/neuro_thesis/cls_thesis/server_files/results/pele_one_v_all_mf_short_data')

disp("done loading")

fres = 250;
fs = 1000;

x_range = sfreqs(fres,fs);
green=[0 1 0];black=[0 0 0];

figure(1);

for region = 1 : length(k_gc_forward)
    %% average over monkey and then over monkeys
%     k_gc_1_ave = mean(k_gc_forward{region}, 1);
%     k_gc_2_ave = mean(k_gc_backward{region}, 1);
%     p_gc_1_ave = mean(p_gc_forward{region}, 1);
%     p_gc_2_ave = mean(p_gc_backward{region}, 1);
%     gc_1_ave = (k_gc_1_ave + p_gc_1_ave) / 2;
%     gc_2_ave = (k_gc_2_ave + p_gc_2_ave) / 2;
    gc_one = cat(1,k_gc_forward{region}, p_gc_forward{region});
    gc_two = cat(1,k_gc_backward{region}, p_gc_backward{region});
    
    %% do some plotting
    % Calculate the peak theta and align
    peak_start = 2;
    peak_end = 5;
    w = 2;
    
    p_info_backward = PeakAlignSpectrum(peak_start,peak_end,w,gc_one);
    p_info_forward = PeakAlignSpectrum(peak_start,peak_end,w,gc_two);
    g_min_f = p_info_forward.ave_loc - w;
    g_max_f = p_info_forward.ave_loc + w;
    g_min_b = p_info_backward.ave_loc - w;
    g_max_b = p_info_backward.ave_loc + w;
    
%     sp1 = subplot(1,3,1);    
    sp1 = subplot(length(p_region_names), 3, (3*region) - 2);
    
    
    % get the std lines
    lower_f_std = p_info_forward.peak - (2*p_info_forward.std);
    upper_f_std = p_info_forward.peak + (2*p_info_forward.std);
    lower_b_std = p_info_backward.peak - (2*p_info_backward.std);
    upper_b_std = p_info_backward.peak + (2*p_info_backward.std);
    in_between_f = [upper_f_std, fliplr(lower_f_std)];
    in_between_b = [upper_b_std, fliplr(lower_b_std)];
    
    wrap_x_range_f = [transpose(x_range(g_min_f:g_max_f)), fliplr(transpose(x_range(g_min_f:g_max_f)))];
    wrap_x_range_b = [transpose(x_range(g_min_b:g_max_b)), fliplr(transpose(x_range(g_min_b:g_max_b)))];
    
    %Plot the std
    hold on
    fill(wrap_x_range_f, in_between_f, green, 'FaceAlpha',0.4, 'LineStyle','none');
    fill(wrap_x_range_b, in_between_b, black, 'FaceAlpha',0.5, 'LineStyle','none');
    
    plot(x_range(g_min_f:g_max_f), p_info_forward.peak, ...
        'DisplayName', 'FF', 'Color', green, 'LineWidth', 1.25);
    hold on
    plot(x_range(g_min_b:g_max_b), p_info_backward.peak, ...
        'DisplayName', 'FB', 'Color', black, 'LineWidth', 1.25);
    xlabel('Theta')
    title(sprintf("%s v.s. All", p_region_names(region)));
    set(gca,'Layer','top');
    max_t_f = max(upper_f_std);
    min_t_f = min(lower_f_std);
    max_t_x_f = g_max_f;
    min_t_x_f = g_min_f;
    max_t_b = max(upper_b_std);
    min_t_b = min(lower_b_std);
    max_t_x_b = g_max_b;
    min_t_x_b = g_min_b;
    hold on

%     plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
%     plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%     legend('FF','FB');legend boxoff;
%     xlim([0 10]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%     ylim([0 zgc]);
%     title(sprintf("%s v.s. all Theta", p_region_names(region)));
%     set(gca,'Layer','top');
    
    % Beta
    % Calculate the peak beta and align
    peak_start = 5;
    peak_end = 15;
    w = 5;
    p_info_backward = PeakAlignSpectrum(peak_start,peak_end,w,gc_one);
    p_info_forward = PeakAlignSpectrum(peak_start,peak_end,w,gc_two);
    g_min_f = p_info_forward.ave_loc - w;
    g_max_f = p_info_forward.ave_loc + w;
    g_min_b = p_info_backward.ave_loc - w;
    g_max_b = p_info_backward.ave_loc + w;
%     sp2 = subplot(1,3,2);

    sp2 = subplot(length(p_region_names), 3, (3*region) - 1);
    
    % get the std lines
    lower_f_std = p_info_forward.peak - (2*p_info_forward.std);
    upper_f_std = p_info_forward.peak + (2*p_info_forward.std);
    lower_b_std = p_info_backward.peak - (2*p_info_backward.std);
    upper_b_std = p_info_backward.peak + (2*p_info_backward.std);
    in_between_f = [upper_f_std, fliplr(lower_f_std)];
    in_between_b = [upper_b_std, fliplr(lower_b_std)];
    
    wrap_x_range_f = [transpose(x_range(g_min_f:g_max_f)), fliplr(transpose(x_range(g_min_f:g_max_f)))];
    wrap_x_range_b = [transpose(x_range(g_min_b:g_max_b)), fliplr(transpose(x_range(g_min_b:g_max_b)))];
    
    %Plot the std
    hold on
    fill(wrap_x_range_f, in_between_f, green, 'FaceAlpha',0.4, 'LineStyle','none');
    fill(wrap_x_range_b, in_between_b, black, 'FaceAlpha',0.5, 'LineStyle','none');
    plot(x_range(g_min_f:g_max_f), p_info_forward.peak, ...
        'DisplayName', 'FF', 'Color', green, 'LineWidth', 1.25);
    hold on
    plot(x_range(g_min_b:g_max_b), p_info_backward.peak, ...
        'DisplayName', 'FB', 'Color', black, 'LineWidth', 1.25);
    xlabel('Beta')
    title(sprintf("%s v.s. All", p_region_names(region)));
    set(gca,'Layer','top');
%     legend
    max_b_f = max(in_between_f);
    min_b_f = min(in_between_f);
    max_b_x_f = g_max_f;
    min_b_x_f = g_min_f;
    max_b_b = max(in_between_b);
    min_b_b = min(in_between_b);
    max_b_x_b = g_max_b;
    min_b_x_b = g_min_b;
    hold on


%     subplot(length(p_region_names), 3, (3*region) - 1)
%     plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
%     plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%     legend('FF','FB');legend boxoff;
%     xlim([10 30]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%     ylim([0 zgc]);
%     title(sprintf("%s v.s. all Beta", p_region_names(region)));
%     set(gca,'Layer','top');
    
    % Gamma
    % Calculate the peak gamma and align
    peak_start = 27;
    peak_end = 46;
    w = 10;
    p_info_backward = PeakAlignSpectrum(peak_start,peak_end,w,gc_one);
    p_info_forward = PeakAlignSpectrum(peak_start,peak_end,w,gc_two);
    g_min_f = p_info_forward.ave_loc - w;
    g_max_f = p_info_forward.ave_loc + w;
    g_min_b = p_info_backward.ave_loc - w;
    g_max_b = p_info_backward.ave_loc + w;
    sp3 = subplot(length(p_region_names), 3, (3*region) - 0);
%     sp3 = subplot(1,3,3);

    % get the std lines
    lower_f_std = p_info_forward.peak - (2*p_info_forward.std);
    upper_f_std = p_info_forward.peak + (2*p_info_forward.std);
    lower_b_std = p_info_backward.peak - (2*p_info_backward.std);
    upper_b_std = p_info_backward.peak + (2*p_info_backward.std);
    in_between_f = [upper_f_std, fliplr(lower_f_std)];
    in_between_b = [upper_b_std, fliplr(lower_b_std)];
    
    wrap_x_range_f = [transpose(x_range(g_min_f:g_max_f)), fliplr(transpose(x_range(g_min_f:g_max_f)))];
    wrap_x_range_b = [transpose(x_range(g_min_b:g_max_b)), fliplr(transpose(x_range(g_min_b:g_max_b)))];
    
    %Plot the std
    hold on
    fill(wrap_x_range_f, in_between_f, green, 'FaceAlpha',0.4, 'LineStyle','none');
    fill(wrap_x_range_b, in_between_b, black, 'FaceAlpha',0.5, 'LineStyle','none');
 
    plot(x_range(g_min_f:g_max_f), p_info_forward.peak, ...
        'DisplayName', 'FF', 'Color', green, 'LineWidth', 1.25);
    hold on
    plot(x_range(g_min_b:g_max_b), p_info_backward.peak, ...
        'DisplayName', 'FB', 'Color', black, 'LineWidth', 1.25);
    xlabel('Gamma')
    title(sprintf("%s v.s. All", p_region_names(region)));
    set(gca,'Layer','top');
    max_g_f = max(in_between_f);
    min_g_f = min(in_between_f);
    max_g_x_f = g_max_f;
    min_g_x_f = g_min_f;
    max_g_b = max(in_between_b);
    min_g_b = min(in_between_b);
    max_g_x_b = g_max_b;
    min_g_x_b = g_min_b;
    hold on
    
    %% Formatting of plots
    max_y = max([max_t_f, max_b_f,max_g_f,max_t_b, max_b_b,max_g_b]);
    min_y = min([min_t_f, min_b_f,min_g_f,min_t_b, min_b_b,min_g_b]);
    ylim(sp1,[min_y-(0.075*min_y) max_y+(0.075*max_y)]);
    ylim(sp2,[min_y-(0.075*min_y) max_y+(0.075*max_y)]);
    ylim(sp3,[min_y-(0.075*min_y) max_y+(0.075*max_y)]);
    xlim(sp1,[x_range(min([min_t_x_f min_t_x_b])) x_range(max([max_t_x_f max_t_x_b]))]);
    xlim(sp2,[x_range(min([min_b_x_f min_b_x_b])) x_range(max([max_b_x_f max_b_x_b]))]);
    xlim(sp3,[x_range(min([min_g_x_f min_g_x_b])) x_range(max([max_g_x_f max_g_x_b]))]);
    sfh1 = subplot(length(p_region_names), 3, (3*region) - 2, 'Parent', gcf);
    sfh2 = subplot(length(p_region_names), 3, (3*region) - 1, 'Parent', gcf);
    sfh3 = subplot(length(p_region_names), 3, (3*region) - 0, 'Parent', gcf);
%     sfh1.Position = [0.11 0.13 0.15 0.805];
%     sfh2.Position = [0.33 0.13 0.25 0.805];
%     sfh3.Position = [0.655 0.13 0.31 0.805];

    % Plotting labels
%     han = axes(figure(1), 'visible', 'off');
%     han.Title.Visible = 'on';
%     han.XLabel.Visible = 'on';
%     han.YLabel.Visible = 'on';
%     ylh = ylabel(han,'Granger-Causal Influence');xlh = xlabel(han, 'Frequency (Hz)');
%     ylh.Position(1) = -.1;
%     xlh.Position(2) = -.07;
%     formatSpec = 'Granger Analysis %s %s:%s';
%     title(han, sprintf(formatSpec, monkey_caps, region1_string,region2_string))

%     subplot(length(p_region_names), 3, (3*region) - 0)
%     plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
%     plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on;
%     legend('FF','FB');legend boxoff;
%     xlim([30 100]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%     ylim([0 zgc]);
%     title(sprintf("%s v.s. all Gamma", p_region_names(region)));
%     set(gca,'Layer','top');

%     subplot(length(p_region_names), 4, (region*4))
%     plot(x_range,gc_2_ave,'Color',blue,'LineWidth',3);hold on;
%     plot(x_range,gc_1_ave,'Color',red,'LineWidth',3);hold on
%     legend('FF','FB');legend boxoff;
%     xlim([0 140]);zgc=1.1*max(max(max(gc_2_ave),max(gc_1_ave)));
%     ylim([0 zgc]);
%     title(sprintf("%s v.s. all Full", p_region_names(region)));
%     set(gca,'Layer','top');
    
end

subplot(length(p_region_names), 3, 3);
legend(['Feed Forward'; 'Feed Back'])


han = axes(figure(1), 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han,'Granger Causal Influence');xlabel(han, 'Frequency (Hz)');
title(han,sprintf('Granger Causal Influence: Each Regions v,s All Other \n'));
