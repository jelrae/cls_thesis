addpath('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\cls_thesis\helper_functions\')
set_paths_roi;

clear all
clc
ft_defaults

lineColor = linspecer(2,'qualitative');

%% Select the monkey and define region interest (hyperparamters)

monkeys = {'kurt' 'pele'};
for m = monkeys
    monkey = m{1};
    % Run script to get ROIs specific to the monkey
    fig_6_ROIS;
    if strcmp(m{1}, 'kurt')
        monkey_caps = 'Kurt';
        % Bad channels giving power in 100's
        a7A(6) = [];
        a7A(5) = [];
    else
        monkey_caps = 'Pele';
    end
    load(sprintf('%s_p_all_AttIn.mat',m{1}));
    load(sprintf('%s_p_all_AttOut.mat',m{1}));
    
    for r1 = 1:length(region_names)-1
        region1_string = region_names{r1};
        eval(sprintf('roi = %s;', region1_string));
        for r2 = r1+1:length(region_names)
            region2_string = region_names{r2};
            eval(sprintf('roi_2 = %s;', region2_string));
            roi_combine = [roi roi_2];
            roi_list_in = cell(0,2);
            roi_list_out = cell(0,2);
            
            for i = roi
                if ismember(i, all_AttIn.label)
                    for j = roi_2
                        if ismember(j, all_AttIn.label)
                            roi_list_in(end+1,:) = {i{1} j{1}};
                        end
                    end
                end
                
                if ismember(i, all_AttOut.label)
                    for j = roi_2
                        if ismember(j, all_AttOut.label)
                            roi_list_out(end+1,:) = {i{1} j{1}};
                        end
                    end
                end
            end

            %% Calculate Fourier on ALL channels Attention In

            pow_cfg = [];
            pow_cfg.average = 'yes';
            pow_cfg.method = 'mtmfft';
            pow_cfg.keeptapers = 'yes';
            pow_cfg.output = 'fourier';
            pow_cfg.channel = roi_combine;

            pow_cfg.keeptrials = 'yes';%'no';
            pow_cfg.pad='maxperlen'; % Padding: not adding zeros
            pow_cfg.flag = 0;

            % for the beta and gamma range
            pow_cfg.tapsmofrq = 4;
            pow_cfg.foilim = [2 140];
            pow_cfg.taper = 'dpss';

            b_g_coh_analysis = ft_freqanalysis(pow_cfg,all_AttIn);

            % for the theta  range
            pow_cfg.tapsmofrq = 1;
            pow_cfg.foilim = [2 12];
            pow_cfg.taper = 'hanning';

            t_coh_analysis = ft_freqanalysis(pow_cfg,all_AttIn);

            % Connectivity (coherence) analysis 
            cfg = [];
            cfg.method = 'coh';
            cfg.channelcmb = roi_list_in;
            b_g_coh_analysis = ft_connectivityanalysis(cfg, b_g_coh_analysis);
            t_coh_analysis = ft_connectivityanalysis(cfg, t_coh_analysis);

            %Start peak align and plotting

            fig1 = figure(1);

            % Calculate the peak theta and align

            peak_start = 2;
            peak_end = 5;
            w = 2;

            p_info = PeakAlignSpectrum(peak_start,peak_end,w,t_coh_analysis.cohspctrm);

            g_min = p_info.ave_loc - w;
            g_max = p_info.ave_loc + w;
            sp1 = subplot(1,3,1);
            plot(t_coh_analysis.freq(g_min:g_max), p_info.peak, ...
                'DisplayName', 'AttIn', 'Color', lineColor(1,:), 'LineWidth', 1.25);
            xlabel('Theta')

            max_t1 = max(p_info.peak);
            max_t1_x = g_max;
            min_t1_x = g_min;

            hold on

            % Calculate the peak beta and align

            peak_start = 5;
            peak_end = 15;
            w = 5;

            p_info = PeakAlignSpectrum(peak_start,peak_end,w,b_g_coh_analysis.cohspctrm);

            g_min = p_info.ave_loc - w;
            g_max = p_info.ave_loc + w;

            sp2 = subplot(1,3,2);
            plot(b_g_coh_analysis.freq(g_min:g_max), p_info.peak, ...
                'DisplayName', 'AttIn', 'Color', lineColor(1,:), 'LineWidth', 1.25);
            xlabel('Beta')

            max_b1 = max(p_info.peak);
            max_b1_x = g_max;
            min_b1_x = g_min;

            hold on

            % Calculate the peak gamma and align

            peak_start = 30;
            peak_end = 46;
            w = 10;

            p_info = PeakAlignSpectrum(peak_start,peak_end,w,b_g_coh_analysis.cohspctrm);

            g_min = p_info.ave_loc - w;
            g_max = p_info.ave_loc + w;

            sp3 = subplot(1,3,3);
            plot(b_g_coh_analysis.freq(g_min:g_max), p_info.peak, ...
                'DisplayName', 'AttIn', 'Color', lineColor(1,:), 'LineWidth', 1.25);
            xlabel('Gamma')

            max_g1 = max(p_info.peak);
            max_g1_x = g_max;
            min_g1_x = g_min;

            hold on

            %% Calculate Fourier on ALL channels Attention Out

            pow_cfg = [];
            pow_cfg.average = 'yes';
            pow_cfg.method = 'mtmfft';
            pow_cfg.keeptapers = 'yes';
            pow_cfg.output = 'fourier';
            pow_cfg.channel = roi_combine;

            pow_cfg.keeptrials = 'yes';%'no';
            pow_cfg.pad='maxperlen'; % Padding: not adding zeros
            pow_cfg.flag = 0;

            % for the beta and gamma range
            pow_cfg.tapsmofrq = 4;
            pow_cfg.foilim = [2 140];
            pow_cfg.taper = 'dpss';

            b_g_coh_analysis = ft_freqanalysis(pow_cfg,all_AttOut);

            % for the theta  range
            pow_cfg.tapsmofrq = 1;
            pow_cfg.foilim = [2 12];
            pow_cfg.taper = 'hanning';

            t_coh_analysis = ft_freqanalysis(pow_cfg,all_AttOut);

            % Connectivity (coherence) analysis 
            cfg = [];
            cfg.method = 'coh';
            cfg.channelcmb = roi_list_out;
            b_g_coh_analysis = ft_connectivityanalysis(cfg, b_g_coh_analysis);
            t_coh_analysis = ft_connectivityanalysis(cfg, t_coh_analysis);

            % Start peak align and plotting
            % Calculate the peak theta and align

            peak_start = 2;
            peak_end = 5;
            w = 2;

            p_info = PeakAlignSpectrum(peak_start,peak_end,w,t_coh_analysis.cohspctrm);

            g_min = p_info.ave_loc - w;
            g_max = p_info.ave_loc + w;
            sp1 = subplot(1,3,1);
            plot(t_coh_analysis.freq(g_min:g_max), p_info.peak, ...
                'DisplayName', 'AttOut', 'Color', lineColor(2,:), 'LineWidth', 1.25);
            xlabel('Theta')

            max_t2 = max(p_info.peak);
            max_t2_x = g_max;
            min_t2_x = g_min;

            legend

            hold on

            % Calculate the peak beta and align

            peak_start = 5;
            peak_end = 15;
            w = 5;

            p_info = PeakAlignSpectrum(peak_start,peak_end,w,b_g_coh_analysis.cohspctrm);

            g_min = p_info.ave_loc - w;
            g_max = p_info.ave_loc + w;

            sp2 = subplot(1,3,2);
            plot(b_g_coh_analysis.freq(g_min:g_max), p_info.peak, ...
                'DisplayName', 'AttOut', 'Color', lineColor(2,:), 'LineWidth', 1.25);
            xlabel('Beta')

            max_b2 = max(p_info.peak);
            max_b2_x = g_max;
            min_b2_x = g_min;

            hold on

            % Calculate the peak gamma and align

            peak_start = 30;
            peak_end = 46;
            w = 10;

            p_info = PeakAlignSpectrum(peak_start,peak_end,w,b_g_coh_analysis.cohspctrm);

            g_min = p_info.ave_loc - w;
            g_max = p_info.ave_loc + w;

            sp3 = subplot(1,3,3);
            plot(b_g_coh_analysis.freq(g_min:g_max), p_info.peak, ...
                'DisplayName', 'AttOut', 'Color', lineColor(2,:), 'LineWidth', 1.25);
            xlabel('Gamma')

            max_g2 = max(p_info.peak);
            max_g2_x = g_max;
            min_g2_x = g_min;

            hold on


            %% Formatting of plots

            max_y = max([max_t1, max_b1,max_g1,max_t2, max_b2,max_g2]);

            ylim(sp1,[0 max_y+.02]);
            ylim(sp2,[0 max_y+.02]);
            ylim(sp3,[0 max_y+.02]);

            xlim(sp1,[t_coh_analysis.freq(min([min_t1_x min_t2_x])) t_coh_analysis.freq(max([max_t1_x max_t2_x]))]);
            xlim(sp2,[b_g_coh_analysis.freq(min([min_b1_x min_b2_x])) b_g_coh_analysis.freq(max([max_b1_x max_b2_x]))]);
            xlim(sp3,[b_g_coh_analysis.freq(min([min_g1_x min_g2_x])) b_g_coh_analysis.freq(max([max_g1_x max_g2_x]))]);

            sfh1 = subplot(1,3,1, 'Parent', gcf);
            sfh2 = subplot(1,3,2, 'Parent', gcf);
            sfh3 = subplot(1,3,3, 'Parent', gcf);

            sfh1.Position = [0.11 0.13 0.15 0.805];
            sfh2.Position = [0.33 0.13 0.25 0.805];
            sfh3.Position = [0.655 0.13 0.31 0.805];

            % Plotting labels
            han = axes(figure(1), 'visible', 'off');
            han.Title.Visible = 'on';
            han.XLabel.Visible = 'on';
            han.YLabel.Visible = 'on';
            ylh = ylabel(han,'Coherence');xlh = xlabel(han, 'Frequency (Hz)');
            ylh.Position(1) = -.1;
            xlh.Position(2) = -.07;
            formatSpec = 'Coherence Analysis %s %s:%s';
            title(han, sprintf(formatSpec, monkey, region1_string,region2_string))

            saveas(fig1, sprintf('results/tbg_%s_%s_%s.jpg', monkey_caps, region1_string, region2_string))
            close(fig1)
            
        end
    end
    
end
