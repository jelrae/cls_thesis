function power_check(data, cnd)
    plot_loglog = 'yes';
    %% Fourier transfrom
%     tapsmofrq = [8 4 2]; %last with 12 instead 0f 9
    tapsmofrq = [8]; %last with 12 instead 0f 9
    foilim = {[20 140],[8 52],[2 15]}; %this is old{[20 140] [5 30] [2 10]};
    pow_cfg = [];
    pow_cfg.average = 'yes';
    pow_cfg.method = 'mtmfft';
    pow_cfg.keeptapers = 'yes';
    pow_cfg.output = 'fourier';

    pow_cfg.keeptrials = 'yes';%'no';
    pow_cfg.pad='maxperlen'; % Padding: not adding zeros
    pow_cfg.flag = 0;
    
    figure(1);
    
    for z=1:length(tapsmofrq)
        pow_cfg.tapsmofrq = tapsmofrq(z);
        pow_cfg.foilim = foilim{z};
        if z==3
            pow_cfg.taper = 'hanning';
        else
            pow_cfg.taper = 'dpss';
        end

        for n = 1:length(cnd)
            pow_cfg.channel = cnd{n};
            disp(pow_cfg.channel);
            tmp_angle = [];
            tmp = [];
            tmp_pow = [];
            sprintf('Condition %s',cnd{n})
            tmp = ft_freqanalysis(pow_cfg,data);
            
            %% Freqdescriptives
            pow2_cfg = [];
            pow2_cfg.psi = 'no'; % Phase slope Index
            pow2_cfg.channel = cnd{n};
            pow2_cfg.jackknife = 'yes';
            pow2_cfg.avgChOI = 'yes';

            %% Obtaining and plotting the power analysis

            % Frequency of interest limits
            tmp.freq;

            % Power
            tmp = ft_freqdescriptives(pow2_cfg, tmp);
            if strcmp(plot_loglog,'yes')
                loglog(tmp.freq,tmp.powspctrm, ...
                    'DisplayName', cnd{n})
            elseif strcmp(plot_loglog,'no')
                plot(tmp.freq,tmp.powspctrm, ...
                    'DisplayName', cnd{n})
            end
            hold on
%             eval(sprintf('powertap%d_%s=ft_freqdescriptives(pow2_cfg, tmp);',tapsmofrq(z),cnd{n}));
        end
    end

    
    %% Plotting
    % Plotting labels
    ylabel('Power')
    xlabel('Frequency (Hz)')
    formatSpec = 'Power analysis for Monkey %s: %s';
    title("Power check for individual regions")
    legend
end