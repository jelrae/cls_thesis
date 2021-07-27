%% Adapted from MVGC toolbox: 
% https://github.com/ajaykarpur/granger/blob/master/lib/mvgc/utils/plot_spw.m
function gc_plot(f, fs)


h = size(f,2);
fres = h-1;
lam = sfreqs(fres,fs)';


xlims = [lam(1) lam(end)];
ylims = [min(f(:)) 1.1*max(f(:))];
plot(lam,f(1,:), 'Color',[0.3 0.7 0.3],'LineWidth',2)
hold on
plot(lam,f(2,:), 'Color',[0.99 0.45 0.1],'LineWidth',2)
set(gca,'box','off');
legend('V1 to V4','V4 to V1');
xlim(xlims);
ylim(ylims);
xlabel('Frequency (Hz)')
ylabel('Granger Causality');