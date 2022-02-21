%% simulate RSA matrix for linear coding of value (or, why does OFC RSA have the structure that it does?)

% see supplementary note of paper for discussion

nUnits = 1000; %number of neurons
baseline_firing_rate = poissrnd(4,nUnits,1); %for each neuron
linear_value_slope = randn(nUnits,1); %how much each neuron correlates with value
value_levels = [-2:2]; %different possible levels of expected value (demeaned)
fr_noise = 10; %noise level

%model firing rates (fr)
for i = 1:nUnits %loop over units
    for v = 1:5 %loop over different levels of value
        fr(i,v)   = baseline_firing_rate(i) + linear_value_slope(i)*value_levels(v) + poissrnd(fr_noise);
        fr(i,v+5) = baseline_firing_rate(i) + linear_value_slope(i)*value_levels(v) + poissrnd(fr_noise);
        fr(i,v+10) = baseline_firing_rate(i) + linear_value_slope(i)*value_levels(v) + poissrnd(fr_noise);
        fr(i,v+15) = baseline_firing_rate(i) + linear_value_slope(i)*value_levels(v) + poissrnd(fr_noise);
    end
end
fr = fr - repmat(mean(fr,2),[1 20]); %demean firing rate 
imagesc(corrcoef(fr)); %calculate RSA matrix across simulated firing rates
colormap('hot');
caxis([-0.3 0.3]);
set(gca,'YTickLabel',[1:5],'YTick',1:20,'XTickLabel',[1:5],'XTick',1:20,'FontSize',10);
xlabel(sprintf(['Prob.              Mag.              Prob.             Mag.\n' ...
    'Left                                      Right']))
set(get(gca,'XLabel'),'FontSize',14);
ylabel(sprintf(['Right                                      Left\n' ...
    'Mag.             Prob.               Mag.             Prob.']))
%ylabel('Left                              Right')
set(get(gca,'YLabel'),'FontSize',14);