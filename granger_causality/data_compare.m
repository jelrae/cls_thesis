clear all;

load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_v1_AttIn.mat')
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_kurt_data\bipolar_post_data\kurt_p_v4_AttIn.mat')
kurt_v1_in = v1_AttIn.trial;
kurt_v4_in = v4_AttIn.trial;
kurt_v1_in = cat(3,kurt_v1_in{:});
kurt_v4_in = cat(3,kurt_v4_in{:});
clear v1_AttIn
clear v4_AttIn
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v1_AttIn.mat')
load('C:\Users\Jordan\Documents\cls_thesis\neuro_thesis\data\monkey_pele_data\bipolar_post_data\pele_p_v4_AttIn.mat')
pele_v1_in = v1_AttIn.trial;
pele_v4_in = v4_AttIn.trial;
pele_v1_in = cat(3,pele_v1_in{:});
pele_v4_in = cat(3,pele_v4_in{:});
clear v1_AttIn
clear v4_AttIn

close all;

chan = 2;

figure('Position',[50,50,1000,400]);
title('5 trials Kurt and Pele, V1 and V4');

x_range = 1:1:501;

subplot(2,2,1);
plot(x_range, squeeze(kurt_v1_in(chan,:,1:10))); hold on;
title('Problem child Kurt V1');

subplot(2,2,2);
plot(x_range, squeeze(kurt_v4_in(chan,:,1:10))); hold on;
title('Problem child Kurt V4');

subplot(2,2,3);
plot(x_range, squeeze(pele_v1_in(chan,:,1:10))); hold on;
title('Pele V1');

subplot(2,2,4);
plot(x_range, squeeze(pele_v4_in(chan,:,1:10))); hold on;
title('Pele V4');