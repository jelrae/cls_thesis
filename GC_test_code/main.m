%
% Example code to compute spectral Granger causality with the MVGC toolbox
%
% The toolbox needs to be in the path. The command 'startup' will then
% generate a command line output and then the MVGC routines can be used.
%
% Jorge Mejias, May 2021
%
% -------------------
 

format short;clear all;
close all;clc;rng(938197);
startup; %to initialize the MVGC toolbox (usually not necessary)
load exampledata.mat;
%exampledata.mat contains the time series that you want to analyze: X0a,
%with the first dimension being the number of cortical areas (2), the
%second being the number of measurements across time (~18k) for each trial,
%and the third being the number of trials (40). You probably won't need to
%worry about other parameters there.


%Granger causality analysis:
fs=round(1/(par.binx*par.dt)); %sampling frequency

f0a=granger(X0a,100,fs,30,1e4); %granger is a (modified) example code from MVGC toolbox


%coherence anaysis:
window=1000;overlap=round(0.5*window);freqdisplay=2:2:100;
s10a=squeeze(X0a(1,:,:));m10a=squeeze(X0a(2,:,:));
for i=1:size(s10a,2)
    [pcoh0a(i,:),fcoh0a(i,:)]=mscohere(s10a(:,i),m10a(:,i),window,overlap,freqdisplay,fs);
end



%plot
figure('Position',[50,50,1000,400]);
blue=[.3 .8 .9];red=[1 .1 .1];resbin=1; 

%GC plot:
subplot(1,2,1);
dt=par.binx*par.dt;
z2to1a=squeeze(f0a(1,2,:));
z1to2a=squeeze(f0a(2,1,:));
frequ0=1:1:length(z1to2a);
nyq=2*length(z1to2a)*dt;
frequ0=frequ0./nyq;
Ntrials=size(fcoh0a,1);

frequ=frequ0(1:resbin:end);
GC1to2a=z1to2a(1:resbin:end);
GC2to1a=z2to1a(1:resbin:end);

plot(frequ,GC1to2a,'Color',blue,'LineWidth',3);hold on;
plot(frequ,GC2to1a,'Color',red,'LineWidth',3);hold on;
set(gca,'FontSize',24,'LineWidth',5,'TickLength',[0.02 0.02])
set(gca,'box','off');legend('FF','FB');legend boxoff;
xlim([0 80]);zgc=1.1*max(max(max(z2to1a),max(z1to2a)));
ylim([0 zgc]);title('V1-V4');set(gca,'Layer','top');
ylabel('Granger causality');xlabel('Frequency (Hz)');

%coherence plot:
subplot(1,2,2)
cohpa=mean(pcoh0a,1);cohpsiga=std(pcoh0a,1);cohfa=mean(fcoh0a,1);
myeb2(cohfa,cohpa,cohpsiga,[.2 .8 .6]);hold on;
set(gca,'FontSize',24,'LineWidth',5,'TickLength',[0.02 0.02]);
set(gca,'box','off');xlabel('Frequency (Hz)');ylabel('Coherence');
xlim([0 80]);zc=1.1*max(max(cohpa+cohpsiga));
ylim([0 zc]);title('V1-V4');



