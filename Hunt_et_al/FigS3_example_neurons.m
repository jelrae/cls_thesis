% make figure S3, with different example single unit responses
% relies upon plot_neurons_average_peristim.m
clear all; close all;
bd = get_basedir;

%% define areas. indices 1 = ACC; 2 = LPFC 3 = OFC

load('all_units_info.mat','all_units');
for i = 1:length(all_units)
    area(i) = all_units{i}.UnitInfo.area;
    average_firing_rate(i) =  all_units{i}.UnitInfo.average_firing_rate;
end
area_index = [3 2 1];
color_index = {[0.2081 0.1663 0.5292] [0.1986 0.7214 0.6310] [0.7763 0.7831 0.0538]};
area_name = { 'ACC' 'DLPFC' 'OFC'};

area([87 find(average_firing_rate<1)]) = nan; %bad cells

%% make all figures

for u = [94 124 200 591 681]; %some units that have clean/illustrative peri-stimulus PSTHs
    if area(u)>=1&&area(u)<=3 %only draw units that come from OFC, ACC, DLPFC
        figure;
        clf;
        out = plot_neurons_average_peristim(u,[]);fields_to_workspace(out);
        this_cell_area = area_name{area(u)};

        
        out_fname = fullfile(bd,sprintf('%s_unit%04.0f.eps',this_cell_area,u));
    
        out_fname2 = fullfile(bd,sprintf('%s_unit%04.0f.eps',this_cell_area,u));
        
        
        %colours for plotting
        cc = [[0.1:0.2:0.9]' 0.5*ones(5,1) [0.9:-0.2:0.1]'];

        ymin = 0; % min(min([Cue1LPmn Cue1RMmn Cue1RPmn Cue1LMmn]))/1.1;
        ymax = max(max([Cue1LPmn Cue1RMmn Cue1RPmn Cue1LMmn]))*1.3;
        
        subplot(2,2,1);hold on;
        title('Cue 1 Left Probability');tidyfig;xlim([-100 400]);
        for i = 1:5;
            p(i) = plot(timeind,Cue1LPmn(i,:),'Color',cc(i,:),'LineWidth',2);
        end
        l = legend(p(5:-1:1),{'Best', ' ', ' ', ' ', 'Worst'},'Location','northwest');
        set(l,'Box','off');
        ylim([ymin ymax]);xlim([-100 400]);ylabel('Firing rate (Hz)');
        tidyfig(gcf,10);
        
        subplot(2,2,2);hold on;
        title('Cue 1 Right Probability');tidyfig;xlim([-100 400]);
        for i = 1:5;
            p(i) = plot(timeind,Cue1RPmn(i,:),'Color',cc(i,:),'LineWidth',2);
        end

        ylim([ymin ymax]);xlim([-100 400]);ylabel('Firing rate (Hz)');
        tidyfig(gcf,10);
        
        subplot(2,2,3);hold on;
        title('Cue 1 Left Magnitude');tidyfig;xlim([-100 400]);
        for i = 1:5;
            p(i) = plot(timeind,Cue1LMmn(i,:),'Color',cc(i,:),'LineWidth',2);
        end

        ylim([ymin ymax]);xlim([-100 400]);ylabel('Firing rate (Hz)');
        tidyfig(gcf,10);
        
        
        subplot(2,2,4);hold on;
        title('Cue 1 Right Magnitude');tidyfig;xlim([-100 400]);
        for i = 1:5;
            p(i) = plot(timeind,Cue1RMmn(i,:),'Color',cc(i,:),'LineWidth',2);
        end

        ylim([ymin ymax]);xlim([-100 400]);ylabel('Firing rate (Hz)');
        tidyfig(gcf,10);
        saveas(gcf,out_fname,'epsc');
        
        figure;clf;
        tmp(1,:,:) = Cue1LPmn; 
        tmp(2,:,:) = Cue1LMmn;
        tmp(3,:,:) = Cue1RPmn;
        tmp(4,:,:) = Cue1RMmn;
        bar(mean(tmp(:,:,501:end),3));
        set(gca,'XTickLabel',{'Left Prob','Left Mag','Right Prob','Right Mag'});
        ylabel('Firing Rate (Hz)');
        tidyfig(gcf,10);
        %saveas(gcf,out_fname2,'epsc');
    end
end