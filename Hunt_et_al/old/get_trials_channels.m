%% Load packages
addpath('D:/MATLAB/fieldtrip-20200515')
ft_defaults

%% Reference data
load('D:/MATLAB/forSiobhan/UvA_internship/Kurt/Data/ku_dataAttIn.mat');
load('D:/MATLAB/forSiobhan/UvA_internship/Kurt/Data/ku_dataAttOut.mat');
%% Get file list
f_dir = 'D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/f_baseline/';
m_dir = 'D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/m_baseline/';
dir_data = 'D:/2020/Amsterdam_internship/gc_hierarchies/Hunt_etal/Data/unprocessed/';
addpath(genpath(fullfile(dir_data)));

% See get_lfps.m: generates the variables stored in filenames
filenames=dir(fullfile(dir_data,'LTH_processed_data_*'));
%% General parameters 
num_chan = 24;
% cues used to extract subset of the data
start_cue = 11; 
end_cue = 12;
%% Select monkey information and parameters



%% Loop through directory
x = 0;
trial = zeros(num_chan, 30);
f_session_3 = {};
f_session_13 = {};
f_session_8 = {};
f_session_2 = {};

m_session_3 = {};
m_session_13 = {};
m_session_8 = {};
m_session_2 = {};

f_all_sessions = {};
m_all_sessions = {};

file_length = length(filenames);
for i=1:file_length
    file = filenames(i).name;
    load(file)
    
    if strcmp(sb, 'F')
        num_sessions = 36;
    elseif strcmp(sb, 'M')
        num_sessions = 49;
    end
%     eval(sprintf('freq_tap%d_%s=tmp;',tapsmofrq(z),cnd{n}));
%     tmpFr.avgChOI = pow2_cfg.avgChOI;
%     eval(sprintf('frdescr_%s_tap%d_%s=tmpFr;',roi(a).name,tapsmofrq(z),cnd{n}));
    
    for session=1:num_sessions
        if session == snum
            for k=1:length(InfoLFP)
                current_trial = InfoLFP{k};
         
                x = x + 1;
                % get lfp baseline
                for t=1:length(TimingStrobe)
                    if k == t
                        for ts=1:size(TimingStrobe{t},1)
                            if TimingStrobe{t}(ts,1) == start_cue
                                start_baseline = TimingStrobe{t}(ts, 2);
                            else
                                start_baseline = NaN;
                            end
                            if TimingStrobe{t}(ts,1) == start_cue && TimingStrobe{t}(ts + 1,1) == end_cue
                                end_baseline = TimingStrobe{t}(ts + 1, 2);
                            else
                                end_baseline = NaN;
                            end

                            if isnan(start_baseline) || isnan(end_baseline)
                                continue
                            elseif not(isnan(start_baseline)) && not(isnan(end_baseline))
                                baseline_lfp = InfoLFP{t}(start_baseline:end_baseline);
                            end
                            % Create loop for varying sessions

                            if snum == 3 && strcmp(sb, 'F')
                                f_session_3{k}(chnum, :) = baseline_lfp(1:500)';
                            elseif snum == 3 && strcmp(sb, 'M')
                                m_session_3{k}(chnum, :) = baseline_lfp(1:500)';

                            elseif snum == 8 strcmp(sb, 'F')
                                f_session_8{k}(chnum, :) = baseline_lfp(1:500)';
                            elseif snum == 8 strcmp(sb, 'M')
                                m_session_8{k}(chnum, :) = baseline_lfp(1:500)';
                                
                            elseif snum == 2 strcmp(sb, 'F')
                                f_session_2{k}(chnum, :) = baseline_lfp(1:500)';
                            elseif snum == 2 strcmp(sb, 'M')
                                m_session_2{k}(chnum, :) = baseline_lfp(1:500)';
                                
                            elseif snum == 13 strcmp(sb, 'F')
                                f_session_13{k}(chnum, :) = baseline_lfp(1:500)';
                            elseif snum == 13 strcmp(sb, 'M')
                                m_session_13{k}(chnum, :) = baseline_lfp(1:500)';    
                            end
                        end
                    end
                end

            end
        end
    end 
end

%% Create cell with time info: time axis [1*Ntime double] per trial
session_of_interest = f_session_2;
% remove empty trials
for tr=1:size(session_of_interest,2)
    if not(iscell(session_of_interest{1, tr}))
        continue
    end
    if isempty(session_of_interest{1, tr}) 
        session_of_interest(tr) = [];
    end
    
end
        
num_trial =(size(session_of_interest,2))- 2;
num_time = cell(num_trial,1);
sample_info = cell(num_trial,1);
for time=1:num_trial
    num_time{time, 1} = linspace(0, size(baseline_lfp(1:500), 1), size(baseline_lfp(1:500), 1));
end
for sinfo=1:num_trial
    sample_info{sinfo, 1} = 1;
    sample_info{sinfo, 2} = size(baseline_lfp(1:500), 1);
end

%% Create label: {num_chanx1 cell}
label = cell(num_chan,1);
for i=1:num_chan
    label{i,1} = int2str(i);
end
%% Create cfg

baseline = [];
baseline.label = label;
baseline.trial = session_of_interest(:, 3:end);
baseline.time = num_time';
baseline.sampleinfo = sample_info;
baseline.dimord = 'chan_time';
baseline.fsample = 1000;

cfg = [];

cfg.keeptapers = 'yes';
cfg.average = 'yes';
cfg.method = 'mtmfft';
cfg.keeptapers = 'yes';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';%'no';
cfg.pad='maxperlen'; % Padding: not adding zeros
cfg.flag = 0;

tmp = ft_freqanalysis(cfg, baseline);