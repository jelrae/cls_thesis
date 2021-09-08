addpath('/home/jordan/common/matlab/fieldtrip-20210411')

%% Set a default trialdef

trialdef.startpath = '/home/jordan/preproc_redo_test/data_kurt_super_raw/';
trialdef.color     = 'both';
trialdef.stim      = {'changeCue'};
trialdef.offset    = 300; %100;
trialdef.trialdur  = 1300;%600;%800; %this is the minumum length of the trial
trialdef.nStim     = 2;
trialdef.target    = 'both';
trialdef.cnd       = {'AttIn' 'AttOut'};
trialdef.redef     = 'yes';
trialdef.timwin    = 0.5;
trialdef.behav     = 'no';
trialdef.electr    = 'yes';
trialdef.trialfun  = 'trialfun_attask7';

if strcmp(trialdef.electr,'yes')
    flag_electro = true;
else
    flag_electro = false;
end

%% Extract parameters from a .txt file
% start_path = trialdef.startpath;
% fid = fopen([start_path filename],'r');
% dat_loc  = fgetl(fid);
% raw_path = fgetl(fid);
% monkey   = fgetl(fid);
% BM = textscan(fid,'%s%s%s','delimiter',',');
% fclose(fid);
% v_session = str2num(char(BM{1}));
% numses = length(v_session);
% v_maxsubses = str2num(char(BM{2}));
% v_task = str2num(char(BM{3}));
% clear BM;

start_path = trialdef.startpath;
raw_path = '/home/jordan/preproc_redo_test/new_raw_data';
monkey = 'ku';
numses = 1;
v_maxsubses = [2];
v_session = [40];
v_task = [006];

%% get data structure
%eval(dat_loc)
addpath('/home/jordan/preproc_redo_test/data_kurt_super_raw')
basic_data_path = '/home/jordan/preproc_redo_test/data_kurt_super_raw';
cfg = [];
cfg.trialfun           = trialdef.trialfun;
cfg.dataformat         = 'combined_ds';
cfg.headerformat       = 'combined_ds';
cfg.trialdef           = trialdef;
cfg.dftfilter          = 'yes';

if strcmp(monkey,'pe')
    cfg.dftfreq            = [50 100 150 200];
    cfg.padding            = 10;
else
     cfg.padding            = 10;
     cfg.dftfreq            = [49.7:0.1:50.4 99.7:0.1:100.4 119.7:0.1:120.3 150];
end
     %cfg.paddir             = 'left';