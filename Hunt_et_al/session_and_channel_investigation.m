clear

addpath('/home/jordan/neuro_thesis/cls_thesis/helper_functions/')
set_paths_roi;

clear all
clc
ft_defaults

%% set up base directory where data/scripts are stored
[bd] = get_basedir;
addpath(genpath(fullfile(bd)));

datadir = fullfile(bd,'neuronal_data');
filenames=dir(fullfile(datadir,'*.mat'));

load('frank_area');frank_area=area_index; clear area_index;
load('miles_area');

%% main loop over each unit from information gathering experiment

new_correct=100; %this is the minimum length of time for a picture to be considered 'viewed', used below

skip_channels = [];

session_and_channels = [];
length(filenames)
for u=1:length(filenames) %u = 'unit'
    tic
    %% load in the data
    load(filenames(u).name);
    u
    filenames(u);


    %% get information about session/channel from filename
    sb=filenames(u).name(1); %subject name, F=frank, M=miles
    % change name(2:5) to 2:4
    snum=str2num(filenames(u).name(2:4)); %session number, 1-49 for Miles, 1-36 for Frank
    chnum=str2num(filenames(u).name(strfind(filenames(u).name,'C')+1:strfind(filenames(u).name,'U')-2)); %channel number

    %% For future use
    %find which row of miles_area/frank_area this corresponds to, and then
    %find out which brain region was recorded from
    if strcmp(sb,'F'),
        str_area=find(frank_area(:,1)==snum & frank_area(:,2)==chnum);
        brain_region=frank_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
    else
        str_area=find(miles_area(:,1)==snum&miles_area(:,2)==chnum);
        brain_region=miles_area(str_area,3); %brain region (1 = ACC; 2 = DLPFC; 3 = OFC; 4 = Other; 5 = Other; 6 = Other)
    end
    
    if ismember(brain_region, [1,2,3])

        %% To just figure out whats in this bitch
        session_and_channels(end+1,:) = [snum, chnum, brain_region];
        
    else
        skip_channels(end+1,:) = [snum, chnum, brain_region];
    end
    
end

unique_channels = unique(session_and_channels(:,2));
unique_sessions = unique(session_and_channels(:,1));

for i=1:length(unique_channels)
    unique_channels(i,2) = sum(session_and_channels(:,2) == unique_channels(i,1));
end