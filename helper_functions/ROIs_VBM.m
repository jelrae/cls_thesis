function rois_out = ROIs_VBM(cfg,monkey)

% Based on definitions of Vezoli, Bastos and Markov

% cfg.rois = list of ROI names that should be returned ex. cfg.rois =
% {'V1', 'V2'}
% cfg.channel = list of bad channel names that should be excluded from the
% ROIs to prevent ugliness down the road (keep in mind that this is the
% same as the normal cfg.channel, so if you put 'all', you will get all of
% the ROI, this is default, so put the 'all' flag and then '-E03-E02' format, so
% the negative sign gets rid of that channel.

cfg.rois = ft_getopt(cfg,'rois','all');
cfg.channel = ft_getopt(cfg,'channel','all');

if strcmp(cfg.rois,'all')
    rois = {'V1','V2', 'a8L', 'V4','DP','TEO', 'a8M', 'a7A','a7APG','a7B','TPt','a5','a1a2','F1','F2','F4'}; %note: this should be updated to reflect newest ROIs from hierarchy paper, 
                                                                                                  %and also grouped according to hierarchy and to cortical system (viz,Ss,M)
else
    rois = cfg.rois;
end

if isequal(monkey,'kurt') || isequal(monkey,1)
    
V1 = {'J02-J01',  'J14-J13', 'J15-J14', 'J16-J15', 'J17-J16', 'J18-J17', ...
        'J19-J18', 'J20-J19', 'J21-J20', 'J22-J21', 'J24-J23', 'K02-K01', ...
        'K03-K02', 'K04-K03', 'K05-K04', 'K06-K05', 'K07-K06', 'K08-K07', ...
        'K09-K08', 'K10-K09', 'K11-K10', 'K13-K12', 'K14-K13', 'K15-K14', ...
        'K16-K15', 'K17-K16', 'K18-K17', 'K19-K18', 'K20-K19', 'K21-K20', 'K22-K21'};

    V2 = {'J03-J02', 'J04-J03', 'J05-J04', 'J06-J05', 'J07-J06', 'J08-J07', 'J09-J08', 'J10-J09', 'J11-J10'};

    V4 = {'H07-H06', 'H08-H07', 'H09-H08', 'H10-H09', 'H12-H11', 'H13-H12', ...
        'H21-H20', 'H22-H21', 'H23-H22', 'H24-H23', 'H25-H24', 'H26-H25', ...
        'H27-H26', 'H28-H27', 'G12-G11', 'G24-G23', 'G25-G24', 'G26-G25', 'J12-J11'};

    DP = {'H02-H01', 'H03-H02', 'H04-H03', 'H05-H04', 'H06-H05', 'H16-H15',...
        'H17-H16', 'H18-H17', 'H19-H18', 'H20-H19'};
    TEO = {'G13-G12', 'G14-G13', 'G27-G26', 'G28-G27', 'H14-H13'};
    a7A = {'F19-F18', 'G04-G03', 'G05-G04', 'G06-G05', 'G17-G16', 'G18-G17', 'G19-G18', 'G20-G19'};
%     a8 = {'A07-A06', 'A08-A07', 'A09-A08','N02-N01', 'N06-N05','N07-N06', 'O04-O03', 'O06-O05', 'O07-O06', 'O08-O07' ... %% noisy channels
%     'N03-N02', 'N04-N03', 'N08-N07', 'O02-O01', 'O03-O02',}; %putative area 8l = O0x channels; 8m = N0x channels || need further verification
    a8M = {'A07-A06', 'N02-N01', 'N06-N05','N07-N06', 'N03-N02', 'N04-N03', 'N08-N07'};
    a8L = {'A08-A07', 'A09-A08', 'O04-O03', 'O06-O05', 'O07-O06', 'O08-O07','O02-O01', 'O03-O02'}; 
    a7APG = {'E07-E06', 'E15-E14', 'E16-E15', 'F06-F05', 'F07-F06', 'F08-F07', ...
       'F20-F19', 'F21-F20', 'F22-F21', 'G08-G07', 'G09-G08', 'G21-G20', 'G22-G21', 'G23-G22'}; 
   a7B = {'D18-D17', 'E08-E07', 'E09-E08', 'E17-E16', 'E18-E17', 'F09-F08', 'F10-F09', 'F23-F22'};
   TPt = {'F11-F10', 'F24-F23', 'F25-F24', 'G10-G09', 'G11-G10'};
   PBc = {'F12-F11', 'F13-F12', 'F14-F13', 'F26-F25', 'F27-F26', 'F28-F27'};  %ParaBelt caudal or TAa
   a5 = {'E04-E03', 'E05-E04', 'E06-E05', 'E11-E10', 'E12-E11', 'E13-E12', ...
       'E14-E13', 'F02-F01', 'F04-F03', 'F05-F04', ...
       'F16-F15', 'F17-F16', 'F18-F17', 'G02-G01', 'G03-G02', 'G16-G15'}; %removed 'F03-F02' due to error
   F1 = {'A11-A10', 'B02-B01', 'B03-B02', 'B04-B03', 'B05-B04', 'B06-B05', ...
       'B07-B06', 'B11-B10', 'B12-B11', 'B13-B12', 'B14-B13', ...
       'B15-B14', 'B16-B15', 'B17-B16', 'C02-C01', 'C03-C02', ...
       'C04-C03', 'C05-C04', 'C06-C05', 'C07-C06', 'C08-C07', 'C11-C10', 'C12-C11', 'C13-C12', 'C14-C13'}; %removed 'B18-B17' due to error
   F2 = {'A02-A01', 'A03-A02', 'A04-A03', 'A05-A04', 'A06-A05', 'A12-A11', ...
       'A13-A12', 'A14-A13', 'L02-L01', 'L03-L02', 'L04-L03', 'L06-L05', ...
       'L06-L05', 'L07-L06', 'L08-L07', 'M02-M01', 'M03-M02', 'M04-M03', 'M06-M05', 'M07-M06'}; %removed 'M08-M07' which is a "bad channel" in kurt
   F4 = {'A15-A14', 'A16-A15', 'A17-A16', 'A18-A17',  'B08-B07', 'B09-B08'};
   a1a2 = {'C09-C08', 'C18-C17', 'C17-C16', 'C16-C15', 'C15-C14', 'D02-D01', ...
       'D03-D02', 'D04-D03', 'D05-D04', 'D06-D05', 'D08-D07', 'D09-D08', ...
       'D11-D10', 'D12-D11', 'D13-D12', 'D14-D13', 'D15-D14', 'D16-D15', 'D17-D16', 'E02-E01', 'E03-E02'}; 
   
    
elseif isequal(monkey,'pele') || isequal(monkey,2)
    
    V1 = {'H16-H15', 'H17-H16', 'H22-H21', 'H23-H22', 'H24-H23', 'H25-H24', 'H26-H25', 'H27-H26', 'H28-H27', ...
        'J02-J01', 'J03-J02', 'J04-J03', 'J05-J04', 'J06-J05', 'J07-J06', 'J08-J07', 'J09-J08', 'J10-J09', 'J11-J10', ...
        'J12-J11', 'J14-J13', 'J15-J14', 'J16-J15', 'J17-J16', 'J18-J17', 'J19-J18', 'J20-J19', 'J21-J20', 'J22-J21', 'J24-J23', ...
        'K02-K01', 'K03-K02', 'K04-K03', 'K05-K04', 'K06-K05', 'K07-K06', 'K08-K07', 'K09-K08', 'K10-K09', 'K11-K10', 'K13-K12',  ...
        'K14-K13', 'K15-K14','K16-K15', 'K17-K16', 'K18-K17', 'K19-K18', 'K20-K19', 'K21-K20', 'K22-K21'};

    V2 = {'H04-H03','H05-H04','H06-H05', 'H07-H06', 'H08-H07', 'H09-H08', 'H10-H09', 'H12-H11', 'H13-H12', 'H14-H13', ...
        'H18-H17', 'H19-H18', 'H20-H19', 'H21-H20'};

    V4 =  { 'F13-F12', 'F25-F24', 'F26-F25', 'F27-F26', 'G06-G05', 'G10-G09', 'G11-G10', 'G12-G11', 'G13-G12',  ...
           'G20-G19', 'G21-G20', 'G22-G21', 'G23-G22', 'G24-G23', ...
           'G25-G24', 'G26-G25', 'G27-G26', 'G28-G27'};  %without putative TEO channels || TEO covering need further verification

    DP = {'G04-G03', 'G05-G04', 'G16-G15', 'G17-G16', 'G18-G17', 'G19-G18', 'H02-H01', 'H03-H02'};
    
    a7A = {'E05-E04', 'E13-E12', 'E14-E13', 'F04-F03', 'F05-F04', 'F17-F16', 'F18-F17', 'F19-F18', 'G02-G01', 'G03-G02'};
    a8 = {'O08-O07', 'N07-N06', 'N08-N07', 'O03-O02', 'O04-O03'}; %putative area 8l = O0x channels; 8m = N0x channels || need further verification
    a8L = {'O08-O07', 'O03-O02', 'O04-O03'};
    a8M = {'N07-N06', 'N08-N07'};
    TEO = {'F28-F27','G14-G13','F14-F13'};
    a7APG ={'D15-D14', 'D16-D15', 'E07-E06', 'E08-E07', 'E06-E05', 'E15-E14', 'E16-E15', 'E17-E16', 'F06-F05', 'F07-F06', 'F08-F07', 'F20-F19', ...
        'F21-F20' ,'F22-F21', 'G08-G07'}; 
    a7B = {'C09-C08', 'C17-C16', 'C18-C17', 'D08-D07', 'D09-D08', 'D17-D16', 'D18-D17', 'E09-E08', 'E18-E17', 'F09-F08', 'F10-F09', 'F23-F22'};
    TPt = {'F11-F10', 'F12-F11', 'F24-F23', 'G09-G08'};
    a5 = {'D03-D02', 'D04-D03', 'D05-D04', 'D06-D05', 'D11-D10', 'D12-D11', 'D13-D12', 'D14-D13', 'E02-E01', 'E03-E02', 'E04-E03', 'E11-E10', 'E12-E11', 'F02-F01', 'F16-F15'};
    a1a2 = {'A18-A17', 'B06-B05', 'B07-B06', 'B08-B07', 'B09-B08', 'B14-B13', 'B15-B14', 'B16-B15', 'B17-B16', 'C03-C02',  ...
        'C04-C03', 'C05-C04', 'C06-C05', 'C07-C06', 'C08-C07', 'C11-C10', 'C12-C11', 'C13-C12', 'C14-C13', 'C15-C14', 'C16-C15', 'D02-D01'}; %removed 'B18-B17' due to error
    F1 = {'A02-A01', 'A03-A02', 'A04-A03', 'A05-A04', 'A06-A05', 'A07-A06', 'A08-A07', 'A09-A08', 'A11-A10', 'A12-A11', 'A13-A12', 'A14-A13', 'A15-A14', 'A16-A15', 'A17-A16', 'B02-B01',  ...
        'B03-B02', 'B04-B03', 'B05-B04', 'B11-B10', 'B12-B11', 'B13-B12', 'C02-C01'};
    F2 = {'L02-L01', 'L03-L02', 'L04-L03', 'L06-L05', 'L07-L06', 'L08-L07', 'M02-M01', 'M03-M02', 'M04-M03', 'M06-M05', 'M07-M06',  ...
        'N02-N01', 'N03-N02', 'N04-N03'}; %removed 'N05-N04' due to error, removed  'M08-M07' due to NaNs
    F4 = {'O02-O01', 'O06-O05', 'O07-O06', 'N06-N05'};
end

for ii=1:length(rois)
    % remove bad channels
    eval([rois{ii} '= ft_channelselection(cfg.channel,' rois{ii} ')']);
    rois_out(ii).name = rois(ii);
    rois_out(ii).label = eval([rois{ii} '(:);']);
end

%load(['C:\DropBox\monkey\',monkey,'\',monkey,'_channel_assignment']);

for k = 1:numel(rois)
    tmp = rois_out(k).label;
    uni = {};
    for m = 1:numel(tmp)
        uni = [uni; tokenize(tmp{m}, '-')'];
    end
    uni = unique(uni);
    rois_out(k).unipolars = uni;
%     headstage = cell(numel(uni),1);
%     for m = 1:numel(uni)
%         ix = strmatch(uni{m}, channel_assignment(:,1));
%         headstage{m} = channel_assignment{ix,4};
%     end
%     rois_out(k).headstage = headstage;
end
