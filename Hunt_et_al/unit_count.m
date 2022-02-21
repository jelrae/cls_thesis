% script to return table of session ID, number of units recorded from
% ACC/OFC/DLPFC

% fsess and msess are recording sessions for Frank and Miles respectively
% rows are recording sessions
% columns are Session ID, number of ACC cells, number of OFC cells, number
% of DLPFC cells - only for cells with firing rate>=1

% for reporting in supplementary information

load('all_units_info.mat');

fcount = 0;
fsess = zeros(100,3);
mcount = 0;
msess = zeros(100,3);

count_all_units = 1; %set to 0 to include units with firing rate >=1 Hz (as in paper), set to 1 to include all units
for i = 1:724
    if all_units{i}.UnitInfo.isfrank
        fcount = fcount + 1;
        if all_units{i}.UnitInfo.average_firing_rate>=1||count_all_units
            switch all_units{i}.UnitInfo.area_name
                case 'ACC'
                    fsess(all_units{i}.UnitInfo.session_number,1) = fsess(all_units{i}.UnitInfo.session_number,1) + 1;
                case 'OFC'
                    fsess(all_units{i}.UnitInfo.session_number,2) = fsess(all_units{i}.UnitInfo.session_number,2) + 1;
                case 'DLPFC'
                    fsess(all_units{i}.UnitInfo.session_number,3) = fsess(all_units{i}.UnitInfo.session_number,3) + 1;
            end
        end
    else
        mcount = mcount + 1;
        if all_units{i}.UnitInfo.average_firing_rate>=1||count_all_units
            switch all_units{i}.UnitInfo.area_name
                case 'ACC'
                    msess(all_units{i}.UnitInfo.session_number,1) = msess(all_units{i}.UnitInfo.session_number,1) + 1;
                case 'OFC'
                    msess(all_units{i}.UnitInfo.session_number,2) = msess(all_units{i}.UnitInfo.session_number,2) + 1;
                case 'DLPFC'
                    msess(all_units{i}.UnitInfo.session_number,3) = msess(all_units{i}.UnitInfo.session_number,3) + 1;
            end
        end
    end
end

fsess = [[1:100]' fsess];
fsess(sum(fsess(:,2:4),2)==0,:) = [];

msess = [[1:100]' msess];
msess(sum(msess(:,2:4),2)==0,:) = [];

