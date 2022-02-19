function p_info = PeakAlignSpectrum(start_f_loc, end_f_loc, w, data)

%%%
%start_f_loc (int) = Location in the array where lower bound freq is
%end_f_loc (int) = Location in the array where upper bound freq is
%data (trials x freq res) = array data over which the peaks will be
%searched for
%%%

%% Initialize the storing structure
p_info = [];
p_info.location = [];
p_info.value = [];
p_info.width = [];
p_info.no_peak = [];

%%Search for the peaks and store locations and sizes
for i =1:size(data,1)
	[pk, pl, pw, pp] = findpeaks(data(i,start_f_loc:end_f_loc));
    if ~isempty(pk)
        [vm, lm] = max(pk);
        p_info.value(end+1) = vm;
        p_info.location(end+1) = start_f_loc - 1 + pl(lm);
        p_info.width(end+1) = pw(lm);
    else
        p_info.value(end+1) = 0;
        p_info.location(end+1) = 0;
        p_info.width(end+1) = 0;
        p_info.no_peak(end+1) = i;
    end
end

%%Set the peak locations for trials without peaks to the average peak location

p_info.ave_loc = round(sum(p_info.location)/nnz(p_info.location));
p_info.location(p_info.location==0) = p_info.ave_loc;

%% Align peaks and average
% w = floor(end_f_loc - start_f_loc)/2;

p_info.peak = zeros(1,2*w+1);

for i = 1:size(p_info.location,2)
    start_loc = p_info.location(i) - w;
    end_loc = p_info.location(i) + w;
    p_info.peak = p_info.peak + data(i,start_loc:end_loc);
end

p_info.peak = p_info.peak./size(data,1);

% g_min = p_info.ave_loc - w;
% g_max = p_info.ave_loc + w;
% plot(coh_in.freq(g_min:g_max), g_peak);

end