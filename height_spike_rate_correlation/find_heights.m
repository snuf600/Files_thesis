function [heights] = find_heights(zData,t,events)
fs = 1/(t(2)-t(1));
onsets = events.onsets(events.onsets<1600 & events.onsets>900)-900;
moving_window = round(fs*5);
valid_points = ones(length(zData),1);%Valid where no spike within window
for i = 1:length(onsets)
    %First find allocated datapoint to timestamp
    [~,idx] = min(abs(t-onsets(i)));
    %Use this index to set window around it to 0
    valid_points(max(idx-moving_window,1:min(length(valid_points),idx+moving_window))) = 0;
end
% Filter out real noisy parts
moving_std = movstd(zData,100);
moving_std_thr = 0.20;
thresholded_std = moving_std > moving_std_thr;
valid_points_std = ones(length(zData),1);
for i = 1:length(zData)
    if i<moving_window
        valid_points_std(i) = sum(ismember(onsets,1:i+moving_window));
    elseif i>length(zData)-moving_window
        valid_points_std(i) = sum(thresholded_std(i-moving_window+1:end));
    else
        valid_points_std(i) = sum(thresholded_std(i-moving_window+1:i+moving_window));
    end
end
valid_points_std = valid_points_std == 0;%Where sum == 0, no spikes are nearby
valid_points = valid_points_std.*valid_points;


heights = valid_points.*zData;
heights(heights == 0) = nan;


end
