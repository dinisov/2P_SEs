function data = filterChannel(data,degrees,filt_width)
%filterChannel Summary of this function goes here
%   Detailed explanation goes here

    % filter data (remove larger trends in time series)
    traceData = squeeze(mean(mean(data,1),2));

%     figure; plot(traceData);
    
    data = data - reshape(sgolayfilt(traceData,degrees,filt_width), [1 1 length(traceData)]);
    
%     figure; plot(squeeze(mean(mean(data,1),2)));

end

