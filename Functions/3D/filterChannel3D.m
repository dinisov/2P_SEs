function data = filterChannel3D(data,degrees,filt_width)
%filterChannel Savitsky-Golay filter the 2P data
%   this is the 3D version
    
    % average out the Z-direction
    data_aux = squeeze(mean(data,3));

    % mean of every volume (time point) over all trials
    traceData = squeeze(mean(mean(data_aux,1),2));
    
    % subtract filter from data and add mean (so it does not go below 0)
    data = data - reshape(sgolayfilt(traceData,degrees,filt_width), [1 1 1 length(traceData)]) + mean(traceData);
    
%     figure; plot(squeeze(mean(mean(mean(data,3),1),2)));

end