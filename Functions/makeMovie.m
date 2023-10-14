function makeMovie(data, filename, isCorr)
%makeMovie Summary of this function goes here
%   Detailed explanation goes here
    VidObj = VideoWriter(filename, 'Uncompressed AVI');

    VidObj.FrameRate = 2; %set your frame rate
    open(VidObj);
    for f = 1:size(data, 3)
        if isCorr
            %if we're dealing with corr coeffs then map [-1 1] to [0 1]
            writeVideo(VidObj, ind2rgb(uint8(255 * (1+data(:,:,f))/2), jet(256)) );
        else
            %otherwise make video of data pre-mapped to [0 1]
            writeVideo(VidObj, ind2rgb(uint8(255 * data(:,:,f)), jet(256)));
        end
    end
    close(VidObj);
       
end