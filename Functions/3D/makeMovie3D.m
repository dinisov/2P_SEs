function makeMovie3D(data, filename, isCorr)
%makeMovie Summary of this function goes here
%   Detailed explanation goes here
    VidObj = VideoWriter(filename, 'Uncompressed AVI');

    VidObj.FrameRate = 1.5; %set your frame rate
    open(VidObj);
    for f = 1:size(data, 4)
        plot3D(data(:,:,:,f), 'off'); 
        
        % set the colormap bounds to [-1 1]
        if isCorr
            caxis([-1 1]);
        end
        
        writeVideo(VidObj, getframe);
        close all;
    end
    close(VidObj);    
end