%

close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

flies = readtable('../2P Record/2P_record');

% the numbers here should be the original size divided by some power of 2
finalSize = [128 128];

%%

% this level is the list of dates
for fly = 31:31%height(flies)
    
    currentFly = flies(fly,:);
    
    if ~currentFly.Exclude
    
        currentDate = char(datetime(currentFly.Date,'Format','dMMMyy'));

        currentFlyDirectory = ['fly' num2str(currentFly.FlyOnDay) '_exp' num2str(currentFly.Block) '_' currentDate];

        currentDirectory = fullfile(mainDirectory,currentDate,currentFlyDirectory);

        disp(currentFlyDirectory);
        
        loadReduceSave([currentDirectory '\green_channel.raw'], currentFly, finalSize);

        loadReduceSave([currentDirectory '\red_channel.raw'], currentFly, finalSize);
        
    end

end

function loadReduceSave(path, fly, finalSize)

    imageSize = [fly.pixelX fly.pixelY];

    nFrames = fly.realFrames;
    
    tic
    disp('Loading data');
    % load green channel
    fid = fopen(path, 'r','b');
    data = fread(fid, 512*512*nFrames, 'uint16');
    fclose(fid);
    toc
    
    %rearrange
    data = permute(reshape(data, [imageSize nFrames]),[2 1 3]);

    rData = zeros([finalSize nFrames]);
    
    tic
    disp('Reducing size');
    %reduce green channel (nothing comes close to this in terms of speed)
    rData = imresize3(data,size(rData),'box');
    toc
    
    tic
    disp('Saving');
    % save green channel; do not compress we care about speed not size
    save([path(1:end-4) '_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat'], 'rData','-v7.3','-nocompression');
    toc
    
end