%

close all; clear;

rdmDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

scratchDirectory = '../../2P Data';

flies = readtable('../../2P Record/2P_record');

%get rid of excluded flies
flies = flies(~logical(flies.Exclude),:);

% the numbers here should be the original size divided by some power of 2
finalSize = [128 128];

%%

% this level is the list of blocks
for fly =107:120%height(flies)
    
    currentFly = flies(fly,:);
    
    currentDate = char(datetime(currentFly.Date,'Format','dMMMyy'));
    currentFlyDirectory = ['fly' num2str(currentFly.FlyOnDay) '_exp' num2str(currentFly.Block) '_' currentDate];
    
    currentRDMDirectory = fullfile(rdmDirectory,currentDate,currentFlyDirectory);
    currentScratchDirectory = fullfile(scratchDirectory,currentDate,currentFlyDirectory);

    disp(currentFlyDirectory);

    loadReduceSave(currentRDMDirectory,currentScratchDirectory, 'green_channel.raw', currentFly, finalSize);
%         loadReduceSave([currentDirectory '\red_channel.raw'], currentFly, finalSize);

end

function loadReduceSave(RDMDirectory, scratchDirectory, file, fly, finalSize)

    fileRDM = fullfile(RDMDirectory, file);
%     fileScratch = fullfile(scratchDirectory, file);
    
    imageSize = [fly.pixelX fly.pixelY];
    nFrames = fly.realFrames;
    
    if ~exist(scratchDirectory,'dir')
       mkdir(scratchDirectory); 
    end
    
    % move file to scratch before loading
%     disp('Copying file from RDM');
%     tic
%     copyfile(fileRDM,fileScratch);
%     toc
    
    tic
    disp('Loading data');
    % load green channel
    fid = fopen(fileRDM, 'r','b');
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
    
%     reducedFileScratch = fullfile(scratchDirectory, ['green_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']);
    reducedFileRDM = fullfile(RDMDirectory, ['green_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']);
    
    tic
    disp('Saving');
    % save green channel; do not compress we care about speed not size
%     save([path(1:end-4) '_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat'], 'rData','-v7.3','-nocompression');
    save(reducedFileRDM, 'rData','-v7.3','-nocompression');
    toc
    
    %move file back to RDM
%     disp('Moving back to RDM');
%     tic
%     copyfile(reducedFileScratch,reducedFileRDM);
%     toc
    
    %delete file in scratch
%     disp('Deleting files in scratch');
%     tic
%     delete(fileScratch);
%     delete(reducedFileScratch);
%     toc
    
    % free up the memory 
    clear('data');
    
end