%

close all; clear;

rdmDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
finalSize = [128 128];

%%

chosenFlies = [25:41];
chosenBlocks = []; % leave empty if aligning all blocks for one fly

flagParamSaveList = who;
flagParamSaveList = [flagParamSaveList;'flagParamSaveList';'fly'];

% this level is the list of blocks
%for fly =[185:187]%height(flies)
for fly = chosenFlies
    clearvars('-except', flagParamSaveList{:});
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == fly,:);

    if ~isempty(chosenBlocks)
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks),:);
    end

    nBlocks = height(thisFlyBlocks);
    
    %currentDate = char(datetime(currentFly.Date,'Format','dMMMyy'));
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    %currentFlyDirectory = ['fly' num2str(currentFly.FlyOnDay) '_exp' num2str(currentFly.Block) '_' currentDate];
    for b = 1:nBlocks
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate]; %Borrowed from pre_process

        currentRDMDirectory = fullfile(rdmDirectory,currentDate,flyID);

        %disp(currentFlyDirectory);
        disp(['Fly: ',flyID]);
        
        codeStartTime = posixtime(datetime('now'));
        loadReduceSave(currentRDMDirectory, 'green_channel.raw', currentBlock, finalSize); %currentBlock corresponds to currentFly

        codeEndTime = posixtime(datetime('now'));
        MET = codeEndTime - codeStartTime;
        disp(['-- Total time to process: ',num2str(MET),'s --']) 
    end
end

function loadReduceSave(RDMDirectory, file, fly, finalSize)

    fileRDM = fullfile(RDMDirectory, file);
    
    imageSize = [fly.pixelX fly.pixelY];
    nFrames = fly.realFrames;

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
    
    reducedFileRDM = fullfile(RDMDirectory, ['green_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']);
    
    tic
    disp('Saving');
    % save green channel; do not compress we care about speed not size
    save(reducedFileRDM, 'rData','-v7.3','-nocompression');
    toc
    
    % free up the memory 
    clear;
    
end