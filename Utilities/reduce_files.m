%

close all; clear;

sourceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Matt\'; %Will also function as sinkDirectory if altSinkDirectory empty

altSinkDirectory = 'D:\group_vanswinderen\Matt\2p\2P_Data\'; %If empty, will save to sourceDirectory, otherwise will save to here

overrideAlreadyDone = 1; %Whether to override already processed block skipping
deletePrecursors = 1; 

blocks = readtable('../../2P Record/2P_record_MattNew');

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
finalSize = [128 128];

%%

chosenFlies = [7];
chosenBlocks = []; % leave empty if aligning all blocks for one fly

flagParamSaveList = who;
%flagParamSaveList = [flagParamSaveList;'flagParamSaveList';'fly'];
flagParamSaveList = [flagParamSaveList;'flagParamSaveList';'fly';'thisFlyBlocks';'nBlocks';'currentDate';'b']; %Note use of semicolon, not comma

% this level is the list of blocks
%for fly =[185:187]%height(flies)
for fly = chosenFlies
    %clearvars('-except', flagParamSaveList{:});
    
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

        currentSourceDirectory = fullfile(sourceDirectory,currentDate,flyID);
        if isempty(altSinkDirectory)
            currentAltSinkDirectory = [];
        else
            currentAltSinkDirectory = fullfile(altSinkDirectory,currentDate,flyID);
        end
        if length( dir([currentSourceDirectory, filesep, ['green_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']]) ) == 0 || overrideAlreadyDone == 1
            %disp(currentFlyDirectory);
            disp(['Fly: ',flyID]);

            codeStartTime = posixtime(datetime('now'));

            %loadReduceSave(currentRDMDirectory, 'green_channel.raw', currentBlock, finalSize);
            %%loadReduceSave(currentSourceDirectory, 'green_channel.raw', currentBlock, finalSize, currentAltSinkDirectory);
            if deletePrecursors == 1
                %QA
                if isempty( dir([fullfile(currentSourceDirectory,'Image_*.raw')]) ) ~= 1
                    delete(fullfile(currentSourceDirectory,'green_channel.raw'))
                    disp(['Precursor deleted'])
                else
                    ['-# Alert: Raw image not detected; Unsafe to delete precursor #-']
                end
            end

            codeEndTime = posixtime(datetime('now'));
            MET = codeEndTime - codeStartTime;
            disp(['-- Total time to process: ',num2str(MET),'s --'])
        else
            ['-# Block already processed #-']
        end
        
        clearvars('-except', flagParamSaveList{:});
    end
end

%function loadReduceSave(RDMDirectory, file, fly, finalSize, altSinkDirectory)
function loadReduceSave(currentSourceDirectory, file, fly, finalSize, currentAltSinkDirectory)

    fileRDM = fullfile(currentSourceDirectory, file);
    
    imageSize = [fly.pixelX fly.pixelY];
    nFrames = fly.realFrames;

    tic
    disp('Loading data');
    % load green channel
    fid = fopen(fileRDM, 'r','b');
    data = fread(fid, 512*512*nFrames, 'uint16'); %Hardcoded?
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
    
    if isempty(currentAltSinkDirectory)
        reducedFileDirectory = currentSourceDirectory;
    else
        reducedFileDirectory = currentAltSinkDirectory;
    end
    %Check if folder existings
    if exist(reducedFileDirectory) == 0
        mkdir(reducedFileDirectory)
        disp(['-# Making folder for reduced data output #-'])
    end
    reducedFileSink = fullfile(reducedFileDirectory, ['green_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']);
    
    tic
    disp('Saving');
    % save green channel; do not compress we care about speed not size
    save(reducedFileSink, 'rData','-v7.3','-nocompression');
    disp(['Data saved to ',reducedFileDirectory])
    toc
    
    % free up the memory 
    clear;
    
end