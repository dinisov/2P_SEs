%

close all; clear;

rdmDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';
blocks = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");
%rdmDirectory = 'I:\RFDG2021-Q4413\Andre\2p_Data\';
%blocks = readtable("I:\RFDG2021-Q4413\Andre\2p_Record\Andre_2P_record.xlsx");
%['Andre data']

%get rid of excluded flies
% blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
finalSize = [128 128];

%%

% chosenFlies = [4 5 6 7 13 20 22 23 38 50 54];
% chosenBlocks = {[1 3],1,2,[1 2],2,1,3,2,2,2,[2 3]};

chosenFlies = [357];
chosenBlocks = {[1]}; % leave empty if reducing all blocks for one fly
    %MUST BE IN FORMAT {[blocks]}
colourSkip = [0]; %Which colour/s to not skip (1-Green, 2-Red)
    % 0 OR [1,2] - Do all, [1] - Do green, [2] - Do red

flagParamSaveList = who;
flagParamSaveList = [flagParamSaveList;'flagParamSaveList';'fly'];

% this level is the list of blocks
%for fly =[185:187]%height(flies)
for fly = 1:length(chosenFlies)
    clearvars('-except', flagParamSaveList{:});
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);

    if ~isempty(chosenBlocks) && ~isempty(chosenBlocks{fly})
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
    end

    nBlocks = height(thisFlyBlocks);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    %currentFlyDirectory = ['fly' num2str(currentFly.FlyOnDay) '_exp' num2str(currentFly.Block) '_' currentDate];
    %for b = 1:nBlocks
    for b = chosenBlocks{fly}
        %currentBlock = thisFlyBlocks(b,:);
        currentBlock = thisFlyBlocks( find( thisFlyBlocks.Block == b ) ,:); %New
        
        %currentBlock
        
        %QA
        if isempty( currentBlock )
            ['## No blocks found for fly #', num2str(chosenFlies(fly)),' block ', num2str(b)]
            crash = yes
        end
        if finalSize(1) > currentBlock.pixelY || finalSize(2) > currentBlock.pixelX
            disp(['-# One or more dimensions of image smaller than requested final size; Using initial size #-'])
            finalSizeActual = [currentBlock.pixelY, currentBlock.pixelX]
        else
            finalSizeActual = finalSize;
        end
        if currentBlock.pixelY ~= currentBlock.pixelX
            ['-# reduction of files not tested on non-square data #-']
            crash = yes
            %Need to confirm order of dims
        end
        
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate]; %Borrowed from pre_process

        currentRDMDirectory = fullfile(rdmDirectory,currentDate,flyID);

        %disp(currentFlyDirectory);
        disp(['Fly: ',flyID]);
        
        codeStartTime = posixtime(datetime('now'));
        if ismember(colourSkip,0) || ismember(colourSkip,1)
            %loadReduceSave(currentRDMDirectory, 'green_channel.raw', currentBlock, finalSize, 1); %currentBlock corresponds to currentFly
            loadReduceSave(currentRDMDirectory, 'green_channel.raw', currentBlock, finalSizeActual, 1, 'green'); %currentBlock corresponds to currentFly
        else
            disp(['-# Skipping green reduction by request #-'])
        end
        if currentBlock.nChannels == 2
            temp = dir( [currentRDMDirectory,filesep,'red_channel.raw'] );
            if ~isempty( temp )
                if ismember(colourSkip,0) || ismember(colourSkip,2)
                    disp(['Reducing and saving red channel data as well'])
                    loadReduceSave(currentRDMDirectory, 'red_channel.raw', currentBlock, finalSizeActual, 1, 'red');      
                else
                    disp(['-# Skipping red reduction by request #-'])
                end
            else
                ['-# Alert: Two channels specified, but red channel data not found #-']
            end
        end
        
        
        codeEndTime = posixtime(datetime('now'));
        MET = codeEndTime - codeStartTime;
        disp(['-- Total time to process: ',num2str(MET),'s --']) 
    end
end

function loadReduceSave(RDMDirectory, file, fly, finalSize,fragments, colour)
    if ~exist('colour', 'var') || ( isempty(colour) )
        colour = 'green'; %Default
    end
        
    %memUsed = nan(7,1);
     if ~isunix
        [memStruct,~] = memory;
        %memUsed(1) = memStruct.MemUsedMATLAB;
        disp(['Stage 1 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
    else
        disp(['Stage 1 (Mem. info unavailable)'])
    end

    fileRDM = fullfile(RDMDirectory, file);
    disp(['Now reading: ',fileRDM])

    imageSize = [fly.pixelX fly.pixelY];
    nFrames = fly.realFrames;

    if fragments == 1
        %Original
        tic
        disp('Loading data');
        % load green channel
        fid = fopen(fileRDM, 'r','b');
        %data = fread(fid, 512*512*nFrames, 'uint16');
        data = fread(fid, fly.pixelX*fly.pixelY*nFrames, 'uint16');
        fclose(fid);
        if ~isunix
            [memStruct,~] = memory;
            %memUsed(2) = memStruct.MemUsedMATLAB;
            disp(['Stage 2 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        else
            disp(['Stage 2 (Mem. info unavailable)'])
        end
        toc
        %[imageSize nFrames]
        %size(data)
        %rearrange
        %data = permute(reshape(data, [imageSize nFrames]),[2 1 3]);
        try
            data = permute(reshape(data, [imageSize nFrames]),[2 1 3]);
        catch ME
            if (strcmp(ME.identifier,'MATLAB:getReshapeDims:notSameNumel'))
                ['## Alert: Error in resizing raw data ##']
                ['Raw data size: ',num2str(size(data))]
                ['Requested new dims: ',num2str([imageSize nFrames])]
                crash = yes
            else
                ['## Unexpected error in resizing ##']
                crash = yes
            end
        end
        if ~isunix
            [memStruct,~] = memory;
            disp(['Stage 3 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        else
            disp(['Stage 3 (Mem. info unavailable)'])
        end
    
        %%rData = zeros([finalSize nFrames]); %Disabled, since imresize3 almost certainly makes a new variable?
        
        tic
        disp('Reducing size');
        %reduce green channel (nothing comes close to this in terms of speed)
        %rData = imresize3(data,size(rData),'box');
        rData = imresize3(data,[finalSize nFrames],'box'); %Adjusted since rData no longer 'pre' defined
        if ~isunix
            [memStruct,~] = memory;
            disp(['Stage 4 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        else
            disp(['Stage 4 (Mem. info unavailable)'])
        end
        clear data
        if ~isunix
            [memStruct,~] = memory;
            disp(['Stage 5 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        else
            disp(['Stage 5 (Mem. info unavailable)'])
        end
        toc
    else
        %Memory-efficient
        disp('Loading data fragmentised');
        nFramesEffective = nFrames/fragments;
        disp(['(',num2str(nFramesEffective),' frames at a time)'])
        rData = nan([finalSize nFrames]); %Pre-allocate
            %Theoretically could be done in a rolling manner, but might be
            %too much effort...
        disp(['rData preallocated'])
        %QA
        if floor(nFramesEffective) ~= nFramesEffective 
            ['-# Alert: Data cannot be integerly split into ',num2str(fragments),' pieces #-']
            crash = yes
                %Could relatively easily write an adjustment so that frames
                %are properly loaded though...
        end
        a = 1;
        fid = fopen(fileRDM, 'r','b'); %As below, need to only do this once
        for frag = 1:fragments
            
            tic
            % load green channel
            %%fid = fopen(fileRDM, 'r','b');
            %%data = fread(fid, 512*512*nFrames, 'uint16');
            [data,nFramesRead] = fread(fid, 512*512*nFramesEffective, 'uint16'); %Iteratively load fractions at a time
            %fclose(fid); %Need to not do this to keep place in file            
            
            if ~isunix
                [memStruct,~] = memory;
                disp(['Fragment ',num2str(frag),' Stage 2 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
            else
                disp(['Fragment ',num2str(frag),' Stage 2 (Mem. info unavailable'])
            end
            
            toc
            %QA for correct number of frames read
            if nFramesRead ~= 512*512*nFramesEffective
                ['-# Alert: Potentially incorrect number of frames read from file #-']
                crash = yes
                %Most likely explanation here is that file contained fewer
                %frames than advertised
            end
            
            %rearrange
            %data = permute(reshape(data, [imageSize nFrames]),[2 1 3]);
            data = permute(reshape(data, [imageSize nFramesEffective]),[2 1 3]);
            
            if ~isunix
                [memStruct,~] = memory;
                disp(['Fragment ',num2str(frag),' Stage 3 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
            else
                disp(['Fragment ',num2str(frag),' Stage 3 (Mem. info unavailable)'])
            end
        
            %%rData = zeros([finalSize nFrames]); %Disabled, since imresize3 almost certainly makes a new variable?
            
            tic
            disp('Reducing size');
            %reduce green channel (nothing comes close to this in terms of speed)
            %rData = imresize3(data,size(rData),'box');
            %rData = imresize3(data,[finalSize nFrames],'box'); %Adjusted since rData no longer 'pre' defined
            rData(: , : , a:a+nFramesEffective-1) = imresize3(data,[finalSize nFramesEffective],'box');
            
            if ~isunix
                [memStruct,~] = memory;
                disp(['Fragment ',num2str(frag),' Stage 4 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
            else
                disp(['Fragment ',num2str(frag),' Stage 4 (Mem. info unavailable)'])
            end

            clear data
            
            if ~isunix
                [memStruct,~] = memory;
                disp(['Fragment ',num2str(frag),' Stage 5 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
            else
                disp(['Fragment ',num2str(frag),' Stage 5 (Mem. info unavailable)'])
            end

            toc
            a = a + nFramesEffective;
        end
        fclose(fid);
    end
    
    %reducedFileRDM = fullfile(RDMDirectory, ['green_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']);
    reducedFileRDM = fullfile(RDMDirectory, [colour,'_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']);
    
    tic
    disp('Saving');
    % save green channel; do not compress we care about speed not size
    save(reducedFileRDM, 'rData','-v7.3','-nocompression');

    if ~isunix
        [memStruct,~] = memory;
        if fragments ~= 1
            disp(['Fragment ',num2str(frag),' Stage 6 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        else
            disp(['Stage 6 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        end
    else
        if fragments ~= 1
            disp(['Fragment ',num2str(frag),' Stage 6 (Mem. info unavailable)'])
        else
            disp(['Stage 6 (Mem. info unavailable)'])
        end
    end

    toc
    
    % free up the memory 
    clear;

    if ~isunix
        [memStruct,~] = memory;
        disp(['Stage 7 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
    else
        disp(['Stage 7 (Mem. info unavailable)'])
    end
    
end