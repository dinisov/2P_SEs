%

close all; clear;

rdmDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");

%get rid of excluded flies
% blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
finalSize = [128 128];

%%

% chosenFlies = [4 5 6 7 13 20 22 23 38 50 54];
% chosenBlocks = {[1 3],1,2,[1 2],2,1,3,2,2,2,[2 3]};

chosenFlies = [228];
chosenBlocks = {[3]}; % leave empty if reducing all blocks for one fly
    %MUST BE IN FORMAT {[blocks]}

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
    for b = 1:nBlocks
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate]; %Borrowed from pre_process

        currentRDMDirectory = fullfile(rdmDirectory,currentDate,flyID);

        %disp(currentFlyDirectory);
        disp(['Fly: ',flyID]);
        
        codeStartTime = posixtime(datetime('now'));
        loadReduceSave(currentRDMDirectory, 'green_channel.raw', currentBlock, finalSize, 3); %currentBlock corresponds to currentFly
        
        codeEndTime = posixtime(datetime('now'));
        MET = codeEndTime - codeStartTime;
        disp(['-- Total time to process: ',num2str(MET),'s --']) 
    end
end

function loadReduceSave(RDMDirectory, file, fly, finalSize,fragments)
    %New argument 'fragments' added, which dictates how many pieces to split loading/etc into
        %For context, a standard 1h recording will require ~200GB of RAM if not fragmentised
        
    %memUsed = nan(7,1);
    [memStruct,~] = memory;
    %memUsed(1) = memStruct.MemUsedMATLAB;
    disp(['Stage 1 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])

    fileRDM = fullfile(RDMDirectory, file);

    imageSize = [fly.pixelX fly.pixelY];
    nFrames = fly.realFrames;

    if fragments == 1
        %Original
        tic
        disp('Loading data');
        % load green channel
        fid = fopen(fileRDM, 'r','b');
        data = fread(fid, 512*512*nFrames, 'uint16');
        fclose(fid);
        [memStruct,~] = memory;
        %memUsed(2) = memStruct.MemUsedMATLAB;
        disp(['Stage 2 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        toc
        
        %rearrange
        data = permute(reshape(data, [imageSize nFrames]),[2 1 3]);
        [memStruct,~] = memory;
        %memUsed(3) = memStruct.MemUsedMATLAB;
        disp(['Stage 3 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
    
        %%rData = zeros([finalSize nFrames]); %Disabled, since imresize3 almost certainly makes a new variable?
        
        tic
        disp('Reducing size');
        %reduce green channel (nothing comes close to this in terms of speed)
        %rData = imresize3(data,size(rData),'box');
        rData = imresize3(data,[finalSize nFrames],'box'); %Adjusted since rData no longer 'pre' defined
        [memStruct,~] = memory;
        %memUsed(4) = memStruct.MemUsedMATLAB;
        disp(['Stage 4 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        clear data
        [memStruct,~] = memory;
        %memUsed(5) = memStruct.MemUsedMATLAB;
        disp(['Stage 5 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
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
            [memStruct,~] = memory;
            %memUsed(2) = memStruct.MemUsedMATLAB;
            disp(['Fragment ',num2str(frag),' Stage 2 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
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
            [memStruct,~] = memory;
            %memUsed(3) = memStruct.MemUsedMATLAB;
            disp(['Stage 3 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
        
            %%rData = zeros([finalSize nFrames]); %Disabled, since imresize3 almost certainly makes a new variable?
            
            tic
            disp('Reducing size');
            %reduce green channel (nothing comes close to this in terms of speed)
            %rData = imresize3(data,size(rData),'box');
            %rData = imresize3(data,[finalSize nFrames],'box'); %Adjusted since rData no longer 'pre' defined
            rData(: , : , a:a+nFramesEffective-1) = imresize3(data,[finalSize nFramesEffective],'box');
            [memStruct,~] = memory;
            %memUsed(4) = memStruct.MemUsedMATLAB;
            disp(['Stage 4 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
            clear data
            [memStruct,~] = memory;
            %memUsed(5) = memStruct.MemUsedMATLAB;
            disp(['Stage 5 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
            toc
            a = a + nFramesEffective;
        end
        fclose(fid);
    end
    
    reducedFileRDM = fullfile(RDMDirectory, ['green_channel_' num2str(finalSize(1)) 'x' num2str(finalSize(2)) '.mat']);
    
    tic
    disp('Saving');
    % save green channel; do not compress we care about speed not size
    save(reducedFileRDM, 'rData','-v7.3','-nocompression');
    [memStruct,~] = memory;
    %memUsed(6) = memStruct.MemUsedMATLAB;
    disp(['Stage 6 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
    toc
    
    % free up the memory 
    clear;
    [memStruct,~] = memory;
    %memUsed(7) = memStruct.MemUsedMATLAB;
    disp(['Stage 7 mem. used: ',num2str(memStruct.MemUsedMATLAB/1000/10000)])
    
end