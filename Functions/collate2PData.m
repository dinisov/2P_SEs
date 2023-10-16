function FLIES = collate2PData(flyRecord, chosenFlies, gridSize, mainDirectory, RDMDirectory, sequenceDirectory, groupedBlocks)
%collate2PData Summary of this function goes here
%   Detailed explanation goes here

% structure with necessary info
FLIES = struct;

for fly = 1:length(chosenFlies)
    
    BLOCKS = struct;
    
    % the blocks corresponding to this fly
    thisFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(fly),:);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
        
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        currentDirectory = fullfile(mainDirectory,currentDate,flyID);
        
        disp(flyID);
        
        % load 128x128 data
        disp('Loading green channel');
        tic; load(fullfile(currentDirectory,'avg_z_green_aligned')); toc;

%         disp('Loading red channel');
%         tic; load(fullfile(currentDirectory,'avg_z_red_aligned')); toc;
        
        % reduce data to desired size
        disp('Reducing data');
        tic;
        % some blocks were not aligned so rData was not renamed
        % green_channel
        if exist('rData','var')
            BLOCKS(b).greenChannel = imresize3(rData,[gridSize size(rData,3)],'box');
            clear('rData');
        else
            BLOCKS(b).greenChannel = imresize3(avg_z_green_aligned,[gridSize size(avg_z_green_aligned,3)],'box');
%         BLOCKS(b).redChannel = imresize3(avg_z_red_aligned,[gridSize size(avg_z_red_aligned,3)],'box');
        end
        toc;
        
        % get sequences
        if exist(fullfile(sequenceDirectory,[flyID '.mat']),'file')
            randomSequence = load(fullfile(sequenceDirectory,[flyID '.mat']),'randomSequence').';
            BLOCKS(b).randomSequence = randomSequence.randomSequence-1;
        else
            BLOCKS(b).randomSequence = csvread(fullfile(sequenceDirectory,'Data_LEDs',[flyID '.csv'])).';
        end
        
        % add block length for convenience
        BLOCKS(b).length = size(BLOCKS(b).greenChannel,3);
        
        % calculate number of volumes per stimulus train
        nSlices = currentBlock.Steps + currentBlock.FlybackFrames;
        nVolTotal = currentBlock.realFrames/nSlices;
        BLOCKS(b).nVol = nVolTotal/currentBlock.BlockLength;

        if ~exist(fullfile(currentDirectory,'brain.jpg'),'file')
            copyfile(fullfile(RDMDirectory,'Gcamp7s_CC',currentDate,flyID,'brain.jpg'),fullfile(currentDirectory,'brain.jpg'));
        end

        BLOCKS(b).brainImage = imread(fullfile(currentDirectory,'brain.jpg'));
        BLOCKS(b).nStimuli = currentBlock.nStimuli;
        
    end
        
    FLIES(fly).BLOCKS = BLOCKS;
    
    %filter parameters
    degrees = 3;
    filt_width = 55;
    
    % filter data (remove larger trends in time series)
    for b = 1:length(BLOCKS)
        traceGreen = squeeze(mean(mean(FLIES(fly).BLOCKS(b).greenChannel,1),2));
%         traceRed = squeeze(mean(mean(FLIES(fly).BLOCKS(b).redChannel,1),2));
        
        BLOCKS(b).greenChannel = BLOCKS(b).greenChannel - reshape(sgolayfilt(traceGreen,degrees,filt_width), [1 1 length(traceGreen)]);
%         BLOCKS(b).redChannel = BLOCKS(b).redChannel - reshape(sgolayfilt(traceRed,degrees,filt_width), [1 1 length(traceRed)]);
    end
    
    % concatenate data from different blocks
    if groupedBlocks   
        %concatenate aligned blocks
        FLIES(fly).greenChannel = cell(1,1,nBlocks); %#ok<*UNRCH>
    %     FLIES(fly).redChannel = cell(1,1,nBlocks);
        FLIES(fly).randomSequence = cell(1,nBlocks);

    %     for b = [1 3]
        for b = 1:nBlocks
            FLIES(fly).greenChannel{b} = FLIES(fly).BLOCKS(b).greenChannel;
    %         FLIES(fly).redChannel{b} = FLIES(fly).BLOCKS(b).redChannel;
            FLIES(fly).randomSequence{b} = FLIES(fly).BLOCKS(b).randomSequence;
        end

        FLIES(fly).greenChannel = cell2mat(FLIES(fly).greenChannel);
    %     FLIES(fly).redChannel = cell2mat(FLIES(fly).redChannel);
        FLIES(fly).randomSequence = cell2mat(FLIES(fly).randomSequence);
        FLIES(fly).nVol = BLOCKS(1).nVol;
%         FLIES(fly).brainImage = BLOCKS(1).brainImage; % all blocks aligned to first so brain image of first block is image for all
    end
    
end

end