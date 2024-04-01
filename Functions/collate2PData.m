function FLIES = collate2PData(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, groupedBlocks)
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
        currentDirectory = fullfile(dataDirectory,currentDate,flyID);
        
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
        
        nBadTrials = 0;
        nBadBlankTrials = 0;

        %remove bad trials and associated frames (this should be put inside a function)
        if ~isempty(currentBlock.removeFrames{1})
            removeFrames = eval(currentBlock.removeFrames{1});
            badTrials = eval(currentBlock.badTrials{1});
            nBadTrials = length(badTrials);

            BLOCKS(b).greenChannel(:,:,removeFrames) = [];
            
            auxRandomSequence = reshape(BLOCKS(b).randomSequence,[currentBlock.nStimuli length(BLOCKS(b).randomSequence)/currentBlock.nStimuli]);
            
            nBadBlankTrials = sum(auxRandomSequence(1,badTrials) == 5);% 5's were used for a blank trial
%             nBadNormalTrials = length(badTrials)-nBadBlankTrials;

            auxRandomSequence(:,badTrials) = [];
            
            BLOCKS(b).randomSequence = auxRandomSequence(:).';

        end
        
        % calculate number of volumes per stimulus train
        nVolTotal = size(BLOCKS(b).greenChannel,3);
        BLOCKS(b).nVol = nVolTotal/(currentBlock.BlockLength+currentBlock.BlankBlocks-nBadTrials);

%         if ~exist(fullfile(currentDirectory,'brain.jpg'),'file')
%             copyfile(fullfile(RDMDirectory,'Gcamp7s_CC',currentDate,flyID,'brain.jpg'),fullfile(currentDirectory,'brain.jpg'));
%         end

        BLOCKS(b).brainImage = imread(fullfile(currentDirectory,'brain.jpg'));
        BLOCKS(b).nStimuli = currentBlock.nStimuli;
        BLOCKS(b).blankBlocks = currentBlock.BlankBlocks-nBadBlankTrials;
        
        % plot before fitlering
        figure; plot(squeeze(mean(mean(BLOCKS(b).greenChannel,1),2)));
        
        % apply a savitsky-golay filter to remove larger trends in data
        BLOCKS(b).greenChannel = filterChannel(BLOCKS(b).greenChannel,3,55);
        
        % plot after filtering
        figure; plot(squeeze(mean(mean(BLOCKS(b).greenChannel,1),2)));
    
        BLOCKS(b).blankImageStack = [];
        
        % if there are blank blocks, split image stack (after filtering)
        % this creates a baseline based on blank blocks (an F for dF/F)
        if currentBlock.BlankBlocks
            BLOCKS(b) = splitStack(BLOCKS(b));
        end
    
    end
    
    FLIES(fly).BLOCKS = BLOCKS;
    
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