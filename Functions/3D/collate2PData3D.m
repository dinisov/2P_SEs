function FLIES = collate2PData3D(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, groupedBlocks)
%collate2PData Summary of this function goes here
%   Detailed explanation goes here

% structure with necessary info
FLIES = struct;

for fly = 1:length(chosenFlies)
    
    BLOCKS = struct;
    
    % the blocks corresponding to this fly
    thisFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(fly),:);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    %nBlocks = height(thisFlyBlocks);
    blockNumbers = thisFlyBlocks.Block.';
    
    for b = blockNumbers%1:nBlocks
        
        %currentBlock = thisFlyBlocks(b,:);
        currentBlock = thisFlyBlocks(blockNumbers==b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        currentDirectory = fullfile(dataDirectory,currentDate,flyID);

        BLOCKS(b).blockNum = b;
        BLOCKS(b).Trim = currentBlock.Trim;
        
        disp(flyID);
        
        % load 128x128 data
        disp('Loading green channel');
        tic; load(fullfile(currentDirectory,'green_channel_aligned'),'green_channel_aligned'); toc;

%         disp('Loading red channel');
%         tic; load(fullfile(currentDirectory,'avg_z_red_aligned')); toc;
        
        % reduce data to desired size (usually from 128x128 to 32x32)
        disp('Reducing data');
        tic;
        % some blocks were not aligned so rData was not renamed
        % green_channel
        if exist('rData','var')
            BLOCKS(b).greenChannel = imresize3(rData,[gridSize size(rData,3)],'box');
            clear('rData');
        else
            % I am guessing I checked no interpolation is going on in z
            % direction?
            BLOCKS(b).greenChannel = imresize3(reshape(green_channel_aligned,[size(green_channel_aligned,[1 2]) prod(size(green_channel_aligned,[3 4]))]),[gridSize size(green_channel_aligned,3)*size(green_channel_aligned,4)],'box');
%         BLOCKS(b).redChannel = imresize3(avg_z_red_aligned,[gridSize size(avg_z_red_aligned,3)],'box');
        end
        toc;
        
        % reshape back into 4D array (X,Y,Z,t)
        BLOCKS(b).greenChannel = reshape(BLOCKS(b).greenChannel,[gridSize size(green_channel_aligned,3) size(green_channel_aligned,4)]);
        
        % get rid of the flyback frames
        BLOCKS(b).greenChannel = BLOCKS(b).greenChannel(:,:,1:currentBlock.Steps,:);
        
        % this is just for testing
%         BLOCKS(b).greenChannel = squeeze(mean(BLOCKS(b).greenChannel,3));
        
        % get sequences (first os for projector second for LEDs)
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

            BLOCKS(b).greenChannel(:,:,:,removeFrames) = [];
            
            auxRandomSequence = reshape(BLOCKS(b).randomSequence,[currentBlock.nStimuli length(BLOCKS(b).randomSequence)/currentBlock.nStimuli]);
            
            nBadBlankTrials = sum(auxRandomSequence(1,badTrials) == 5);% 5's were used for a blank trial
%             nBadNormalTrials = length(badTrials)-nBadBlankTrials;

            auxRandomSequence(:,badTrials) = [];
            
            BLOCKS(b).randomSequence = auxRandomSequence(:).';

        end
        
        % calculate number of volumes per stimulus train
%         nSlices = currentBlock.Steps + currentBlock.FlybackFrames;
        nVolTotal = size(BLOCKS(b).greenChannel,4);
        BLOCKS(b).nVol = nVolTotal/(currentBlock.BlockLength+currentBlock.BlankBlocks-nBadTrials);

        BLOCKS(b).brainImage = imread(fullfile(currentDirectory,'brain.jpg'));
        BLOCKS(b).nStimuli = currentBlock.nStimuli;
        BLOCKS(b).blankBlocks = currentBlock.BlankBlocks-nBadBlankTrials;
        
        % plot before fitlering
        figure; plot(squeeze(mean(mean(mean(BLOCKS(b).greenChannel,3),1),2)));
        
        % apply a savitsky-golay filter to remove larger trends in data
        BLOCKS(b).greenChannel = filterChannel3D(BLOCKS(b).greenChannel,3,55);
        
        % plot after filtering
        figure; plot(squeeze(mean(mean(mean(BLOCKS(b).greenChannel,3),1),2)));
    
        BLOCKS(b).blankImageStack = [];
        
        % if there are blank blocks, split image stack (after filtering)
        % this creates a baseline based on blank blocks (an F for dF/F)
        if currentBlock.BlankBlocks
            BLOCKS(b) = splitStack3D(BLOCKS(b));
        end
    
    end
    
    FLIES(fly).BLOCKS = BLOCKS;
    
    % concatenate data from different blocks (THIS NEEDS WORK)
    if groupedBlocks   
        %concatenate aligned blocks
        FLIES(fly).greenChannel = cell(1,1,nBlocks); %#ok<*UNRCH>
    %     FLIES(fly).redChannel = cell(1,1,nBlocks);
        FLIES(fly).randomSequence = cell(1,nBlocks);

    %     for b = [1 3]
        for b = blockNumbers%1:nBlocks
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