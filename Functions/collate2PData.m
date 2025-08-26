function FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory, options) %New
%collate2PData Summary of this function goes here
%   Detailed explanation goes here

%{{
%Function form
arguments
    flyRecord table
    chosenFlies double
    chosenBlocks cell
    gridSize double
    dataDirectory string
    sequenceDirectory string
    options.separateByState double = 0
    options.doRolling double = 0
    options.alternateUseCase double = [] %Whether to act normally (empty), use averaged data (0), use unaligned data (1), or use 1 plane of 3D data (3),
    options.reqZ = {} %When alternateUseCase 3 (Single Z-plane), specifies which Z plane/s to use (Note that flyback frames are included here, so be careful)
end
%}
%{
%Non-function form
options.separateByState = 0;
options.doRolling = 0;
options.alternateUseCase = 3;
options.reqZ = chosenZ;
%}


% structure with necessary info
FLIES = struct;

%Arguments
separateByState = options.separateByState;
doRolling = options.doRolling;
alternateUseCase = options.alternateUseCase;
reqZ = options.reqZ;

%Pre-check for some circumstances
if alternateUseCase == 3 && isempty(reqZ)
    ['## Alert: Singular Z requested, but Z planes of interest apparently empty! ##']
    crash = yes
end

for fly = 1:length(chosenFlies)
    
    BLOCKS = struct;
    
    % the blocks corresponding to this fly
    thisFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(fly),:);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    if isempty(chosenBlocks)
        blockNumbers = thisFlyBlocks.Block.';
    else
        blockNumbers = chosenBlocks{fly};
        disp(['Using manual block specification (',num2str(blockNumbers),' for this fly)'])
    end
    
    for b = blockNumbers

        %currentBlock = thisFlyBlocks(blockNumbers==b,:)
        currentBlock = thisFlyBlocks(thisFlyBlocks.Block==b,:); %Note: This modification may have unintended consequences
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        currentDirectory = fullfile(dataDirectory,currentDate,flyID);

        disp(flyID)
        currentBlock

        BLOCKS(b).flyNum = chosenFlies(fly);
        BLOCKS(b).flyID = flyID;
        BLOCKS(b).blockNum = b;
        BLOCKS(b).Trim = currentBlock.Trim;
        %Append four-coord trim if applicable
        if any( strcmp('TrimCoords',currentBlock.Properties.VariableNames) ) %Check if field existing
            if ~isnan(currentBlock.TrimCoords{1})
                BLOCKS(b).TrimCoords = str2num( currentBlock.TrimCoords{1} );
                disp('Trim coords acquired')
            else
                disp(['NaN/empty trim coords'])
            end
        end
        %Append phase shift if applicable
        if any( strcmp('ArbSyncVolShift',currentBlock.Properties.VariableNames) ) %Check if field existing
            if ~isnan(currentBlock.ArbSyncVolShift)
                BLOCKS(b).ArbSyncVolShift = currentBlock.ArbSyncVolShift;
                disp(['Arbitrary phase shift specified (',num2str(BLOCKS(b).ArbSyncVolShift),' vols)'])
            else
                disp(['No phase shift specified'])
            end
        end
        
        %Two rolling flags
        if doRolling == 1
            BLOCKS(b).doRolling = 1; %Whether to do rolling analysis
        else
            BLOCKS(b).doRolling = 0;
        end
        if currentBlock.BlockLength == -1 %Critically relies on this as a consistent standard
            BLOCKS(b).isRolling = 1; %Whether data is actually rolling
        else
            BLOCKS(b).isRolling = 0;
        end

        %Singular Z if applicable
        if ~isempty(reqZ) 
            if iscell(reqZ)
                thisReqZ = reqZ{fly}( find( blockNumbers == b ) ); %Note: Will behave strange if same block called twice
                disp(['Using individual value for singular Z for this fly-block (Z=',num2str(thisReqZ),')'])
            elseif numel(reqZ) == 1 && ~iscell(reqZ)
                thisReqZ = reqZ;
                disp(['Using one value for singular Z across all flies/blocks (Z=',num2str(thisReqZ),')'])
            end
        end
        
        % load 128x128 data
        disp('Loading green channel');
        %if isempty(useUnaligned) || useUnaligned == 0
        if isempty(alternateUseCase) || alternateUseCase == 0
            tic; load(fullfile(currentDirectory,'avg_z_green_aligned')); toc;
        elseif alternateUseCase == 1
            disp(['-# Using unaligned data for analysis #-'])
            tic; load(fullfile(currentDirectory,'avg_z_green_unaligned')); toc;
        elseif alternateUseCase == 3
            disp(['-# Using single Z-plane data for analysis #-'])
            tic; load(fullfile(currentDirectory,'green_channel_aligned')); toc;
        else
            ['Unspecified alternate use case']
            crash = yes
        end

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
            %if isempty(useUnaligned) || useUnaligned == 0
            if isempty(alternateUseCase) || alternateUseCase == 0 %Use aligned average
                BLOCKS(b).greenChannel = imresize3(avg_z_green_aligned,[gridSize size(avg_z_green_aligned,3)],'box');
            elseif alternateUseCase == 1 %Use unaligned
                BLOCKS(b).greenChannel = imresize3(avg_z_green_unaligned,[gridSize size(avg_z_green_unaligned,3)],'box');
            elseif alternateUseCase == 3 %Use 1 plane of aligned
                BLOCKS(b).greenChannel = imresize3( squeeze(green_channel_aligned(:,:,thisReqZ,:)),[gridSize size(green_channel_aligned,4)],'box'); %Note slightly different size call
                    %Squeeze necessary to prevent crash
                disp(['Using only requested Z-plane ',num2str(thisReqZ),' of ',num2str(size(green_channel_aligned,3)),...
                    ' (',num2str(currentBlock.Steps),' real, ',num2str(currentBlock.FlybackFrames),' flyback)'])
                BLOCKS(b).singularZ = thisReqZ;
            end
            %         BLOCKS(b).redChannel = imresize3(avg_z_red_aligned,[gridSize size(avg_z_red_aligned,3)],'box');
        end
        toc;
        
        % get sequences
        if BLOCKS(b).isRolling == 0
            if exist(fullfile(sequenceDirectory,[flyID '.mat']),'file')
                randomSequence = load(fullfile(sequenceDirectory,[flyID '.mat']),'randomSequence').';
                BLOCKS(b).randomSequence = randomSequence.randomSequence-1;
            else
                BLOCKS(b).randomSequence = csvread(fullfile(sequenceDirectory,'Data_LEDs',[flyID '.csv'])).';
            end
        else
            BLOCKS(b).randomSequence = []; %Defer sequence calcs until derived by syncMaster
        end
        
        %Get behavioural data (if requested)
        if separateByState == 1
            if exist([currentDirectory,filesep,'behavSequence.mat']) ~= 0
                load([currentDirectory,filesep,'behavSequence.mat']);
                if isfield(savStruct,'acInac')
                    behavSequence = savStruct.acInac.thisInacBinaryInterp';
                    disp(['-- Behavioural data loaded --'])
                    %QA
                    if BLOCKS(b).isRolling == 0 && size(behavSequence,2) ~= size( BLOCKS(b).randomSequence )
                            %Only try calc if sequence info not deferred
                        ['## Critical desynchronisation between behav and random sequence data! ##']
                        crash = yes
                    end
                else
                    disp(['-# Activity/Inactivity calculations not performed! #-'])
                    behavSequence = [];
                end
            else
                disp(['-# Behavioural data not found! #-'])
                behavSequence = [];
            end
            BLOCKS(b).behavSequence = behavSequence;
        end        
        
        nBadTrials = 0;
        nBadBlankTrials = 0;

        %remove bad trials and associated frames (this should be put inside a function)
        if BLOCKS(b).isRolling == 0 && ~isempty(currentBlock.removeFrames{1})
            removeFrames = eval(currentBlock.removeFrames{1});
            badTrials = eval(currentBlock.badTrials{1});
            nBadTrials = length(badTrials);

            BLOCKS(b).greenChannel(:,:,removeFrames) = [];
            
            auxRandomSequence = reshape(BLOCKS(b).randomSequence,[currentBlock.nStimuli length(BLOCKS(b).randomSequence)/currentBlock.nStimuli]);
            
            nBadBlankTrials = sum(auxRandomSequence(1,badTrials) == 5);% 5's were used for a blank trial
%             nBadNormalTrials = length(badTrials)-nBadBlankTrials;

            auxRandomSequence(:,badTrials) = [];
            
            BLOCKS(b).randomSequence = auxRandomSequence(:).';
            
            %Replicate for behavSequence, if applicable
            if separateByState == 1
                auxBehavSequence = reshape(BLOCKS(b).behavSequence,[currentBlock.nStimuli length(BLOCKS(b).behavSequence)/currentBlock.nStimuli]);
                auxBehavSequence(:,badTrials) = [];
                    %NOTE: UNTESTED
            end
            disp(['-# Bad trials have been removed #-'])            
        end
        %QA to warn user if frames would have been removed if not rolling
        if BLOCKS(b).isRolling == 1 && ~isempty(currentBlock.removeFrames{1})
           ['-# Caution: Frames requested to be removed but this is not (currently) supported for rolling #-'] 
        end
                
        % calculate number of volumes per stimulus train
        nVolTotal = size(BLOCKS(b).greenChannel,3);
        BLOCKS(b).nVol = nVolTotal/(currentBlock.BlockLength + currentBlock.BlankBlocks-nBadTrials);

%         if ~exist(fullfile(currentDirectory,'brain.jpg'),'file')
%             copyfile(fullfile(RDMDirectory,'Gcamp7s_CC',currentDate,flyID,'brain.jpg'),fullfile(currentDirectory,'brain.jpg'));
%         end

        BLOCKS(b).brainImage = imread(fullfile(currentDirectory,'brain.jpg'));
        BLOCKS(b).nStimuli = currentBlock.nStimuli;
        BLOCKS(b).blankBlocks = currentBlock.BlankBlocks-nBadBlankTrials;
        
        % plot before filtering
        figure; 
        plot(squeeze(mean(mean(BLOCKS(b).greenChannel,1),2)));
        title(['Pre SG-filtered raw data - ',strrep(flyID,'_',' ')])
        xlabel(['Frame no.'])
        ylabel(['Intensity'])
        
        % apply a savitsky-golay filter to remove larger trends in data
        BLOCKS(b).greenChannel = filterChannel(BLOCKS(b).greenChannel,3,55);
        
        % plot after filtering
        figure; 
        plot(squeeze(mean(mean(BLOCKS(b).greenChannel,1),2)));
        title(['Post SG-filtered raw data - ',strrep(flyID,'_',' ')])
        xlabel(['Frame no.'])
        ylabel(['Intensity'])
    
        BLOCKS(b).blankImageStack = [];
        
        % if there are blank blocks, split image stack (after filtering)
        % this creates a baseline based on blank blocks (an F for dF/F)
        if BLOCKS(b).isRolling == 0 && currentBlock.BlankBlocks
            BLOCKS(b) = splitStack(BLOCKS(b));
                %Note: splitStack not designed for deferred sequence
        end
            
    end

    %Retroactively remove any empty blocks (e.g. discontinous block numbers)
    for b = size(BLOCKS,2):-1:1
        if isempty( BLOCKS(b).flyNum )
            BLOCKS(b) = [];
            disp(['-# Block ',num2str(b),' retroactively removed from analysis due to emptiness #-'])
        end
    end
    
    FLIES(fly).BLOCKS = BLOCKS;
    
    % concatenate data from different blocks
%     if groupedBlocks   
%         %concatenate aligned blocks
%         FLIES(fly).greenChannel = cell(1,1,nBlocks); %#ok<*UNRCH>
%     %     FLIES(fly).redChannel = cell(1,1,nBlocks);
%         FLIES(fly).randomSequence = cell(1,nBlocks);
% 
%     %     for b = [1 3]
%         for b = 1:nBlocks
%             FLIES(fly).greenChannel{b} = FLIES(fly).BLOCKS(b).greenChannel;
%     %         FLIES(fly).redChannel{b} = FLIES(fly).BLOCKS(b).redChannel;
%             FLIES(fly).randomSequence{b} = FLIES(fly).BLOCKS(b).randomSequence;
%         end
% 
%         FLIES(fly).greenChannel = cell2mat(FLIES(fly).greenChannel);
%     %     FLIES(fly).redChannel = cell2mat(FLIES(fly).redChannel);
%         FLIES(fly).randomSequence = cell2mat(FLIES(fly).randomSequence);
%         FLIES(fly).nVol = BLOCKS(1).nVol;
% %         FLIES(fly).brainImage = BLOCKS(1).brainImage; % all blocks aligned to first so brain image of first block is image for all
%     end
    
end

end