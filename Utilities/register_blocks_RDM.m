% register blocks of SE trials within and across blocks
% currently does the Z averaged stacks but can be extended to align the
% original green/red channels (i.e. align each slice in the Z stack)

close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

% scratchDirectory = '../../2P Data';

blocks = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");

%get rid of excluded flies
% blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
imageSize = [-1 -1]; % <value> -> Requested size, -1 -> Automatically derive size from loaded data 

chosenFlies = [324];

% leave empty if aligning all blocks for one fly
chosenBlocks = {[3]};
    %FORMAT MUST BE {[<block/s>]} 

% chosenFlies = [4 5 6 7 13 20 22 23 38 50 54];
% 
% chosenBlocks = {[1 3],1,2,[1 2],2,1,3,2,2,2,[2 3]};

altnVolSelectionMode = 2; %(Only applies for bendy data [Indicated by negative nVol in flyRecord])
    %Whether to 1 - Select the first 1% of total volume counts (e.g. 6500 volumes -> 7 volumes) or 2 - Select 1%/18 volumes (Min) equally spaced
    %Parameter for how many max frames defined in function down below because CBF adding even more arguments

overrideAlign = []; %Whether to override block Align specification (0 - No alignment, 1 - Align, Empty - No override)
    
%%

% this level is flies just in case we 
for fly = 1:length(chosenFlies)

    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);

    if ~isempty(chosenBlocks) && ~isempty(chosenBlocks{fly})
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
    end
    
    %align inside each block
    for b = 1:height(thisFlyBlocks)
        currentBlock = thisFlyBlocks(b,:);
        if all( imageSize ~= -1 ) && imageSize(1) > currentBlock.pixelY || imageSize(2) > currentBlock.pixelX
            disp(['-# One or more dimensions of image smaller than requested final size; Using initial size #-'])
            imageSizeActual = [currentBlock.pixelY, currentBlock.pixelX]
        else
            imageSizeActual = imageSize;
        end
        %alignBlock(thisFlyBlocks(b,:), imageSize, mainDirectory);
        %alignBlock(currentBlock, imageSizeActual, mainDirectory);
        alignBlock(currentBlock, imageSizeActual, mainDirectory, 'green', altnVolSelectionMode, overrideAlign);
        if currentBlock.nChannels == 2
            %alignBlock(currentBlock, imageSizeActual, mainDirectory, 'red'); %No check for existence
            alignBlock(currentBlock, imageSizeActual, mainDirectory, 'red', altnVolSelectionMode, overrideAlign);
        end
    end

    %no need to align across blocks if only one block
%         if height(thisFlyBlocks) > 1
%             alignAcrossBlocks(thisFlyBlocks, scratchDirectory);
%         end

end
    
function alignBlock(block, imageSize, mainDirectory, colour, altnVolSelectionMode, overrideAlign)

    %For non-function operations
    %{
    block = currentBlock
    imageSize = imageSizeActual
    mainDirectory
    colour = 'green'
    %}

    if ~exist('colour', 'var') || ( isempty(colour) )
        colour = 'green'; %Default
    end
    if ~exist('altnVolSelectionMode', 'var') || ( isempty(altnVolSelectionMode) )
        altnVolSelectionMode = 1; %Default
        disp(['-# nVol selection mode not specified; Using default #-'])
    end
    
    if altnVolSelectionMode == 2
        eqSpaceMaxAllowable = 128; %How many volumes max to allow for use
            %Note: Time spent registering is more a factor of image size and number of volumes than how many frames go into the mean image
    end

    if ~exist('overrideAlign', 'var') || ( isempty(overrideAlign) )
        overrideAlign = [];
    else
        disp(['Align action overriding requested (',num2str(overrideAlign),')'])
    end

    % slices within each volume including flyback 
    nSlices = block.Steps + block.FlybackFrames;
    %QA
    if block.FlybackFrames == 99
        ['## Alert: Extreme likelihood of standin flybackFrames value in flyRecord ##']
        crash = yes
    end
    
    %total number of volumes recorded
    nVolTotal = block.realFrames/nSlices;
    
    % number of volumes recorded after each train of stimuli
    nVol = nVolTotal/(block.BlockLength + block.BlankBlocks);
    if nVol < 0
        disp(['-# Probable rolling nature detected in nVol #-'])
    end

    currentDate = char(datetime(block.Date,'Format','dMMMyy'));
    currentBlockDirectory = ['fly' num2str(block.FlyOnDay) '_exp' num2str(block.Block) '_' currentDate];
    currentDirectory = fullfile(mainDirectory,currentDate,currentBlockDirectory);
    disp(currentDirectory);
    
%     disp('Copying original file');
%     tic;
%     if ~exist(fullfile(currentDirectory,'green_channel_128x128'),'file')
%         copyfile(fullfile(mainDirectory,currentDate,currentBlockDirectory,'green_channel_128x128.mat'),fullfile(currentDirectory,'green_channel_128x128.mat'));
%     end
%     toc;

    %thisFile = dir( [fullfile(currentDirectory,'green_channel_*x*.mat')] );
    thisFile = dir( [fullfile(currentDirectory,[colour,'_channel_*x*.mat'])] );

    if ( isempty(overrideAlign) & block.Align == 1 ) || overrideAlign == 1 %Might not work fully
        oldFile = ['avg_z_',colour,'_aligned.mat'];
    else
        oldFile = ['avg_z_',colour,'_unaligned.mat'];
    end

    %if exist(fullfile(currentDirectory,'green_channel_128x128.mat'),'file') && ~exist(fullfile(currentDirectory,'avg_z_green_aligned.mat'),'file')
    %if ~isempty( thisFile ) && ~exist(fullfile(currentDirectory,'avg_z_green_aligned.mat'),'file')
    if ~isempty( thisFile ) && ~exist(fullfile(currentDirectory,oldFile),'file') %~exist(fullfile(currentDirectory,['avg_z_',colour,'_aligned.mat']),'file')
    
        % load red and green channels
        %disp('Loading green channel');
        disp(['Loading ',colour,' channel']);
        %tic; green_channel = load(fullfile(currentDirectory,'green_channel_128x128')); toc;
        tic; colour_channel = load([ thisFile.folder,filesep,thisFile.name ]); toc; %Will probs crash if >1 file; Replaces "green_channel"

    %     disp('Loading red channel');
    %     tic; red_channel = load(fullfile(currentDirectory,'red_channel_128x128')); toc;
    
        %Automatically derive image size if requested
        if any( imageSize == -1 )
            disp(['Automatically deriving image size'])
            imageSize = size( colour_channel.rData, [1,2] )
        end

        %hyperstack the green and red channels (pixelX,pixelY,nSlices,time)
        %[imageSize nSlices nVolTotal]
        colour_channel = reshape(colour_channel.rData,[imageSize nSlices nVolTotal]);
    %     rc_hstack = reshape(red_channel.rData,[imageSize nSlices nVolTotal]);
    
        disp(['Image size: ', num2str(imageSize),', nSlices: ', num2str(nSlices), ', nVolTotal: ', num2str(nVolTotal)])

        % average over the volume
        avg_z_colour = squeeze(mean(colour_channel,3)); %Replaces "avg_z_green"
    %     avg_z_red = squeeze(sum(rc_hstack,3));
    
        disp(['Averaged image dimensions: ', num2str( size(avg_z_colour) )])

        if overrideAlign == 1 || ( isempty(overrideAlign) && block.Align )

            % z-average aligned
            avg_z_colour_aligned = zeros(size(avg_z_colour)); %Replaces "avg_z_green_aligned"
        %     avg_z_red_aligned = zeros(size(avg_z_red));

            % full stack aligned
            colour_channel_aligned = zeros(size(colour_channel)); %Replaces "green_channel_aligned"

            %make a reference image for registering (mean of first recording of nVol)
            if nVol > 0
                refImage = mean(avg_z_colour(:,:,1:nVol),3); %Old block calcs
            else
                if altnVolSelectionMode == 1
                    disp(['Using first 1%/ ', num2str(ceil(size(avg_z_colour,3)*0.001)),' averaged volumes as reference'])
                    refImage = mean(avg_z_colour(:,:, 1:ceil(size(avg_z_colour,3)*0.001) ),3); %Use first 1% of total frames as reference
                        %Note: Might have issues with very short recordings, etc
                else
                    minVolCount = min( [ceil(size(avg_z_colour,3)*0.001), eqSpaceMaxAllowable] );
                    disp(['Using min of ',num2str([ceil(size(avg_z_colour,3)*0.001), eqSpaceMaxAllowable]),' (', num2str( minVolCount ),') equally spaced as ref'])
                    refVolInds = floor( linspace( 1, size(avg_z_colour,3) , minVolCount ) );
                    %QA
                    if any( refVolInds < 1 ) || any( refVolInds > size(avg_z_colour,3) ) || numel( unique(refVolInds) ) ~= numel( refVolInds )
                        ['## Alert: Either sub-zero vol, overmax vol, or non-unique vols requested for reference ##']
                        crash = yes
                    end
                    refImage = mean(avg_z_colour(:,:, refVolInds ),3);
                    
                end
            end

            [opt,metric]=imregconfig('multimodal');

        %     opt.MaximumIterations = 300;
            opt.InitialRadius = 1e-3;

            %register green channel
            disp('Aligning stacks');
            volMarkers = floor(linspace(1,nVolTotal,10)); %Used for progress reports
            tic;
            parfor vol = 1:nVolTotal

                im_trans = imregtform(avg_z_colour(:,:,vol),refImage,'translation',opt,metric);
                R = imref2d(size(refImage));

                % apply transformation to avg image
                avg_z_colour_aligned(:,:,vol) = imwarp(avg_z_colour(:,:,vol),im_trans,'OutputView',R, 'SmoothEdges', false,'interp','nearest'); %#ok<*PFOUS>
        %         avg_z_red_aligned(:,:,i) = imwarp(avg_z_red(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true);

                % apply transformation to each slice in z direction (can this be done all at once for a volume?)
                for z = 1:nSlices
                    colour_channel_aligned(:,:,z,vol) = imwarp(colour_channel(:,:,z,vol),im_trans,'OutputView',R, 'SmoothEdges', false,'interp','nearest'); %#ok<*PFOUS>
                end
                
                %Report
                if any( vol == volMarkers ) %Note: Given parfor nature, these may all hit at same time, depending on number of parallel workers
                    disp([ num2str(find( vol == volMarkers, 1, 'first' )),'% marker reached' ])                    
                end

            end
            toc;
            disp(['(',num2str(nVolTotal/toc),' vol/s)'])

            %disp('Saving AVG green channel aligned');
            disp(['Saving AVG ',colour,' channel aligned']);
            %tic; save(fullfile(currentDirectory,'avg_z_green_aligned'),'avg_z_green_aligned','-v7.3','-nocompression'); toc;
            temp = struct;
            temp.(['avg_z_',colour,'_aligned']) = avg_z_colour_aligned; %'Rename'
            tic; save(fullfile(currentDirectory,['avg_z_',colour,'_aligned']),'-struct','temp','-v7.3','-nocompression'); toc;
            %disp('Saving full green channel aligned');
            disp(['Saving full ',colour,' channel aligned']);
            %tic; save(fullfile(currentDirectory,'green_channel_aligned'),'green_channel_aligned','-v7.3','-nocompression'); toc;
            temp = struct;
            temp.([colour,'_channel_aligned']) = colour_channel_aligned; %'Rename'
            tic; save(fullfile(currentDirectory,[colour,'_channel_aligned']),'-struct', 'temp','-v7.3','-nocompression'); toc;
        
        else
            
            disp(['-# Not aligning, by request #-'])

            disp(['Saving AVG ',colour,' channel UNaligned']);
            temp = struct;
            temp.(['avg_z_',colour,'_unaligned']) = avg_z_colour; %Use original
            tic; save(fullfile(currentDirectory,['avg_z_',colour,'_unaligned']),'-struct','temp','-v7.3','-nocompression'); toc;
            disp(['Saving full ',colour,' channel UNaligned']);
            temp = struct;
            temp.([colour,'_channel_unaligned']) = colour_channel; %Ditto
            tic; save(fullfile(currentDirectory,[colour,'_channel_unaligned']),'-struct', 'temp','-v7.3','-nocompression'); toc;

        end

    %     disp('Saving green channel before alignment');
    %     tic; save(fullfile(currentDirectory,'avg_z_green'),'avg_z_green','-v7.3','-nocompression'); toc;

    %     disp('Saving red channel');
    %     tic; save(fullfile(currentDirectory,'avg_z_red_aligned'),'avg_z_red_aligned'); toc;

        %delete files
    %     disp('Deleting files');
    %     tic;
    %     if exist(fullfile(currentDirectory,'green_channel_128x128.mat'),'file')
    %         delete(fullfile(currentDirectory,'green_channel_128x128.mat'));
    %     end
    %     toc;
    
    else
        disp(['-# Either no data found, or registered data already existing #-'])
        disp(['Old data search: ',oldFile])
    end
    
end

function alignAcrossBlocks(thisFlyBlocks, baseDirectory)

    %first block should have the same date as all others
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    currentBlockDirectory = ['fly' num2str(thisFlyBlocks.FlyOnDay(1)) '_exp' num2str(thisFlyBlocks.Block(1)) '_' currentDate];
    currentDirectory = fullfile(baseDirectory,currentDate,currentBlockDirectory);
    disp(currentDirectory);
    
    % mean of first block as reference for other blocks
    disp('Loading avg green channel (block 1)');
    tic; load(fullfile(currentDirectory,'avg_z_green_aligned'),'avg_z_green_aligned'); toc;
    
    %make a reference image for registering the remaining blocks
    refImage = mean(avg_z_green_aligned,3);
    
    clear avg_z_green_aligned
    
    for b = 2:height(thisFlyBlocks)
        
        currentBlockDirectory = ['fly' num2str(thisFlyBlocks.FlyOnDay(b)) '_exp' num2str(thisFlyBlocks.Block(b)) '_' currentDate];
        currentDirectory = fullfile(baseDirectory,currentDate,currentBlockDirectory);
        disp(currentDirectory);

        disp(['Loading avg green channel (block ' num2str(b) ')']);
        tic; avg_z_green_aligned = load(fullfile(currentDirectory,'avg_z_green_aligned'),'avg_z_green_aligned'); toc;
%         disp(['Loading avg red channel (block ' num2str(b) ')']);
%         tic; load(fullfile(currentDirectory,'avg_z_red_aligned'),'avg_z_red_aligned'); toc;
        
        avg_z_green_aligned = avg_z_green_aligned.avg_z_green_aligned;

        nVolTotal = size(avg_z_green_aligned,3);

%         avg_z_green_block_aligned = zeros(size(avg_z_green_aligned));
%         avg_z_red_block_aligned = zeros(size(avg_z_red_aligned));

        [opt,metric]=imregconfig('multimodal');
        opt.InitialRadius = 1e-3;
%         opt.MaximumIterations = 300;

        disp('Aligning stacks');
        tic;
        parfor i = 1:nVolTotal
            im_trans = imregtform(avg_z_green_aligned(:,:,i),refImage,'translation',opt,metric);
            R = imref2d(size(refImage));
            avg_z_green_aligned(:,:,i) = imwarp(avg_z_green_aligned(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true); %#ok<*PFOUS>
%             avg_z_red_aligned(:,:,i) = imwarp(avg_z_red_aligned(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true);
        end
        toc;

        disp('Saving green channel');
        tic; save(fullfile(currentDirectory,'avg_z_green_aligned'),'avg_z_green_aligned'); toc;

%         disp('Saving red channel');
%         tic; save(fullfile(currentDirectory,'avg_z_red_aligned'),'avg_z_green_aligned'); toc;
    
    end
    
end