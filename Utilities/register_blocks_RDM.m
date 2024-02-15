% register blocks of SE trials within and across blocks
% currently does the Z averaged stacks but can be extended to align the
% original green/red channels (i.e. align each slice in the Z stack)

close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

% scratchDirectory = '../../2P Data';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
imageSize = [128 128];

chosenFlies = 85:120;

% leave empty if aligning all blocks for one fly
chosenBlocks = [];

%%

% this level is flies just in case we 
for fly = chosenFlies

    thisFlyBlocks = blocks(blocks.Fly == fly,:);

    if ~isempty(chosenBlocks)
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks),:);
    end
    
    %align inside each block
    for b = 1:height(thisFlyBlocks)
        alignBlock(thisFlyBlocks(b,:), imageSize, mainDirectory);
    end

    %no need to align across blocks if only one block
%         if height(thisFlyBlocks) > 1
%             alignAcrossBlocks(thisFlyBlocks, scratchDirectory);
%         end

end
    
function alignBlock(block, imageSize, mainDirectory)

    % slices within each volume including flyback 
    nSlices = block.Steps + block.FlybackFrames;
    
    %total number of volumes recorded
    nVolTotal = block.realFrames/nSlices;
    
    % number of volumes recorded after each train of stimuli
    nVol = nVolTotal/(block.BlockLength + block.BlankBlocks);

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
    
    % load red and green channels
    disp('Loading green channel');
    tic; green_channel = load(fullfile(currentDirectory,'green_channel_128x128')); toc;

%     disp('Loading red channel');
%     tic; red_channel = load(fullfile(currentDirectory,'red_channel_128x128')); toc;

    %hyperstack the green and red channels (pixelX,pixelY,nSlices,time)
    green_channel = reshape(green_channel.rData,[imageSize nSlices nVolTotal]);
%     rc_hstack = reshape(red_channel.rData,[imageSize nSlices nVolTotal]);
    
    % average over the volume
    avg_z_green = squeeze(mean(green_channel,3));
%     avg_z_red = squeeze(sum(rc_hstack,3));

    if block.Align

        % z-average aligned
        avg_z_green_aligned = zeros(size(avg_z_green));
    %     avg_z_red_aligned = zeros(size(avg_z_red));

        % full stack aligned
        green_channel_aligned = zeros(size(green_channel));
    
        %make a reference image for registering (mean of first recording of nVol)
        refImage = mean(avg_z_green(:,:,1:nVol),3);

        [opt,metric]=imregconfig('multimodal');

    %     opt.MaximumIterations = 300;
        opt.InitialRadius = 1e-3;

        %register green channel
        disp('Aligning stacks');
        tic;
        parfor vol = 1:nVolTotal
            
            im_trans = imregtform(avg_z_green(:,:,vol),refImage,'translation',opt,metric);
            R = imref2d(size(refImage));
            
            % apply transformation to avg image
            avg_z_green_aligned(:,:,vol) = imwarp(avg_z_green(:,:,vol),im_trans,'OutputView',R, 'SmoothEdges', false,'interp','nearest'); %#ok<*PFOUS>
    %         avg_z_red_aligned(:,:,i) = imwarp(avg_z_red(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true);
    
            % apply transformation to each slice in z direction (can this be done all at once for a volume?)
            for z = 1:nSlices
                green_channel_aligned(:,:,z,vol) = imwarp(green_channel(:,:,z,vol),im_trans,'OutputView',R, 'SmoothEdges', false,'interp','nearest'); %#ok<*PFOUS>
            end
    
        end
        toc;
        
        disp('Saving AVG green channel aligned');
        tic; save(fullfile(currentDirectory,'avg_z_green_aligned'),'avg_z_green_aligned','-v7.3','-nocompression'); toc;
        disp('Saving full green channel aligned');
        tic; save(fullfile(currentDirectory,'green_channel_aligned'),'green_channel_aligned','-v7.3','-nocompression'); toc;
    
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