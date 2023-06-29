%

close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

scratchDirectory = '../../2P Data';

flies = readtable('../../2P Record/2P_record');

%get rid of excluded flies
flies = flies(~logical(flies.Exclude),:);

% the numbers here should be the original size divided by some power of 2
imageSize = [128 128];

%%

dates = unique(flies.Date);

for d = 16:16
   
    % all blocks for a particular date
    thisDateBlocks = flies(flies.Date == dates(d),:);
    
    thisDateFlies = unique(thisDateBlocks.Fly);
    
    %for each fly this date
    for i = 2:length(thisDateFlies)
        
        thisFlyBlocks = thisDateBlocks(thisDateBlocks.Fly == thisDateFlies(i),:);
        
        %align inside each block
%         for b = 1:1%height(thisFlyBlocks)
%             alignBlock(thisFlyBlocks(b,:), imageSize, scratchDirectory);
%         end
            
        %no need to align across blocks if only one block
        if height(thisFlyBlocks) > 1
            alignAcrossBlocks(thisFlyBlocks, scratchDirectory);
        end

    end
    
end
    
function alignBlock(block, imageSize, baseDirectory)

    % slices within each volume including flyback 
    nSlices = block.Steps + block.FlybackFrames;
    
    %total number of volumes recorded
    nVolTotal = block.realFrames/nSlices;
    
    % number of volumes recorded after each train of stimuli
%     nVol = nVolTotal/block.BlockLength;

    currentDate = char(datetime(block.Date,'Format','dMMMyy'));
    currentBlockDirectory = ['fly' num2str(block.FlyOnDay) '_exp' num2str(block.Block) '_' currentDate];
    currentDirectory = fullfile(baseDirectory,currentDate,currentBlockDirectory);
    disp(currentDirectory);
    % load red and green channels
    disp('Loading green channel');
    tic; green_channel = load(fullfile(currentDirectory,'green_channel_128x128')); toc;

    disp('Loading red channel');
    tic; red_channel = load(fullfile(currentDirectory,'red_channel_128x128')); toc;

    %hyperstack the green and red channels (pixelX,pixelY,nSlices,time)
    gc_hstack = reshape(green_channel.rData,[imageSize nSlices nVolTotal]);
    rc_hstack = reshape(red_channel.rData,[imageSize nSlices nVolTotal]);
    
    % average over the volume
    avg_z_green = squeeze(sum(gc_hstack,3));
    avg_z_red = squeeze(sum(rc_hstack,3));
    
    avg_z_green_aligned = zeros(size(avg_z_green));
    avg_z_red_aligned = zeros(size(avg_z_red));
    
    %make a reference image for registering (mean of first 11 time points)
    refImage = mean(avg_z_green(:,:,1:11),3);

    [opt,metric]=imregconfig('multimodal');
    
%     opt.MaximumIterations = 300;
    opt.InitialRadius = 3e-3;

    %register both green and red channel
    disp('Aligning stacks');
    tic;
    parfor i = 1:nVolTotal
        im_trans = imregtform(avg_z_green(:,:,i),refImage,'translation',opt,metric);
        R = imref2d(size(refImage));
        avg_z_green_aligned(:,:,i) = imwarp(avg_z_green(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true); %#ok<*PFOUS>
        avg_z_red_aligned(:,:,i) = imwarp(avg_z_red(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true);
    end
    toc;
    
    disp('Saving green channel');
    tic; save(fullfile(currentDirectory,'avg_z_green_aligned'),'avg_z_green_aligned'); toc;
    
    disp('Saving red channel');
    tic; save(fullfile(currentDirectory,'avg_z_red_aligned'),'avg_z_red_aligned'); toc;
    
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
        disp(['Loading avg red channel (block ' num2str(b) ')']);
        tic; load(fullfile(currentDirectory,'avg_z_red_aligned'),'avg_z_red_aligned'); toc;
        
        avg_z_green_aligned = avg_z_green_aligned.avg_z_green_aligned;

        nVolTotal = size(avg_z_green_aligned,3);

%         avg_z_green_block_aligned = zeros(size(avg_z_green_aligned));
%         avg_z_red_block_aligned = zeros(size(avg_z_red_aligned));

        [opt,metric]=imregconfig('multimodal');

        disp('Aligning stacks');
        tic;
        parfor i = 1:nVolTotal
            im_trans = imregtform(avg_z_green_aligned(:,:,i),refImage,'translation',opt,metric);
            R = imref2d(size(refImage));
            avg_z_green_aligned(:,:,i) = imwarp(avg_z_green_aligned(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true); %#ok<*PFOUS>
            avg_z_red_aligned(:,:,i) = imwarp(avg_z_red_aligned(:,:,i),im_trans,'OutputView',R, 'SmoothEdges', true);
        end
        toc;

        disp('Saving green channel');
        tic; save(fullfile(currentDirectory,'avg_z_green_aligned'),'avg_z_green_aligned'); toc;

        disp('Saving red channel');
        tic; save(fullfile(currentDirectory,'avg_z_red_aligned'),'avg_z_green_aligned'); toc;
    
    end
    
end