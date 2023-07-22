close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_swinderen\Dinis\2P SEs\Functions\');

close all; clear;

% mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';
mainDirectory = '../2P Data';

%where the sequence data is located (stimulus files)
% sequenceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\SE_2P_data\Data_LEDs';
sequenceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\SE_2P_data';

% scratchDirectory = '../../2P Data';

blocks = readtable('../2P Record/2P_record');

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [32 32];

flyList = unique(blocks.Fly);

chosenFlies = 30:35;

%% collate, reduce, filter and concatenate pre-aligned data

% structure with necessary info
FLIES = struct;

for fly = 1:length(chosenFlies)
    
    BLOCKS = struct;
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
        
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        currentDirectory = fullfile(mainDirectory,currentDate,flyID);
        
        disp(flyID);
        
        % add brain image (for now just the Z project of 1st block)
        if b == 1
            FLIES(fly).brainImage = imread(fullfile(currentDirectory,'brain.jpg'));
        end
        
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
    
    %concatenate aligned blocks
    FLIES(fly).greenChannel = cell(1,1,nBlocks);
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
    
end

%% analyse SEs
% separates images according to preceding sequence of stimuli and
% calculates mean images as a function of the sequence

%structure to put results in
R = struct;

disp('Analysing SEs...');
tic;
for fly = 1:length(FLIES)
    thisFly = FLIES(fly);
    for b = 1:length(thisFly.BLOCKS)
        thisBlock = thisFly.BLOCKS(b);
        R(fly).BLOCK(b) = analyseSequentialEffectsTwoPhoton2(thisBlock.greenChannel,thisBlock.randomSequence,thisBlock.nVol);
    end
    R(fly).ALL = analyseSequentialEffectsTwoPhoton2(thisFly.greenChannel,thisFly.randomSequence,thisFly.nVol);
end
toc;

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS per volume
% load slrp_lrpr.mat

load six_hertz.mat

disp('Calculating R and R^2 per volume');
tic;
% all blocks for each fly (separation just for organisation purposes)
for fly = 1:length(FLIES)
    meanDataSeq = R(fly).ALL.meanDataSeq;
    for vol = 1:FLIES(fly).nVol
        thisVolData = permute(squeeze(meanDataSeq(vol,:,:,:)),[2,3,1]);
        R(fly).ALL.r2Vol(vol) = calculateR2(thisVolData,six_hertz);
        R(fly).ALL.rVol(vol) = calculateR(thisVolData,six_hertz);
    end
end

% for each block of each fly
for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)   
        meanDataSeq = R(fly).BLOCK(b).meanDataSeq;
        for vol = 1:FLIES(fly).BLOCKS(b).nVol
            thisVolData = permute(squeeze(meanDataSeq(vol,:,:,:)),[2,3,1]);
            R(fly).BLOCK(b).r2Vol(vol) = calculateR2(thisVolData,six_hertz);
            R(fly).BLOCK(b).rVol(vol) = calculateR(thisVolData,six_hertz);
        end    
    end
end
toc;

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS collapsed across volumes (time)
disp('Calculating R and R^2 collapsed across volumes (time)');
tic;
% all blocks for each fly
for fly = 1:length(FLIES)
    thisFlyData = permute(squeeze(mean(R(fly).ALL.meanDataSeq,1)),[2,3,1]);
    R(fly).ALL.r2 = calculateR2(thisFlyData,six_hertz);
    R(fly).ALL.r = calculateR(thisFlyData,six_hertz);
end

% for each block
for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)
        thisBlockData = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,1)),[2,3,1]);
        R(fly).BLOCK(b).r2 = calculateR2(thisBlockData,six_hertz);
        R(fly).BLOCK(b).r = calculateR(thisBlockData,six_hertz);
    end
end
toc;

%% calculate mass t-tests for AAAA vs AAAR and RRRR vs RRRA (per volume)

disp('Calculating t-tests for AAAA vs AAAR and RRRR vs RRRA (per volume)');
tic;
% all blocks for each fly (separation just for organisation purposes)
for fly = 1:length(FLIES)
    for vol = 1:FLIES(fly).nVol
        thisVolData = squeeze(R(fly).ALL.dataSeq(vol,:,:,:,:));
        [R(fly).ALL.AAAAvsAAARVol(vol), R(fly).ALL.RRRRvsRRRAVol(vol)]  = calculateTtests(thisVolData);
    end
end

% for each block of each fly
for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)   
        for vol = 1:FLIES(fly).BLOCKS(b).nVol
            thisVolData = squeeze(R(fly).BLOCK(b).dataSeq(vol,:,:,:,:));
            [R(fly).BLOCK(b).AAAAvsAAARVol(vol), R(fly).BLOCK(b).RRRRvsRRRAVol(vol)]  = calculateTtests(thisVolData);
        end    
    end
end
toc;

%% calculate mass t-tests for AAAA vs AAAR and RRRR vs RRRA (collapsed across volumes/time)
disp('Calculating t-tests for AAAA vs AAAR and RRRR vs RRRA');
tic;
% all blocks for each fly
for fly = 1:length(FLIES)
    thisFlyData = squeeze(mean(R(fly).ALL.dataSeq,1));
    [R(fly).ALL.AAAAvsAAAR, R(fly).ALL.RRRRvsRRRA]  = calculateTtests(thisFlyData);
end

% for each block
for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)
        thisBlockData = squeeze(mean(R(fly).BLOCK(b).dataSeq,1));
        [R(fly).BLOCK(b).AAAAvsAAAR, R(fly).BLOCK(b).RRRRvsRRRA]  = calculateTtests(thisBlockData);
    end
end
toc;

%% plotting

outputDirectory = '../2P Results';

% plot results per block
for fly = 1:length(FLIES)
    subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))]);
    if ~exist(subDirectory,'dir')
       mkdir(subDirectory); 
    end
    plotFly(R(fly), FLIES(fly).brainImage, subDirectory);
end

%% fit and plot some seq eff profiles of interest

% options = optimset('Algorithm','interior-point','FinDiffType','central');
%
% %10 Feb both
% % profiles = [23 19;28 14;23 5;23 19;26 4];
% 
% %10 Feb exp3
% % profiles = [13 18;23 28;21 12; 24 31; 26 3];
% 
% %10 Feb exp4
% profiles = [23 18; 28 17; 23 6; 28 14; 10 14];
% 
% type = [1 2 3 4 4 5];
% 
% names = {'slrp','lrpr','weird','combined','ephys'};
% 
% for i = 1:size(profiles,1)
%     
%     seq_eff_pattern = squeeze(meanVolData(profiles(i,1),profiles(i,2),:));
%    
%     switch type(i)
%         
%         case 1
%   
%            %fit only to slrp or lrpr
%            [x,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 0 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
%            seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%          
%         case 2
%             
%            [x,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 0 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
%            seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%             
%         case 3
%             
%            [x,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
%            seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%             
%         case 4
%             
%            %fit overall
%             [x,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf -inf 0 -inf],[inf inf 0 inf],[],options);
%             seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%             
%         case 5
%             
%             [x,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,six_hertz,seq_eff_pattern),[0 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
%             seq_eff_fit =  x(3)*six_hertz + x(4);
% 
%     end
%     
%     figure; create_seq_eff_plot(seq_eff_pattern,seq_eff_fit);
%     saveas(gcf,[ names{type(i)} '_best.png'])
%     
% end