close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_swinderen\Dinis\2P SEs\Functions\');

close all; clear;

% mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';
mainDirectory = '../2P Data';

outputDirectory = '../2P Results';

RDMDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\';

%where the sequence data is located (stimulus files)
% sequenceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\SE_2P_data\Data_LEDs';
sequenceDirectory = fullfile(RDMDirectory,'RPiData/');

% scratchDirectory = '../../2P Data';

blocks = readtable('../2P Record/2P_record');

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [32 32];

flyList = unique(blocks.Fly);

chosenFlies = 43:43;
% chosenFlies = flyList;%do not choose any flies

%whether to analyse grouped blocks
groupedBlocks = 0;

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

%% analyse SEs
% separates images according to preceding sequence of stimuli and
% calculates mean images as a function of the sequence

%structure to put results in
R = struct;

disp('Analysing SEs...');
tic;
for fly = 1:length(FLIES)
    disp(['Fly ' num2str(fly)]);
    thisFly = FLIES(fly);
    for b = 1:length(thisFly.BLOCKS)
        disp(['Block ' num2str(b)]);
        thisBlock = thisFly.BLOCKS(b);
        R(fly).BLOCK(b) = analyseSequentialEffectsTwoPhoton2(thisBlock.greenChannel,thisBlock.randomSequence,thisBlock.nVol);
        
        thisBlockDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]);
        if ~exist(thisBlockDirectory,'dir')
            mkdir(thisBlockDirectory); 
        end
        meanDataSeq = R(fly).BLOCK(b).meanDataSeq;
        save(fullfile(thisBlockDirectory,'results'),'meanDataSeq');
    end
    
    % add brain images to structure
    for b = 1:length(thisFly.BLOCKS)
        R(fly).BLOCK(b).brainImage = FLIES(fly).BLOCKS(b).brainImage;
    end
    
    if groupedBlocks
        disp(['Fly ' num2str(fly) ' grouped blocks']);
        R(fly).ALL = analyseSequentialEffectsTwoPhoton2(thisFly.greenChannel,thisFly.randomSequence,thisFly.nVol);
        R(fly).ALL.brainImage = FLIES(1).BLOCKS(b).brainImage;
        save(fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],'All','results'),'results');%UNTESTED
    end
end
toc;

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS per volume
load slrp_lrpr.mat

load six_hertz.mat

disp('Calculating R and R^2 per volume');
tic;

if groupedBlocks
    % grouped blocks for each fly (separation just for organisation purposes)
    for fly = 1:length(FLIES)
        meanDataSeq = R(fly).ALL.meanDataSeq;
        for vol = 1:FLIES(fly).nVol
            thisVolData = permute(squeeze(meanDataSeq(vol,:,:,:)),[2,3,1]);
            R(fly).ALL.r2Vol(vol) = calculateR2(thisVolData,six_hertz);
            R(fly).ALL.rVol(vol) = calculateR(thisVolData,six_hertz);
        end
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

%% make movies of fits over time

trim = 3;

newGridSize = gridSize-2*trim;

for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)
        subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'Movies');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        
        %alternation component movie
        movieMatrix = zeros([newGridSize FLIES(fly).BLOCKS(b).nVol]);
        for vol = 1:FLIES(fly).BLOCKS(b).nVol
            movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_slrp(trim+1:end-trim,trim+1:end-trim);
        end
        makeMovie(movieMatrix,fullfile(subDirectory,'ALT_fit.avi'),true);
        
        %repetition component movie
        movieMatrix = zeros([newGridSize FLIES(fly).BLOCKS(b).nVol]);
        for vol = 1:FLIES(fly).BLOCKS(b).nVol
            movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_lrpr(trim+1:end-trim,trim+1:end-trim);
        end
        makeMovie(movieMatrix,fullfile(subDirectory,'REP_fit.avi'),true);
        
        %weird component movie
        movieMatrix = zeros([newGridSize FLIES(fly).BLOCKS(b).nVol]);
        for vol = 1:FLIES(fly).BLOCKS(b).nVol
            movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_weird(trim+1:end-trim,trim+1:end-trim);
        end
        makeMovie(movieMatrix,fullfile(subDirectory,'WEIRD_fit.avi'),true);
    end
end

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS collapsed across volumes (time)
% disp('Calculating R and R^2 collapsed across volumes (time)');
% tic;
% 
% if groupedBlocks
%     % all blocks for each fly
%     for fly = 1:length(FLIES)
%         thisFlyData = permute(squeeze(mean(R(fly).ALL.meanDataSeq,1)),[2,3,1]);
%         R(fly).ALL.r2 = calculateR2(thisFlyData,six_hertz);
%         R(fly).ALL.r = calculateR(thisFlyData,six_hertz);
%     end
% end
% 
% % for each block
% for fly = 1:length(FLIES)
%     for b = 1:length(FLIES(fly).BLOCKS)
%         thisBlockData = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,1)),[2,3,1]);
%         R(fly).BLOCK(b).r2 = calculateR2(thisBlockData,six_hertz);
%         R(fly).BLOCK(b).r = calculateR(thisBlockData,six_hertz);
%     end
% end
% toc;

%% calculate mass t-tests for AAAA vs AAAR and RRRR vs RRRA (per volume)

% disp('Calculating t-tests for AAAA vs AAAR and RRRR vs RRRA (per volume)');
% tic;
% 
% if groupedBlocks
%     % all blocks for each fly (separation just for organisation purposes)
%     for fly = 1:length(FLIES)
%         for vol = 1:FLIES(fly).nVol
%             thisVolData = squeeze(R(fly).ALL.dataSeq(vol,:,:,:,:));
%             [R(fly).ALL.AAAAvsAAARVol(vol), R(fly).ALL.RRRRvsRRRAVol(vol)]  = calculateTtests(thisVolData);
%         end
%     end
% end
% 
% % for each block of each fly
% for fly = 1:length(FLIES)
%     for b = 1:length(FLIES(fly).BLOCKS)   
%         for vol = 1:FLIES(fly).BLOCKS(b).nVol
%             thisVolData = squeeze(R(fly).BLOCK(b).dataSeq(vol,:,:,:,:));
%             [R(fly).BLOCK(b).AAAAvsAAARVol(vol), R(fly).BLOCK(b).RRRRvsRRRAVol(vol)]  = calculateTtests(thisVolData);
%         end    
%     end
% end
% toc;

%% calculate mass t-tests for AAAA vs AAAR and RRRR vs RRRA (collapsed across volumes/time)
% disp('Calculating t-tests for AAAA vs AAAR and RRRR vs RRRA');
% tic;
% 
% if groupedBlocks
%     % all blocks for each fly
%     for fly = 1:length(FLIES)
%         thisFlyData = squeeze(mean(R(fly).ALL.dataSeq,1));
%         [R(fly).ALL.AAAAvsAAAR, R(fly).ALL.RRRRvsRRRA]  = calculateTtests(thisFlyData);
%     end
% end
% 
% % for each block
% for fly = 1:length(FLIES)
%     for b = 1:length(FLIES(fly).BLOCKS)
%         thisBlockData = squeeze(mean(R(fly).BLOCK(b).dataSeq,1));
%         [R(fly).BLOCK(b).AAAAvsAAAR, R(fly).BLOCK(b).RRRRvsRRRA]  = calculateTtests(thisBlockData);
%     end
% end
% toc;

%% RRRR-RRRA and AAAA-AAAR (collapsed and over time videos)

trim = 3;

% newGridSize = gridSize-2*trim;

for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)
        subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'Oddballs');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
       thisBlockData = squeeze(mean(R(fly).BLOCK(b).meanDataSeq,1));
       AAAA = squeeze(thisBlockData(16,trim+1:end-trim,trim+1:end-trim)); AAAR = squeeze(thisBlockData(8,trim+1:end-trim,trim+1:end-trim)); 
       RRRR = squeeze(thisBlockData(1,trim+1:end-trim,trim+1:end-trim)); RRRA = squeeze(thisBlockData(9,trim+1:end-trim,trim+1:end-trim));
       
       figure; imagesc(AAAA); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAA.png')); 
       figure; imagesc(RRRR); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRR.png'));
       figure; imagesc(RRRA); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRA.png'));
       figure; imagesc(AAAR); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAR.png'));
       figure; imagesc((AAAA-AAAR)); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAAminusAAAR.png'));
       figure; imagesc((RRRR-RRRA)); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRRminusRRRA.png'));
       close all;
       
       %videos of differences over time
       
       %AAAA-AAAR
       AAAAminusAAAR_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,16,trim+1:end-trim,trim+1:end-trim)),[2 3 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,8,trim+1:end-trim,trim+1:end-trim)),[2 3 1]);
       RRRRminusRRRA_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,1,trim+1:end-trim,trim+1:end-trim)),[2 3 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,9,trim+1:end-trim,trim+1:end-trim)),[2 3 1]);
       
       AAAAminusAAAR_t = normalize(AAAAminusAAAR_t,3,'range',[0 1]);
       makeMovie(AAAAminusAAAR_t,fullfile(subDirectory,'AAAAminusAAAR_transient.avi'),false);
       
       RRRRminusRRRA_t = normalize(RRRRminusRRRA_t,3,'range',[0 1]);
       makeMovie(RRRRminusRRRA_t,fullfile(subDirectory,'RRRRminusRRRA_transient.avi'),false);
       
    end
end

%% L vs R analysis (t-tests, L-R, L and R, L and R movies)

trim = 3;

newGridSize = gridSize-2*trim;

disp('Calculating L vs R');
tic;
% for each block
for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)
        subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'LvsR');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        thisBlockData = squeeze(mean(R(fly).BLOCK(b).dataSeqIso,1));
        R(fly).BLOCK(b).LvsR  = calculateTtestsLR(thisBlockData);
        figure; imagesc(R(fly).BLOCK(b).LvsR.h); saveas(gcf,fullfile(subDirectory,'h.png'));
        figure; imagesc(R(fly).BLOCK(b).LvsR.p); colormap(hot(256));saveas(gcf,fullfile(subDirectory,'p.png'));
        
        LRDiff = R(fly).BLOCK(b).LvsR.LRDiff(trim+1:end-trim,trim+1:end-trim); %LRDiff = normalize(LRDiff,'range',[0 1]);
        figure; imagesc(LRDiff); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'LminusR.png'));
        
        meanL = R(fly).BLOCK(b).LvsR.meanL(trim+1:end-trim,trim+1:end-trim);% meanL = normalize(meanL,'range',[0 1]);
        figure; imagesc(meanL); colormap(jet(256));colorbar; saveas(gcf,fullfile(subDirectory,'meanL.png'));
        
        meanR = R(fly).BLOCK(b).LvsR.meanR(trim+1:end-trim,trim+1:end-trim);% meanR = normalize(meanR,'range',[0 1]);
        figure; imagesc(meanR); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'meanR.png'));
        close all;
        
        %make movies for L and R transients and difference
        left = permute(squeeze(mean(mean(R(fly).BLOCK(b).dataSeqIso(:,1:2:31,trim+1:end-trim,trim+1:end-trim,:),5),2)),[2 3 1]);
        left = normalize(left,3,'range',[0 1]);
        makeMovie(left,fullfile(subDirectory,'left_transient.avi'),false);
        
        right = permute(squeeze(mean(mean(R(fly).BLOCK(b).dataSeqIso(:,2:2:32,trim+1:end-trim,trim+1:end-trim,:),5),2)),[2 3 1]);
        right = normalize(right,3,'range',[0 1]);
        makeMovie(right,fullfile(subDirectory,'right_transient.avi'),false);
        
        diffLR = left-right;
        diffLR = normalize(diffLR,3,'range',[0 1]);
        makeMovie(right,fullfile(subDirectory,'difference_transient.avi'),false);
         
    end
end
toc;

%% plotting

% plot results per block
% for fly = 1:length(FLIES)
%     subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))]);
%     if ~exist(subDirectory,'dir')
%        mkdir(subDirectory); 
%     end
%     plotFly(R(fly), groupedBlocks, subDirectory);
% end

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

function makeMovie(data, filename, isCorr)

    VidObj = VideoWriter(filename, 'Uncompressed AVI');

    VidObj.FrameRate = 2; %set your frame rate
    open(VidObj);
    for f = 1:size(data, 3)
        if isCorr
            %if we're dealing with corr coeffs then map [-1 1] to [0 1]
            writeVideo(VidObj, ind2rgb(uint8(255 * (1+data(:,:,f))/2), jet(256)) );
        else
            %otherwise make video of data pre-mapped to [0 1]
            writeVideo(VidObj, ind2rgb(uint8(255 * data(:,:,f)), jet(256)));
        end
    end
    close(VidObj);
       
end