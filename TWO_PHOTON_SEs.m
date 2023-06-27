close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_swinderen\Dinis\2P SEs\Functions\');

flyID = 'fly1_exp1_15Jun23';

% where the 2P data is located
twoPDirectory = ['\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\15Jun23\' flyID];

%where the sequence data is located (stimulus setup files)
sequenceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\SE_2P_data';

imageSize = [512 512];

squareSize = 4;

gridSize = imageSize/squareSize;

I = imread('cc_2.jpg');

%% load 2P data

n_blocks = 1;

% structure with necessary info
D = struct;

% two_p_filenames = {'AVG_aligned_exp3_24Mar23.raw','AVG_aligned_exp4_24Mar23.raw'};
% 
% stimuli_filenames = {'SE_fly1_exp3_24Mar2023_1000_trials.mat','SE_fly1_exp4_24Mar2023_1000_trials.mat'};

% two_p_filenames = {'AVG_Z_aligned_fly2_exp2_23Feb'};

% green channel file

% red channel (if it exists)

% either (1) the full length of recording (number of slices irrespective of volume or trial) 
% or (2) a desired number of slices to load (hopefully a multiple of nVol)
two_p_length = [55000];

nVol = [11];% number of volumes in each recording epoch

%% load and reduce data
tic
for b = 1:1%n_blocks

    % load 2P data
    fid = fopen([twoPDirectory '\green_channel.raw'], 'r','b');

    data = fread(fid, 512*512*two_p_length(b), 'uint16');

    fclose(fid);

    data = permute(reshape(data, [imageSize two_p_length(b)]),[2 1 3]);

    % load sequence data
%     randomSequence = load(stimuli_filenames{b},'randomSequence');
    randomSequence = csvread(fullfile(sequenceDirectory,'Data_LEDs',[flyID '.csv'])).';

    %% average images in a grid pattern

    rData = zeros(size(data,1)/squareSize,size(data,2)/squareSize,size(data,3));

%     avgFun = @(x) mean(x.data(:));
% 
%     for i = 1:size(data,3)
%         rData(:,:,i) = blockproc(data(:,:,i),[squareSize, squareSize],avgFun);
% %         rData(:,:,i) = imresize(data(:,:,i),[size(data,1)/squareSize, size(data,2)/squareSize],'box');
%     end

    % man oh man this is faster
    rData = imresize3(data,size(rData),'box');
    
    %rData = normalize(rData,3);
%     rData = rData(:,:,randperm(3060));
    
    D(b).randomSequence = randomSequence;
    D(b).rData = rData;
    D(b).nVol = nVol(b);

end
toc
%% restrict to some blocks

% D = D(2);

%% analyse SEs

R = analyseSequentialEffectsTwoPhoton(D,false);

%% save results

% save('data_two_photon_24Mar23_exp4','R');

%% load results as a shortcut

% load('data_two_photon_24Mar23_exp4','R');

%% calculate fit to SLRP and LRPR for nVol time steps

load slrp_lrpr.mat

load six_hertz.mat

% slrp = rand(16,1);
% lrpr = rand(16,1);
% weird = rand(16,1);
% six_hertz = rand(16,1);

meanDataSeq = R.meanDataSeq;

%%  r and r^2 for each volume in time series

for vol = 1:nVol(1)
    
    thisVolData = permute(squeeze(meanDataSeq(vol,:,:,:)),[2,3,1]);
    
    R2Vol(vol) = calculateR2(thisVolData,six_hertz); %#ok<SAGROW>
    
    RVol(vol) = calculateR(thisVolData,six_hertz); %#ok<SAGROW>
    
end

%% r and r^2 average across all volumes (time steps)

meanVolData = permute(squeeze(sum(meanDataSeq,1)),[2,3,1]);

R2All = calculateR2(meanVolData,six_hertz);

RAll = calculateR(meanVolData,six_hertz);

%% plot AAAA-AAAR and RRRR-RRRA

outputDirectory = 'Results/';

plotAR(meanDataSeq,10,I,outputDirectory);

%% plot r_squared for the average of all volumes

plotR2(R2All, outputDirectory);

%% plot r for the average across all volumes

plotR(RAll, outputDirectory);

%% r plots thresholded by significance level (not overlayed on fly's brain)

% figure; h = imagesc(r_slrp); colorbar; set(h,'alphadata',p_slrp); set(gca,'color','k'); caxis([-1 1]);
% figure; h = imagesc(r_lrpr); colorbar; set(h,'alphadata',p_lrpr); set(gca,'color','k'); caxis([-1 1]);
% figure; h = imagesc(r_weird); colorbar; set(h,'alphadata',p_weird); set(gca,'color','k'); caxis([-1 1]);
% figure; h = imagesc(r_ephys); colorbar; set(h,'alphadata',p_ephys); set(gca,'color','k'); caxis([-1 1]);

%% plot r overlayed in fly's brain (thresholded by significance)

signLevel = 0.05;

plotBrain(RAll, I, signLevel, outputDirectory);

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