close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');

nSlices = 2;

imageSize = [512 512];

%% load data

n_blocks = 1;

% structure with necessary info
D = struct; 

two_p_filenames = {'two_slice_stack_not_aligned.raw'};

stimuli_filenames = {'SE_fly1_exp3_10feb2023_1000_trials.mat'};

two_p_length = [3073];

for b = 1:1%n_blocks

    % load 2P data
    fid = fopen(two_p_filenames{b}, 'r','b');

    data = fread(fid, 512*512*nSlices*two_p_length(b), 'uint16');

    fclose(fid);

    data = permute(reshape(data, [imageSize nSlices two_p_length(b)]),[2 1 3 4]);

    % get rid of the first 6 volumes (there are two triggers before stimuli start)
%     data = data(:,:,14:end);

    % load sequence data
    randomSequence = load(stimuli_filenames{b},'randomSequence');

    % average images in a grid pattern

%     squareSize = 16;
% 
%     gridSize = imageSize/squareSize;
% 
%     rData = zeros(size(data,1)/squareSize,size(data,2)/squareSize,size(data,3));
% 
%     avgFun = @(x) mean(x.data(:));
% 
%     for i = 1:size(data,3)
%         rData(:,:,i) = blockproc(data(:,:,i),[squareSize, squareSize],avgFun);
%     end
%     
%     rData = normalize(rData,3);
%     
%     D(b).randomSequence = randomSequence.randomSequence-1;
%     D(b).rData = rData;
%     D(b).nVol = 3;

end

%% restrict to some blocks

% D = D(2);

%% analyse SEs

R = analyseSequentialEffectsTwoPhoton(D,false);

%% save results

% save('data_two_photon','R');

%% calculate fit to SLRP and LRPR for nVol time steps

load slrp_lrpr.mat

% slrp = rand(16,1);
% lrpr = rand(16,1);

options = optimset('Algorithm','interior-point','FinDiffType','central');

nVol = D.nVol;

meanDataSeq = R.meanDataSeq;

r_squared_slrp = zeros(gridSize);
r_squared_lrpr = zeros(gridSize);
r_squared_overall = zeros(gridSize);
r_squared_weird = zeros(gridSize);

test = zeros(16,16);

% for each volume in time series
for vol = 3:3%nVol
    
    thisVolData = permute(squeeze(meanDataSeq(vol,:,:,:)),[2,3,1]);
    
    for xInd = 1:gridSize

        for yInd = 1:gridSize
            
            test(xInd,yInd) = thisVolData(xInd,yInd,1);

            seq_eff_pattern = squeeze(thisVolData(xInd,yInd,:));

            %fit only to slrp or lrpr
           [x_slrp,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
           [x_lrpr,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
           [x_weird,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);

           %fit overall
           [x_overall,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);

           sse_total = sum((seq_eff_pattern-mean(seq_eff_pattern)).^2);

           r_squared_slrp(xInd,yInd) = 1-(sse_slrp/sse_total);
           r_squared_lrpr(xInd,yInd) = 1-(sse_lrpr/sse_total);
           r_squared_weird(xInd,yInd) = 1-(sse_weird/sse_total);
           r_squared_overall(xInd,yInd) = 1-(sse_overall/sse_total);

        end

    end
    
end

%% plot r_squared

figure; title('SLRP');
imagesc(r_squared_slrp); colorbar; colormap('hot');
xlabel('X'); ylabel('Y');
saveas(gcf,'slrp.png');

figure; title('LRPR');
imagesc(r_squared_lrpr); colorbar; colormap('hot');
xlabel('X'); ylabel('Y');
saveas(gcf,'lrpr.png');

figure; title('WEIRD');
imagesc(r_squared_weird); colorbar; colormap('hot');
xlabel('X'); ylabel('Y');
saveas(gcf,'weird.png');

figure; title('COMBINED ');
imagesc(r_squared_overall); colorbar; colormap('hot');
xlabel('X'); ylabel('Y');
saveas(gcf,'combined.png');