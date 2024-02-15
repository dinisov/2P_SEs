function R = componentFits2P3D(R, FLIES, groupedBlocks)
%componentFits2P Summary of this function goes here
%   Detailed explanation goes here
load slrp_lrpr.mat

load six_hertz.mat

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS per volume

disp('Calculating R and R^2 per volume');
tic;

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

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS collapsed across volumes (time)
disp('Calculating R and R^2 collapsed across volumes (time)');

% for each block
for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)
        thisBlockData = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,1)),[2,3,1]);
        R(fly).BLOCK(b).r2 = calculateR2(thisBlockData,six_hertz);
        R(fly).BLOCK(b).r = calculateR(thisBlockData,six_hertz);
    end
end

if groupedBlocks
    % all blocks for each fly
    for fly = 1:length(FLIES)
        thisFlyData = permute(squeeze(mean(R(fly).ALL.meanDataSeq,1)),[2,3,1]);
        R(fly).ALL.r2 = calculateR2(thisFlyData,six_hertz);
        R(fly).ALL.r = calculateR(thisFlyData,six_hertz);
    end
end
toc;

end