function R = ttests2P(R, groupedBlocks)
%ttests2P Summary of this function goes here
%   Detailed explanation goes here

%% calculate mass t-tests for AAAA vs AAAR and RRRR vs RRRA (per volume)
disp('Calculating t-tests for AAAA vs AAAR and RRRR vs RRRA (per volume)');
tic;

if groupedBlocks
    % all blocks for each fly (separation just for organisation purposes)
    for fly = 1:length(R)
        for vol = 1:R(fly).nVol
            thisVolData = squeeze(R(fly).ALL.dataSeq(vol,:,:,:,:));
            [R(fly).ALL.AAAAvsAAARVol(vol), R(fly).ALL.RRRRvsRRRAVol(vol)]  = calculateTtests(thisVolData);
        end
    end
end

% for each block of each fly
for fly = 1:length(R)
    for b = [R(fly).BLOCK.blockNum]  
        for vol = 1:R(fly).BLOCK(b).nVol
            thisVolData = squeeze(R(fly).BLOCK(b).dataSeq(vol,:,:,:,:));
            [R(fly).BLOCK(b).AAAAvsAAARVol(vol), R(fly).BLOCK(b).RRRRvsRRRAVol(vol)]  = calculateTtests(thisVolData);
        end    
    end
end
toc;

%% calculate mass t-tests for AAAA vs AAAR and RRRR vs RRRA (collapsed across volumes/time)
disp('Calculating t-tests for AAAA vs AAAR and RRRR vs RRRA');
tic;

if groupedBlocks
    % all blocks for each fly
    for fly = 1:length(R)
        thisFlyData = squeeze(mean(R(fly).ALL.dataSeq,1));
        [R(fly).ALL.AAAAvsAAAR, R(fly).ALL.RRRRvsRRRA]  = calculateTtests(thisFlyData);
    end
end

% for each block
for fly = 1:length(R)
    for b = [R(fly).BLOCK.blockNum]  
        thisBlockData = squeeze(mean(R(fly).BLOCK(b).dataSeq,1));
        [R(fly).BLOCK(b).AAAAvsAAAR, R(fly).BLOCK(b).RRRRvsRRRA]  = calculateTtests(thisBlockData);
    end
end
toc;
end

