function LvsR = calculateTtestsLR(dataSeqIso)
%calculateTtests Summary of this function goes here
%   Detailed explanation goes here

    LvsR = struct;
    
    gridSize = [size(dataSeqIso,2) size(dataSeqIso,3)];
    
    h_LvsR = zeros(gridSize);
    p_LvsR = zeros(gridSize);
    
    % because of parfor BS
    xSize = gridSize(1);
    ySize = gridSize(2);
    
    L = squeeze(sum(dataSeqIso(1:2:31,:,:,:),1));
    R = squeeze(sum(dataSeqIso(2:2:32,:,:,:),1));
    
    parfor xInd = 1:xSize
        for yInd = 1:ySize
            
           thisL = squeeze(L(xInd,yInd,:));
           thisR = squeeze(R(xInd,yInd,:));

           [h,p] = ttest2(thisL(logical(thisL)),thisR(logical(thisR)));
           h_LvsR(xInd,yInd) = h; 
           p_LvsR(xInd,yInd) = p;

        end
    end
    
    meanR = sum(R,3)./sum(R~=0,3);
    meanL = sum(L,3)./sum(L~=0,3);
    
    meanR(isnan(meanR)) = 0;
    meanL(isnan(meanL)) = 0;
    
    LvsR.h = h_LvsR;
    LvsR.p = p_LvsR;
    LvsR.LRDiff = meanL-meanR;
    LvsR.meanL = meanL;
    LvsR.meanR = meanR;
%     LvsR.dataL = L;
%     LvsR.dataR = R;