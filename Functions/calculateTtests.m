function [AAAAvsAAAR, RRRRvsRRRA] = calculateTtests(dataSeq)
%calculateTtests Summary of this function goes here
%   Detailed explanation goes here

    AAAAvsAAAR = struct;
    RRRRvsRRRA = struct;
    
    gridSize = [size(dataSeq,2) size(dataSeq,3)];
    
    h_AAAAvsAAAR = zeros(gridSize);
    p_AAAAvsAAAR = zeros(gridSize);
    
    h_RRRRvsRRRA = zeros(gridSize);
    p_RRRRvsRRRA = zeros(gridSize);
    
    % because of parfor BS
    xSize = gridSize(1);
    ySize = gridSize(2);
    
    AAAA = squeeze(dataSeq(16,:,:,:));
    AAAR = squeeze(dataSeq(9,:,:,:));
    RRRR = squeeze(dataSeq(1,:,:,:));
    RRRA = squeeze(dataSeq(8,:,:,:));
    
    parfor xInd = 1:xSize
        for yInd = 1:ySize
            
           thisAAAA = squeeze(AAAA(xInd,yInd,:));
           thisAAAR = squeeze(AAAR(xInd,yInd,:));

           [h,p] = ttest2(thisAAAA(logical(thisAAAA)),thisAAAR(logical(thisAAAR)));
           h_AAAAvsAAAR(xInd,yInd) = h; 
           p_AAAAvsAAAR(xInd,yInd) = p;
           
           thisRRRR = squeeze(RRRR(xInd,yInd,:));
           thisRRRA = squeeze(RRRA(xInd,yInd,:));
           
           [h,p] = ttest2(thisRRRR(logical(thisRRRR)),thisRRRA(logical(thisRRRA)));
           h_RRRRvsRRRA(xInd,yInd) = h;
           p_RRRRvsRRRA(xInd,yInd) = p;

        end
    end
    
    AAAAvsAAAR.h = h_AAAAvsAAAR;
    AAAAvsAAAR.p = p_AAAAvsAAAR;
    
    RRRRvsRRRA.h = h_RRRRvsRRRA;
    RRRRvsRRRA.p = p_RRRRvsRRRA;

end

