function R = analyse2PBlock3D(block)

    %% sort data according to previous sequence

%     n_back = 5; 
    
    % yields 6D matrix with (vol,seq,pixelX,pixelY,Z,trial)
    [dataSeq, dataSeqIso] = sortSEs2P3D(block.greenChannel, block.randomSequence, block.nVol, block.nStimuli);
    
    % mean across sixth (trial) dimension (much faster than nan mean)
    meanDataSeq = sum(dataSeq,6)./sum(dataSeq~=0,6);
    
    % this fixes the possibility that all images in a stack have zero
    % values for a particular pixel (usually from alignment)
    meanDataSeq(isnan(meanDataSeq)) = 0;
    
    %number of images for each type of sequence (why (16,16) pixel?)
%     nData = sum(dataSeq(1,:,16,16,1,:)~=0, 6);

    % put all results into a neat structure
    R = struct;

    % calculate mean overall transient (will only be used when a blank baseline is not available)
    % the 3D version of this is done by analogy with 2D but is untested
    R.meanTransient = mean(reshape(block.greenChannel,[size(block.greenChannel,1)...
    size(block.greenChannel,2) size(block.greenChannel,3) block.nVol size(block.greenChannel,4)/block.nVol]),5);

    % calculate mean of blank stack
    if ~isempty(block.blankImageStack)
        meanBlankTransient = mean(reshape(block.blankImageStack,[size(block.blankImageStack,1)...
            size(block.blankImageStack,2) size(block.blankImageStack,3) block.nVol size(block.blankImageStack,4)/block.nVol]),5);
        R.meanBlankTransient = meanBlankTransient;
    else
        R.meanBlankTransient = [];
    end
  
    R.dataSeq = dataSeq;
    R.dataSeqIso = dataSeqIso;% this still holds the 32 sequences
    R.meanDataSeq = meanDataSeq;
%     R.nData = nData;
    
end