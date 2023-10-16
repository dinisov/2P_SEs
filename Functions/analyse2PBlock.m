function R = analyse2PBlock(block)

    %% sort data according to previous sequence

    n_back = 5; 
    
    % yields 5D matrix with (vol,seq,pixelX,pixelY,trial)
    [dataSeq, dataSeqIso] = sortSEs2P(block.greenChannel, block.randomSequence, block.nVol, block.nStimuli);
    
    % mean across fifth (trial) dimension (much faster than nan mean)
    meanDataSeq = sum(dataSeq,5)./sum(dataSeq~=0,5);
    
    % this fixes the possibility that all images in a stack have zero
    % values for a particular pixel (usually from alignment)
    meanDataSeq(isnan(meanDataSeq)) = 0;
    
    %number of images for each type of sequence
    nData = sum(dataSeq(1,:,16,16,:)~=0, 5);
    
    % reorder according to the literature
    meanDataSeq = meanDataSeq(:,seq_eff_order(n_back),:,:);
%     semERPs = semERPs(:,seq_eff_order(n_back));
    
    % we can make a figure showing the average image as a function of the
    % sequence, for each time point
%     figure;
%     plot(meanERPs);

    % put all results into a neat structure
    R = struct;
     
    R.dataSeq = dataSeq;
    R.dataSeqIso = dataSeqIso;% this still holds the 32 sequences
    R.meanDataSeq = meanDataSeq;
    R.nData = nData;
    
end