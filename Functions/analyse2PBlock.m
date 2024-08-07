function R = analyse2PBlock(block)

    %% sort data according to previous sequence

%     n_back = 5; 

    %disp( ['block.randomSequence - ', class(block.randomSequence)] )
    %disp( ['block.greenChannel - ', class(block.greenChannel)] )
    %load('I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\3Jun24\fly1_exp1_3Jun24\behavSequence.mat')
    %behavSequence = savStruct.acInac.thisInacBinaryInterp;

    % yields 5D matrix with (vol,seq,pixelX,pixelY,trial)
    %[dataSeq, dataSeqIso] = sortSEs2P(block.greenChannel, block.randomSequence, block.nVol, block.nStimuli);
    [dataSeq, dataSeqIso, dataSeqBehav] = sortSEs2P(block.greenChannel, block.randomSequence, block.nVol, block.nStimuli, 'behavSequence', block.behavSequence);
        %We request dataSeqBehav even if we won't use it (It's empty in such cases)
    
    % mean across fifth (trial) dimension (much faster than nan mean)
    meanDataSeq = sum(dataSeq,5)./sum(dataSeq~=0,5);
    
    % this fixes the possibility that all images in a stack have zero
    % values for a particular pixel (usually from alignment)
    meanDataSeq(isnan(meanDataSeq)) = 0;
    
    %number of images for each type of sequence
    nData = sum(dataSeq(1,:,16,16,:)~=0, 5);

    % put all results into a neat structure
    R = struct;

    % calculate mean overall transient
    R.meanTransient = mean(reshape(block.greenChannel,[size(block.greenChannel,1)...
    size(block.greenChannel,2) block.nVol size(block.greenChannel,3)/block.nVol]),4);

    % calculate mean of blank stack
    if ~isempty(block.blankImageStack)
        meanBlankTransient = mean(reshape(block.blankImageStack,[size(block.blankImageStack,1)...
            size(block.blankImageStack,2) block.nVol size(block.blankImageStack,3)/block.nVol]),4);
        R.meanBlankTransient = meanBlankTransient;
    else
        R.meanBlankTransient = [];
    end
  
    R.dataSeq = dataSeq;
    R.dataSeqIso = dataSeqIso;% this still holds the 32 sequences
    R.meanDataSeq = meanDataSeq;
    R.nData = nData;
    
    % these avoid later passing FLIES without need
    R.blockNum = block.blockNum;
    R.Trim = block.Trim;
    R.nVol = block.nVol;
    
    %Tack on behav-separated data if applicable
    if ~isempty(dataSeqBehav)
        R.dataSeqBehav = dataSeqBehav; %Note: Structure, not double array    
    end
    
end