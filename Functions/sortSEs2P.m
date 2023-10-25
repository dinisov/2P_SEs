function [dataSeq, dataSeqIso, blankStack] = sortSEs2P(imageStack, randomSequence, nVol, nStimuli)
%this function sorts ERPs according to the past sequence of events and 

    nBack = 5;
    nSeq = 2^nBack;
    
    % better to pre-allocate but 5th dim will be too long
    dataSeq = zeros([nVol nSeq/2 size(imageStack)]);
    dataSeqIso = zeros([nVol nSeq size(imageStack)]);

    % groups (1,32),(2,31),(3,30), etc, as representing the same pattern
    % avoids very costly flip() operations later
%     auxSeq = [1:16 16:-1:1];
    [~,order] = sort(seq_eff_order(5));
    auxSeq = [order fliplr(order)];
    
    % index of blank stimuli/blocks
    indBlank = randomSequence == 5;
    
    % remove blank stimuli from random sequence
    randomSequence = randomSequence(~indBlank);
    sequenceLength = length(randomSequence);
    
    % construct index of frames after blank blocks
    indBlankStack = reshape(indBlank,[5 length(indBlank)/5]).';
    indBlankStack = repmat(indBlankStack(:,1),[1 nVol]).';
    
    % split between stimuli image and blank image stacks
    blankStack = imageStack(:,:,indBlankStack(:));
    imageStack = imageStack(:,:,~indBlankStack(:));
    
    % vol here indexes the number of volumes (time points) collected
    % per trial
    for vol = 1:nVol
        %sort images according to sequence
        for n = 0:(sequenceLength/nStimuli-mod((nStimuli+3),7))% using mod is a hack; check later for vaues different from 1 or 5 (unlikely to be used)
            % decimal value of binary sequence of length n_back
            seq = bin2dec(num2str(randomSequence((n*nStimuli + 1):(n*nStimuli + (nBack-1))))) + 1;

            % stack images for each vol and seq along 5th dimension (separated by pattern)
            dataSeq(vol, auxSeq(seq),:, :, n+1) = imageStack(:,:,n*nVol + vol);
            
            % for the isomers (consumes a lot of memory)
            dataSeqIso(vol,seq,:, :, n+1) = imageStack(:,:,n*nVol + vol);
        end
    end

    %%

    %calculate mean and SEM for outlier calculations
%         meanERP = mean(all_erps, 2);
%         STDs = std(all_erps, [], 2);
% 
%         n_sd = 3;

    % remove ERPs beyond n_sd (broadcasting here)
%         outliers = all_erps < (meanERP - n_sd*STDs) | all_erps > (meanERP + n_sd*STDs);
% 
%         disp(['Data lost due to outliers: ' num2str(nnz(sum(outliers))/length(outliers)*100) '%']);

%         good_erps = ~logical(sum(outliers));

    % remove ERP outliers
%         ERPS = ERPS(:,:,good_erps);
%         seqPHOT = seqPHOT(:,:,good_erps);
%         badTrials = badTrials(good_erps);

%         if isfield(blocks,'focusPeaks')
%             focusPeaks = focusPeaks(good_erps);
%             blocks(b).focusPeaks = focusPeaks;
%         end

%     blocks(b).rDataSeq = rDataSeq;
    
end