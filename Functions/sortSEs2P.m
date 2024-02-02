function [dataSeq, dataSeqIso] = sortSEs2P(imageStack, randomSequence, nVol, nStimuli)
%this function sorts ERPs according to the past sequence of events and 

    nBack = 5;
    nSeq = 2^nBack;

    sequenceLength = length(randomSequence);

    % better to pre-allocate but 5th dim will be too long
    dataSeq = zeros([nVol nSeq/2 size(imageStack)]);
    dataSeqIso = zeros([nVol nSeq size(imageStack)]);

    % groups (1,32),(2,31),(3,30), etc, as representing the same pattern
    % avoids very costly flip() operations later
%     auxSeq = [1:16 16:-1:1];
    [~,order] = sort(seq_eff_order(5).');
    auxSeq = [order fliplr(order)];
    
    % vol here indexes the number of volumes (time points) collected
    % per trial
    for vol = 1:nVol
        %sort images according to sequence
        for n = 0:(sequenceLength/nStimuli-mod((nStimuli+4),8))% using mod is a hack; check later for vaues different from 1 or 5 (unlikely to be used)
            % decimal value of binary sequence of length n_back
            seq = bin2dec(num2str(randomSequence((n*nStimuli + 1):(n*nStimuli + nBack)))) + 1;

            % stack images for each vol and seq along 5th dimension (separated by pattern)
            dataSeq(vol, auxSeq(seq),:, :, n+1) = imageStack(:,:,n*nVol + vol);
            
            % for the isomers (consumes a lot of memory)
            dataSeqIso(vol,seq,:, :, n+1) = imageStack(:,:,n*nVol + vol);
        end
    end
    
end