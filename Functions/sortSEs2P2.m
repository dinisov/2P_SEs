function dataSeq = sortSEs2P2(imageStack,randomSequence, nVol)
%this function sorts ERPs according to the past sequence of events and 

    nBack = 5;
    nSeq = 2^nBack;

    sequenceLength = length(randomSequence);

    % better to pre-allocate but 5th dim will be too long
    dataSeq = zeros([nVol nSeq/2 size(imageStack)]);

    % groups (1,32),(2,31),(3,30), etc, as representing the same pattern
    % avoids very costly flip() operations later
    auxSeq = [1:16 16:-1:1];
    
    % vol here indexes the number of volumes (time points) collected
    % per trial
    for vol = 1:nVol
        %sort images according to sequence
        for n = 0:(sequenceLength/nBack-1)
            % decimal value of binary sequence of length n_back
            seq = auxSeq(bin2dec(num2str(randomSequence(n*nBack + 1:n*nBack + (nBack-1)))) + 1);

            % stack images for each vol and seq along 5th dimension
            dataSeq(vol, seq,:, :, n+1) = imageStack(:,:,n*nVol + vol);
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