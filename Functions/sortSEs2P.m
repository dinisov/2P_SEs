function [dataSeq, dataSeqIso, dataSeqBehav] = sortSEs2P(imageStack, randomSequence, nVol, nStimuli, options)
%this function sorts ERPs according to the past sequence of events and 
arguments
    imageStack double
    randomSequence double
    nVol double
    nStimuli double
    options.nBack (1,1) {mustBeNumeric} = 5
    options.behavSequence double = []
end

    %nBack = 5; %Now included as part of arguments
    %nSeq = 2^nBack;
    nSeq = 2^options.nBack;

    sequenceLength = length(randomSequence);

    % better to pre-allocate but 5th dim will be too long
    dataSeq = zeros([nVol nSeq/2 size(imageStack)]);
    dataSeqIso = zeros([nVol nSeq size(imageStack)]);
    if isfield( options, 'behavSequence' ) && ~isempty(options.behavSequence)
        hasBehav = 1;
        %behavSeq = zeros([nVol nSeq/2 size(imageStack)]); %Parallel data array to dataSeq
        %behavSeq = zeros([nVol nSeq/2 1,1, size(imageStack,3)]); %More intelligent singleton where applicable
            %              vol,  seq, Y,X,         t
        behavSeq = zeros([ size(imageStack,3) , 1 ]); %Omega brain single-axis
            %(Easier for using as a index later)
        dataSeqBehav = struct; %Prepare (but not preallocate)
    else
        hasBehav = 0;
        dataSeqBehav = []; %Class type change, empty to boot
    end

    % groups (1,32),(2,31),(3,30), etc, as representing the same pattern
    % avoids very costly flip() operations later
%     auxSeq = [1:16 16:-1:1];
    %[~,order] = sort(seq_eff_order(5).');
    [~,order] = sort(seq_eff_order(options.nBack).'); %Is this meant to be the same as nBack or was it just coincidence?
    auxSeq = [order fliplr(order)];
    
    % vol here indexes the number of volumes (time points) collected
    % per trial
    tic
    for vol = 1:nVol
        %sort images according to sequence
        for n = 0:(sequenceLength/nStimuli-mod((nStimuli+4),8))% using mod is a hack; check later for vaues different from 1 or 5 (unlikely to be used)
            % decimal value of binary sequence of length n_back
            %seq = bin2dec(num2str(randomSequence((n*nStimuli + 1):(n*nStimuli + nBack)))) + 1;
            seq = bin2dec(num2str(randomSequence((n*nStimuli + 1):(n*nStimuli + options.nBack)))) + 1;
            if hasBehav == 1
                bData = median( options.behavSequence((n*nStimuli + 1):(n*nStimuli + options.nBack)) ); %Get the average behavioural state
                %behavSeq(vol, auxSeq(seq),:, :, n+1) = bData; %Save a single value into this data array to represent the behavioural state
                behavSeq( n+1 ) = bData; %Save a single value into this data array to represent the behavioural state
            end
            
            % stack images for each vol and seq along 5th dimension (separated by pattern)
            dataSeq(vol, auxSeq(seq),:, :, n+1) = imageStack(:,:,n*nVol + vol);
            
            % for the isomers (consumes a lot of memory)
            dataSeqIso(vol,seq,:, :, n+1) = imageStack(:,:,n*nVol + vol);
        end
    end
    disp([num2str(toc),'s to assemble dataSeq/dataSeqIso'])
    %size(behavSeq)
    
    %slice out from arrays if behavData present
    if hasBehav == 1
        nonEmptyInds = find( squeeze( sum(dataSeq, [1:4]) ) ~= 0 ); %Do only once to save time
        tic
        unBehavs = unique( behavSeq );
        %QA
        if length(unBehavs) > 4
            ['-# Alert: Potentially anomalous number of unique behavioural states #-']
            crash = yes
            %Considering how long it takes to process a state, this is a safeguard
        end
        for unI = 1:length(unBehavs)
            disp(['Searching for behaviour state ',num2str(unBehavs(unI))])
            thisInds = find( behavSeq == unBehavs(unI) );
            disp(['Found ',num2str(length(thisInds)),' instances'])
            %dataSeqBehav(unI).dataSeq = dataSeq( :,:,:,:, thisInds  );
            thisInds = intersect( thisInds, nonEmptyInds ); %Reduce, since dataSeq is len 5000 unnecessarily
            disp(['(Reduced to ',num2str(length(thisInds)),' non-zero instances)'])
            dataSeqReduced = dataSeq( :,:,:,:, thisInds );
            dataSeqReduced(isnan(dataSeqReduced)) = 0; %Just in case
            meanDataSeqReduced = mean( dataSeqReduced, 5 ); %Slightly different to how done in analyse2PBlock, but conceptually same?
            
            dataSeqBehav(unI).state = unBehavs(unI); %What state this referred to
            dataSeqBehav(unI).dataSeqReduced = dataSeqReduced;
            dataSeqBehav(unI).meanDataSeqReduced = meanDataSeqReduced;
            dataSeqBehav(unI).behavSeq = behavSeq; %Save information (with every state, but w/e)
            %disp([num2str(size(dataSeq)),' -> ',num2str(size(dataSeqBehav(unI).dataSeq))])
            disp([num2str(size(dataSeq)),' -> ',num2str(size(dataSeqBehav(unI).dataSeqReduced))])
        end
        disp([num2str(toc),'s to assemble dataSeqBehav for ',num2str(unI),' unique state/s'])
    end
    
end