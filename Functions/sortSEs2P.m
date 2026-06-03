function [dataSeq, dataSeqIso, dataSeqBehav] = sortSEs2P(imageStack, randomSequence, nVol, nStimuli, options)
%this function sorts ERPs according to the past sequence of events
%Mk 1 - Dinis original functionality (Note: dataSeq hardcoded to be constructed as block design, even from rolling data)
%   2 - Modifications to use rollStruct code for rolling
%   3 - Reversion to use faux-block nature for rolling

arguments
    imageStack double
    randomSequence double
    nVol double %Note: Currently, for rolling this = nVol per nBack (e.g. 5) train of stimuli (e.g. 11/stimuli * 5 stimuli = nVol of 55 etc), whereas for block this = inter-stimulus period (e.g. 24)
    nStimuli double
    options.nBack (1,1) {mustBeNumeric} = 5
    options.behavSequence double = []
    options.doRolling (1,1) {mustBeNumeric} = 0 %Whether to additionally calculate a rolling version of dataSeq; Disabled currently
    options.skipIso (1,1) {mustBeNumeric} = 0 %Whether to not calculate isomer array (Return empty)
end

%for non-function running
%{
imageStack = block.greenChannel;
randomSequence = block.randomSequence;
nVol = block.nVol;
nStimuli = block.nStimuli;
options.doRolling = block.doRolling;
options.skipIso = 1; %Forced default
options.nBack = 5; %Forced default
%}

    %Reminder: This is the 2D version, so imageStack only has 3 effective dimensions (XYT), because Z is flattened

    %nBack = 5; %Now included as part of arguments
    %nSeq = 2^nBack;
    nSeq = 2^options.nBack;

    skipIso = options.skipIso;

    doRolling = options.doRolling;

    sequenceLength = length(randomSequence);

    %QA in case of bad nStimuli etc
    if (sequenceLength/nStimuli-mod((nStimuli+4),8)) <= 0
        ['## Error: Likely anomalous sequence sorting second index calculated (',num2str((sequenceLength/nStimuli-mod((nStimuli+4),8))),...
            ') ##']
    end

    % better to pre-allocate but 5th dim will be too long
    %dataSeq = zeros([nVol nSeq/2 size(imageStack)]);
    %if ~doRolling
    maxN = (sequenceLength/nStimuli-mod((nStimuli+4),8)) + 1; %Theoretical maximum position a volume of data will ever be assigned to (5th axis)
    dataSeq = zeros([nVol nSeq/2 size(imageStack, [1,2]) maxN ]); %Use maximum n (See below) rather than image frames as last (5th) axis length
        %Note: May cause issues if anything later assumes dataSeq 5th axis size == nFrames etc
    %else %Instead of being pre-allocated here, dataSeq for rolling will be calculated below, once framesPer is known
    %    dataSeq = zeros([nVol nSeq/2 size(imageStack, [1,2]) maxN ]);
    if ~skipIso
        %dataSeqIso = zeros([nVol nSeq size(imageStack)]);
        dataSeqIso = zeros([nVol nSeq size(imageStack, [1,2]) maxN]); %Even more crucial for this to be as small as possible...
    else
        dataSeqIso = [];
    end
    %end
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
    
    %Rolling-specific indexing
    if isfield( options, 'doRolling' ) && options.doRolling == 1 %Note: Currently deprecated?
        starterN = options.nBack - 1; %Rolling, due to faux blocklike/expanded nature starts at nBack-1 (e.g. randomSeq pos 21 -> ...)
    else
        starterN = 0; %Block design, per original Dinis construction starts at 0 (+1 etc)
    end

    %Whether to tack on rolling analysis
        %Might make an actual switch later (i.e. Block OR Rolling)
    %{
    %(Disabled)
    if isfield( options, 'doRolling' ) && options.doRolling == 1
        doRolling = 1;
        rollStruct = struct;
        %%rollStruct.rollSeq = nan([nSeq/2 size(imageStack)]); %Excessive last dim probably warranted here, at least for assembly (Unlike for dataSeq)
            %SXYT
            %Disabled, for memory reasons (Particularly with fast/long recordings)
        rollStruct.rollSeqReduced = nan([nSeq/2 size(imageStack)]);
    else
        doRolling = 0;
        rollStruct = [];
        %rollSeq = [];
    end
    %}

    % groups (1,32),(2,31),(3,30), etc, as representing the same pattern
    % avoids very costly flip() operations later
%     auxSeq = [1:16 16:-1:1];
    %[~,order] = sort(seq_eff_order(5).');
    [~,order] = sort(seq_eff_order(options.nBack).'); %Is this meant to be the same as nBack or was it just coincidence?
    auxSeq = [order fliplr(order)];
    
    % vol here indexes the number of volumes (time points) collected
    % per trial
    %if ~doRolling
        %disp(['~~ Sorting for non-rolling nature ~~'])
        disp(['~~ Sorting ~~'])
        tic
        for vol = 1:nVol
            %sort images according to sequence
            %for n = 0:(sequenceLength/nStimuli-mod((nStimuli+4),8))% using mod is a hack; check later for vaues different from 1 or 5 (unlikely to be used)
            for n = starterN:(sequenceLength/nStimuli-mod((nStimuli+4),8)) %Slight modification to allow for correct start indexing 
                    %For block design the second value here will be the number of blocks(?), whilst for rolling it will be the sequenceLength / 4 (More or less)
                        %Note: Even with rolling this would have still been calculated in a block-like nature (No 1-by-1 iteration/overlap) [I think]
                            %i.e. Even if the sequence was continuous, at length 1385, it would have yielded 277 'blocks'
                % decimal value of binary sequence of length n_back
                %seq = bin2dec(num2str(randomSequence((n*nStimuli + 1):(n*nStimuli + nBack)))) + 1;
                seq = bin2dec(num2str(randomSequence((n*nStimuli + 1):(n*nStimuli + options.nBack)))) + 1;
                if hasBehav == 1
                    bData = median( options.behavSequence((n*nStimuli + 1):(n*nStimuli + options.nBack)) ); %Get the average behavioural state
                    %behavSeq(vol, auxSeq(seq),:, :, n+1) = bData; %Save a single value into this data array to represent the behavioural state
                    behavSeq( n+1 ) = bData; %Save a single value into this data array to represent the behavioural state
                end
                
                % stack images for each vol and seq along 5th dimension (separated by pattern)
                %Quick new QA to make sure new elements aren't being added
                if n+1 > size( dataSeq, 5 )
                    ['## Alert: dataSeq 5th axis insertion position (',num2str(n+1),') larger than pre-allocated size (',num2str(size( dataSeq, 5 )),') ##']
                    crash = yes %Not crucial to crash here, but probably means new/better code needs to be written
                end
                dataSeq(vol, auxSeq(seq),:, :, n+1) = imageStack(:,:,n*nVol + vol); %Unmodified
                    %Note: With rolling, the first nBack-1 elements will remain as zeros as allocated?
                    %Secondary note: Due to indexing here, dataSeq is indeed assembled in the 'Dinis' order
                
                % for the isomers (consumes a lot of memory)
                if ~skipIso
                    dataSeqIso(vol,seq,:, :, n+1) = imageStack(:,:,n*nVol + vol);
                end
            end
            %['vol:',num2str(vol),', max n:',num2str(n), 'last frame: ',num2str(n*nVol + vol)]
        end

        %Matt testatory plot to show mean transients for all seqs
        temp = nanmean( dataSeq , [3,4,5] );
        figure
        for seq = 1:size(dataSeq,2)
            subplot( ceil(sqrt(size(dataSeq,2))) , ceil(sqrt(size(dataSeq,2))), seq )
            plot( temp(:,seq) )
            title(['Seq #',num2str(seq)])
        end
        set(gcf,'Name','Mean seq transient')

        if skipIso
            disp([num2str(toc),'s to assemble dataSeq'])
        else
            disp([num2str(toc),'s to assemble dataSeq (and dataSeqIso)'])
        end
        %size(behavSeq)
    
    %Rolling, if applicable
    %elseif doRolling == 1
        %The below code is a mix of the original implementation of rollStruct and a modification for that to be the means to acquire rolling data
            %It has been disabled since the bulk of Dinis' 2p code assumes a blocklike nature, but in theory rollStruct could be revived
        %{
        tic
        disp(['~~ Sorting for rolling nature ~~'])
        disp(['Calculating rolling sequence'])
        %['imageStack size:',num2str(size(imageStack))]
        %['max projected ind:',num2str( (sequenceLength/nVol)*nVol + nVol )]
        %['max previous block n:',num2str( n )]
        %['max previous block ind:',num2str( n*nVol + vol )]
        %['sequence length: ',num2str(sequenceLength)]
        %['dataSeq size: ', num2str(size(dataSeq))]
        
        %kc4k
        
        %sizImgStack = 4800
        %sequenceLength = 4000

        %framesPer = sizImgStack / sequenceLength; %Not necessarily whole number
        framesPer = size(imageStack,3) / sequenceLength; %Not necessarily whole number
        disp(['Imaging frames/stimulus event: ',num2str(framesPer)])

        %blirg = framesPer:framesPer:sizImgStack; %For some reason, calculating like this gives float issues
        %frameIndices = 0:framesPer:sizImgStack;
        frameIndices = 0:framesPer:size(imageStack,3);
        frameIndices = frameIndices(2:end); %Nominally equal in size to sequenceLength
        %QA
        if size(frameIndices,2) ~= sequenceLength
            ['## Critical failure in frame indices calculations for rolling sequence! ##']
            crash = yes
        end
        %stimInds = nan( sequenceLength , ceil( sizImgStack / sequenceLength ) ); %In theory should be max number of frames per stim event
        minFrameCount = floor( framesPer );
        disp(['Min. frames/event: ',num2str(minFrameCount)])
        %QA
        if minFrameCount == 0
            ['-# Alert: Cannot calculate rolling with sub-1 volumes/stimulus event #-']
            crash = yes
        end

        dataSeq = zeros([minFrameCount nSeq/2 size(imageStack, [1,2]) sequenceLength ]); %Might be slightly too large, on account of skipping first nBack elements?
        seqSeq = zeros( minFrameCount , nSeq/2 , 1,1, sequenceLength  ); %QA tracker, for sequence correctness
        seqTrack = zeros(1,nSeq/2); %Similar to seqSeq, but simpler
        if ~skipIso
            dataSeqIso = zeros([minFrameCount nSeq size(imageStack, [1,2]) sequenceLength ]); %Even more crucial for this to be as small as possible...
            seqSeqIso = zeros( minFrameCount , nSeq , 1,1, sequenceLength  );
        else
            dataSeqIso = [];
            seqSeqIso = [];
        end


        %mirakuru
        
        %for vol = 1:nVol %'Vol' (i.e. Timepoint) only applicable if using block design?
            %e.g. 4800 frames / 6 'volumes' (i.e. timepoints) = 800 nonblank stimulus blocks
        %for n = 1:size(frameIndices)
        rollStruct.indTracker = ones( 1, nanmax(auxSeq) ); %Effectively a count of how many frames collected per sequence total
        rollStruct.framePos = cell(1,nanmax(auxSeq)); %For rollSeqReduced, list of frame positions of each element, split by sequence
        k = 0;
        tic
        for n = options.nBack+1:size(frameIndices,2) %Skip first nBack+1 elements (e.g. 5+1)
            %seq = bin2dec(num2str(randomSequence((n*nStimuli + 1):(n*nStimuli + options.nBack)))) + 1;
            seq = bin2dec(num2str(randomSequence(n-options.nBack:n-1))) + 1;
                %**Should** equate to the same as above (i.e. n = 6 (First loop) -> 6-5= 1   : 6-1= 5   
            %disp([num2str(k),' - ', num2str(seq)])
            %n
            %randomSequence(n-options.nBack:n-1)
            % stack images for each vol and seq along 5th dimension (separated by pattern)

            %%rollSeq(vol, auxSeq(seq),:, :, n) = imageStack(:,:,n*nVol + vol); %Note subtly different matrix insertion position
                %NOTE: MATH NOT VERIFIED FOR IMAGESTACK INDICES

            startEnd = [ frameIndices(n)-framesPer , frameIndices(n) ];
            rangi = [ceil(startEnd(1)):floor(startEnd(2))];
            framesToUse = rangi(end-minFrameCount+1:end);
            %disp( ['Stim event #',num2str(i),char(10),'start/stop: ', num2str( startEnd )] )
            %disp( ['Selected frames: ', num2str( framesToUse )] )

            dataSeq(:, auxSeq(seq),:, :, n-1) = permute( imageStack(:,:, [framesToUse] ) , [3,4,1,2] ); %'Extra' dim 4 in permute is cheat way to add singleton dim
                %Note n-1 positioning
            seqSeq( :, auxSeq(seq), 1,1, n-1 ) = repmat( auxSeq(seq), minFrameCount , 1 ); %Replicate for size time
            if ~skipIso
                dataSeqIso(:, seq,:, :, n-1) = permute( imageStack(:,:, [framesToUse] ) , [3,4,1,2] );
                seqSeqIso( :, seq, 1,1, n-1 ) = repmat( seq, minFrameCount , 1 );
            end
            seqTrack( auxSeq(seq) ) = seqTrack( auxSeq(seq) ) + 1;

            %%rollStruct.rollSeq(auxSeq(seq),:, :, [framesToUse]) = imageStack(:,:, [framesToUse] ); %Place at any random location
                %Disabled this, on account of inefficient memory use
            smartPos = [rollStruct.indTracker( auxSeq(seq) ) : rollStruct.indTracker( auxSeq(seq) )+minFrameCount-1 ]; %Note: Hardcoded always be min frame count 
            %smartPos
            rollStruct.indTracker( auxSeq(seq) ) = rollStruct.indTracker( auxSeq(seq) ) + minFrameCount;
            rollStruct.rollSeqReduced(auxSeq(seq),:, :, smartPos  ) = imageStack(:,:, [framesToUse] ); % Place in reduced manner
                %Verified math
            rollStruct.framePos{ auxSeq(seq) } = [ rollStruct.framePos{ auxSeq(seq) } , framesToUse ]; %Probably inefficient

            %%rollSeqWide = ??? %Intended to be a block-style version 
                %(i.e. Instead of collecting the min. number of frames between events, collect frames for nBack)
            k = k + 1;
        end
        %Trim rollSeq (if using reduced)
        %nansum( rollSeq( :,:,:, nanmax(indTracker):end ), 'all' ) %Quick calc to reveal if any actual data in the soon-to-be-reduced portion
        rollStruct.rollSeqReduced( :,:,:, nanmax(rollStruct.indTracker):end ) = [];
        %disp(['rollSeq reduced'])
        disp(['dataSeq (and rollSeqReduced) for rolling assembled in ',num2str(toc),'s'])
        
        %end 
        %disp([num2str(toc),'s to assemble rollSeq'])
        %disp([num2str(toc),'s to assemble rollSeqReduced'])
        %disp( ['size rollSeq: ',num2str(size(rollStruct.rollSeq)),' /reduced: ',num2str(size(rollStruct.rollSeqReduced))] )
        disp( ['size reduced: ',num2str(size(rollStruct.rollSeqReduced))] )
        %k
        disp( ['last seq #:',num2str(seq),' (#',num2str(auxSeq(seq)),' non-iso), startEnd: ',num2str(startEnd),', and framesToUse: ',num2str(framesToUse)] )
        disp( ['collected ', num2str(k-1),' events'] )
        %disp( ['seq dist.: ',num2str(rollStruct.indTracker/minFrameCount)] )
        disp( ['seq dist.: ',num2str(seqTrack)] )
        %Reduce rollSeq to only non-empty elements
        %hard
        %}
    %end
    %varien    
    
    
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
        d = 1;
        for unI = 1:length(unBehavs)
            disp(['Searching for behaviour state ',num2str(unBehavs(unI))])
            thisInds = find( behavSeq == unBehavs(unI) );
            disp(['Found ',num2str(length(thisInds)),' instances'])
            %dataSeqBehav(unI).dataSeq = dataSeq( :,:,:,:, thisInds  );
            thisInds = intersect( thisInds, nonEmptyInds ); %Reduce, since dataSeq is len 5000 unnecessarily
            disp(['(Reduced to ',num2str(length(thisInds)),' non-zero instances)'])
            %QA for empty
            if isempty( thisInds )
                ['-# State data empty; Deleting #-']
                continue
            end
            dataSeqReduced = dataSeq( :,:,:,:, thisInds );
            dataSeqReduced(isnan(dataSeqReduced)) = 0; %Just in case
            meanDataSeqReduced = mean( dataSeqReduced, 5 ); %Slightly different to how done in analyse2PBlock, but conceptually same?
            
            dataSeqBehav(d).state = unBehavs(unI); %What state this referred to
            dataSeqBehav(d).dataSeqReduced = dataSeqReduced;
            dataSeqBehav(d).meanDataSeqReduced = meanDataSeqReduced;
            dataSeqBehav(d).behavSeq = behavSeq; %Save information (with every state, but w/e)
            dataSeqBehav(d).numInstances = size(dataSeqBehav(d).dataSeqReduced,5); %Useful for later
            %disp([num2str(size(dataSeq)),' -> ',num2str(size(dataSeqBehav(unI).dataSeq))])
            disp([num2str(size(dataSeq)),' -> ',num2str( size(dataSeqBehav(d).dataSeqReduced) )])
            disp([ '(',num2str(dataSeqBehav(d).numInstances),' instances)' ])
            d = d + 1;
        end
        %disp([num2str(toc),'s to assemble dataSeqBehav for ',num2str(unI),' unique state/s'])
        disp([num2str(toc),'s to assemble dataSeqBehav for ',num2str(d-1),' unique state/s'])
    end
    
end