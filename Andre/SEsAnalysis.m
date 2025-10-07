% SEs analyis by Andre Spencer

% version+mark - date decommissioned, description
% V1...Mk0 - (2/9/25), script made, core functionality up to and including cross-fly analyses
% V1...Mk1 - (), sequence orders are now based on what side last stimulus was presented on rather than seq_eff_order (for side-sepped profiles - averaging pending) 

testing = 0;
sgFilter = 1; % option to apply a savitzky-golay filter to account for drift
doNorm = 1; %option to normalise flurorescence (against whole experiment average [for ROI])
baselineCorrect = 1; %optiont to perform a baseline correction - based on mean of blanks
individualPlots = 0; %option to create individual plots
ampSpecs = 1; % number of points to use when calculating max/min for amplitude
doSmooth = 1; % option to smooth data for better amplitude calcs
sideTraces = 1;
loopProfiles = 1;
fullProfiles = 0; %option to plot 16-history based profiles
limitProfiles = 1; %option to plot 4-hitsory based profiles (RRRR, AAAR, RRRA, AAAA)
useMaxOnly = 1;

%plotting options
doDots = 1;
doLines = 0;

ROI = 1;
regionName = 'leftMB'; %leftMB, background, 

%gathering or generating ROI maps
roiFlies = struct;
if ~ROI
    regionName = 'full field';    
else
    roiStruct = struct;
    roiCount = 1;

    %feels like bad practice to have near identical loop outside of real loop, but I want to select all ROIs before completing analysis
    for roiFly = 1:size(chosenFlies, 2)
        roiFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(roiFly),:);
        for roiBlock = 1:size(FLIES(roiFly).BLOCKS, 2)
            currentDate = char(datetime(roiFlyBlocks.Date(1),'Format','dMMMyy'));
            roiDirectory = [dataDirectory, currentDate, filesep, 'ROI'];
            % roiStructureDirectory = [roiDirectory, filesep, regionName];
            roiFileName = [FLIES(roiFly).BLOCKS(roiBlock).flyID, '_', regionName, '_roiSpecs.mat'];
            
            if exist([roiDirectory, filesep, roiFileName], 'file') == 2 
                load([roiDirectory, filesep, roiFileName]);
                roiMap = roiShortCut.roiMap;
                disp(['pre-existing ROI data loaded for fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum )])
            else
                disp('ROI data not found')
                data = cell(1);
                data{1} = nanmean(FLIES(roiFly).BLOCKS(roiBlock).greenChannel, 3);

                figure
                set(gcf, 'Name', ['ROI generation: fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum)], 'NUmberTitle', 'off' )
                % subplot(1,2,1)
                imagesc( nanmean( data{1} , 3) )
                title(['Mean for fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum)])
                %figure
                % subplot(1,2,2)
                % imagesc( nanstd( data{1} , [], 3) )
                % title(['SD'])
                %Make ROI
                disp(['Please draw polygon ROI on plot'])
                roy = drawpolygon;

                %Find coordinates
                allX = repmat( [1:size(data{1},2)]', size(data{1},2), 1 );
                allY = repmat( [1:size(data{1},2)], size(data{1},2), 1 );
                allY = reshape( allY, size(allY,1)*size(allY,2), 1 );

                %Find what inside
                [isInty] = isinterior( polyshape(roy.Position), allY, allX );
                isIntyRes = reshape( isInty, size(data{1},1) , size(data{1},2) );

                %Plot
                figure
                set(gcf, 'Name', ['ROI specs: fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum)], 'NUmberTitle', 'off' )
                imagesc( isIntyRes )

                %Rename
                roiMap = isIntyRes;
                roiPos = [NaN]; %Standin value to simplify lower loops/etc

                roiShortCut = struct;
                roiShortCut.roiMap = roiMap;

                if exist(roiDirectory, 'dir') ~= 7
                    mkdir(roiDirectory)
                    disp(['ROI file made at ', roiDirectory])
                end
                save([roiDirectory, filesep, roiFileName], 'roiShortCut')
                disp('ROI specs saved')

            end
            roiStruct(roiCount).roiMap = roiMap;
            roiCount = roiCount + 1;
        end
    end
    roiCount = 1; %resetting roiCount for use in main loop
end

% haha

%generating easy access cell for x-tick titles
% may be adapted in future to a selectable window size
SETitles = cell(1,size(seq_eff_order(5), 1));
% SETitles = [];
axisOrder = seq_eff_order(5) - 1; %5 hardcoded here (maybe blocklength? but could be interesting to add option of varying historylength
% for hist = 1:(size(historyTraces, 4) -1)/2 %-1 for blank, /2 for sidedness
for hist = 1:size(axisOrder, 1)
    currentSequence = dec2bin(axisOrder(hist), 5);
    thisHistory = [];
    for i = 2:size(currentSequence, 2)
        if currentSequence(i) == currentSequence(i-1)
            thisHistory = [thisHistory, 'R'];
        else
            thisHistory = [thisHistory, 'A'];
        end
    end
    SETitles{hist} = string(thisHistory);
end

%using seq_eff_order to generate side-separated dependeny orders
seq = double(dec2bin(seq_eff_order(5)-1, 5)) - 48; %sequence histories from seq_eff_order as a matrix of numbers
leftEnd = []; %will store sequences ending on left
for s = 1:size(seq, 1)
    currentSeq = seq(s,:);
    if currentSeq(end) ~= 0
        currentSeq(currentSeq == 1) = 2;
        currentSeq(currentSeq == 0) = 1;
        currentSeq(currentSeq == 2) = 0;
    end
    leftEnd = [leftEnd; currentSeq];
end

rightEnd = leftEnd; % will store sequences ending on the right
rightEnd(rightEnd == 1) = 2; % changing rights to an arbitrary number so they're position maintained
rightEnd(rightEnd == 0) = 1; % inverting left for right
rightEnd(rightEnd == 2) = 0; % re-inserting rights in correct syntax ('1')

%generating index
leftIndex = []; %will store index for to present sequences that end on the left
rightIndex = []; %^^ for ending on right
for h = 1:size(leftEnd, 1)
    leftBin = num2str(leftEnd(h,:)); %chose bin for 'binomial'
    leftBin(isspace(leftBin)) = [];
    leftIndex = [leftIndex, bin2dec(leftBin)];

    rightBin = num2str(rightEnd(h, :));
    rightBin(isspace(rightBin)) = [];
    rightIndex = [rightIndex, bin2dec(rightBin)];
end
leftIndex = leftIndex + 1;
rightIndex = rightIndex + 1;

SEStruct = struct;
structCount = 1;

%% volume collection

for fly = 1:size(FLIES, 2) %[1,2,5,7] - 'good' ATR flies?
    currentFly = FLIES(fly).BLOCKS;

    for b = 1:size(currentFly, 2)
        currentBlock = currentFly(b);

        %loading non-Asgardian data
        blockRecords = flyRecord(flyRecord.Fly == currentBlock.flyNum & flyRecord.Block == currentBlock.blockNum, :); %finding fly and specifying block in fly record
        % randomSequence = currentBlock.randomSequence;
        currentDate = char(datetime(blockRecords.Date(1),'Format','dMMMyy')); %dnt think this needs (1) index but left just in case
        currentDate_2 = char(datetime(blockRecords.Date(1),'Format','ddMMyy')); %2nd date uses different format
        flyID = [currentDate_2, '_f', num2str(blockRecords.FlyOnDay), '_b', num2str(blockRecords.Block)];

        btDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btData.csv'];
        btAncillaryDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btAncillary.mat'];
        matDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_MAT.mat'];

        %Loading files
        BTData = csvread(btDataDirectory);
        load(btAncillaryDirectory);
        load(matDataDirectory);

        %data gathering
        rawVols = currentBlock.greenChannel;
        volTimes = currentBlock.volTimes;
        volTimes(1,:) = volTimes(1,:) - volTimes(1,1);
        %maybe volTimes correction here - but need to verify
        sequenceLength = BTData(find(BTData(:,8) == -1, 1, 'last'), 5); %determinig sequencelength to use based on last full imaging period in BT
        eventSize = matSave.blockLength + matSave.pauseBlockLength + 1; % duration (in i) of stim and imaging periods , + 2 for leeway

        smoothPhot = smooth( FLIES(fly).BLOCKS(b).photData , 0.032*30000 );
        randomSequence = matSave.randomSequence;

        %applying ROI
        if ROI
            roiMap = roiStruct(roiCount).roiMap;
            for i = 1:size(rawVols, 3)
                thisIVols = squeeze(rawVols(:,:,i));
                thisIVols(~roiMap) = nan;
                rawVols(:,:,i) = thisIVols;
            end
        end

        %normalisation
        if doNorm
            f0 = nanmean(nanmean(rawVols, [1,2])); % nanmean('all')??
            rawVols = (rawVols - f0)/f0;
        end

        %smoothing (to aid amplitude selection)
        smoothVols = smooth(squeeze(nanmean(rawVols, [1,2])), 8); %8 currently subjective because using GCaMP6S at 6 steps(128x128)

        % data collection (all Realms)
        iTimes = []; %will contain time when each i-element changed (from BTData)
        iCount = 0;
        for l = 1:size(BTData(:, 6), 1)
            if BTData(l, 5) ~= iCount
                iTimes = [iTimes, BTData(l,6)];
                iCount = BTData(l, 5);
            end
        end
        %iTimes correction - to make first stimulus onset 0s (obtained from voltime row 2)
        iTimes = iTimes - BTData(volTimes(2,1), 6);

        %separating vols into i-fragments
         tic
         iCount = 1;
         lastVol = 0;
         volLength = 30; %standard length for 6-step 128x128 recordings (my default for battery/imaging)
         shiftAmount = 0;

         allIVols = nan(64, 64, volLength, sequenceLength); %[pixelX, pixelY, volumes, i]
         allIExtras = nan(3, volLength, sequenceLength); %same as allIVols but for time(1) photIndex(2) smoothedData(3) (separated by rows)
         if doSmooth
             allIExtras = nan(3, volLength, sequenceLength); %adding smoothed data to allIExtras as element 3
         end
         thisIVols = nan(64, 64, volLength, 2); % vols for a single i

         for vol = 1:size(volTimes, 2)
             if vol <= abs(shiftAmount) && shiftAmount < 0
                 thisIVols(:,:, vol - lastVol, 1) = nan(64, 64);
                 thisIVols(1,1, vol - lastVol, 2) = nan;
                 thisIVols(1,2, vol - lastVol, 2) = nan;
                 thisIVols(1,3, vol - lastVol, 2) = nan;
             else
                 thisIVols(:, :, vol - lastVol, 1) = rawVols(:, :, vol + shiftAmount); %actual vols
                 thisIVols(1,1, vol - lastVol, 2) = volTimes(1, vol + shiftAmount); %time
                 thisIVols(1,2, vol - lastVol, 2) = smoothPhot(volTimes(3, vol + shiftAmount)); %photodiode data
                 thisIVols(1,3, vol - lastVol, 2) = smoothVols(vol + shiftAmount); %smoothData
             end

             if iCount == sequenceLength
                 %something to collect last few vols
                 if vol == size(volTimes, 2)
                     continue %put in here because having issues with final i containing orders of magnitude more vols than every other i
                     allIVols( :, :, :, iCount) = thisIVols(:,:,:,1);
                     allIExtras(:, :, iCount) = thisIVols(1,[1,2],:,2);
                 end
             elseif volTimes(1, vol) >= iTimes(iCount+1) % could i+1 be changed to i+volAdjust?? <- Nope
                 % allIVols( :, :, :, 1, iCount ) = thisIVols(:, :, :, 1);
                 % allIVols( :, :, :, 2, iCount ) = thisIVols(:, :, :, 2);
                 allIVols(:, :, :, iCount ) = thisIVols(:,:,:,1); %maybe simpler/faster than above but unsure of functionality
                 allIExtras(:, :, iCount) = squeeze(thisIVols(1, [1:3], :, 2));
                 lastVol = vol;

                 % iCount = iCount + 1;
                 iCount = find(iTimes<=volTimes(1,vol), 1, 'last');
                 % iCount = find(iTimes>=volTimes(1,i), 1, 'first'); % FORBIDDEN! - don't use this line, it corrupts timing
                 thisIVols = nan(size(thisIVols));
             end
         end
         toc

         %chopping randomSequence into sequence order
         splicedRandomSequence = nan(2, eventSize, ceil(sequenceLength/(eventSize-1))); % randomSequence spliced into stim+imaging periods | [randomSequence(1)/i(2), length of stim+imaging in i, number of stim+imaging events]
         eventCount = 0; % in theory (if using BTData index, first val of loop should always be first event)
         prevFinalI = BTData( find( BTData(:,8)~=-1, 1, 'first'), 5) - 1;
        
         %using imaging spaces to splice
         for i = BTData(volTimes(2,1), 5):sequenceLength %(BTData(volTimes(2,1), 5) = volTimes offset
            % stimStyle = randomSequence(i);
            if randomSequence(i-1) == -1 && randomSequence(i) ~= -1
                eventCount = eventCount + 1;
                prevFinalI = i-1;
            end
            
            splicedRandomSequence(1, i - prevFinalI, eventCount) = randomSequence(i);
            splicedRandomSequence(2, i - prevFinalI, eventCount) = i;
         end

         %
         seqOrder = nan(1,size(splicedRandomSequence, 3));
         for seq = 1:size(splicedRandomSequence, 3)
            currentSeq = splicedRandomSequence(1, [1:matSave.blockLength], seq);
            currentSeq(currentSeq==2) = 0; %converting sequence into binary code
            if ismember(-2, currentSeq)
                seqOrder(seq) = -2;
            elseif isnan(currentSeq)
                % seqOrder(seq) = nan;
                continue
            else
                seqOrder(seq) = bin2dec(num2str(currentSeq));
            end
         end 

         seqOrder(isnan(seqOrder)) = [];
        
         % splicedRandomSequence(isnan(splicedRandomSequence)) = [];

         %using sliced ranomSequence to restitch allIVols
         eventTraces = nan(3, 40, size(splicedRandomSequence, 3));
         % eventTraces = nan(size(splicedRandomSequence));

         tic
         for eventNum = 1:size(splicedRandomSequence,3)
             if find(isnan(splicedRandomSequence(:,:,eventNum)), 1, 'first')
                 continue
             else
                 lastVol = 1;
                 for i = splicedRandomSequence(2,:,eventNum)
                     thisIVols = squeeze(nanmean(allIVols(:,:,:,i), [1,2]))';
                     currentSize = find(isnan(thisIVols), 1, 'first') - 1;

                     if doSmooth
                         eventTraces(1,lastVol:currentSize + lastVol - 1, eventNum) = allIExtras(3, 1:currentSize, i);
                     else
                         eventTraces(1,lastVol:currentSize + lastVol - 1, eventNum) = thisIVols(1:currentSize);
                     end
                     % eventTraces(2,lastVol:currentSize + lastVol - 1, eventNum) = repmat(randomSequence(i), 1, currentSize);
                     eventTraces(2,lastVol:currentSize + lastVol - 1, eventNum) = allIExtras( 2, 1:currentSize, i);
                     eventTraces(3,lastVol:currentSize + lastVol - 1, eventNum) = allIExtras( 1, 1:currentSize, i);

                     lastVol = lastVol + currentSize;
                 end
             end
         end
         toc

         % to plot traces on same timescale
         eventDurationTimes = nan(size(eventTraces, [3, 2]));
         for event = 1:size(eventTraces, 3)
            eventDurationTimes(event, :) = eventTraces(3,:,event) - eventTraces(3,1,event);
         end

         %applying filters
         filteredEventTraces = eventTraces;
         filteredEventTraces(1,:,:) = nan;
         sgFactor = sgolayfilt(squeeze(nanmean(rawVols, [1,2]))', 1, 59);
         filtIndex = 1;

         % savitzky-golay filter
         if sgFilter
             for e = 1:size(eventTraces, 3)
                 currentEvent = eventTraces(1,:,e);
                 for i = 1:find(isnan(currentEvent), 1, 'first') - 1
                     filteredEventTraces(1,i,e) = currentEvent(i) - sgFactor(filtIndex);
                     filtIndex = filtIndex+1;
                 end
             end
             eventTraces = filteredEventTraces;
         end

         % historyTraces = nan(size(eventTraces, 1), size(eventTraces, 2), ceil(size(seqOrder,2)/max(seqOrder)+1), max(seqOrder)+2); %+2 because 0 can't be an indice and blanks
         historyTraces = nan(size(eventTraces, 1) + 1, size(eventTraces, 2), ceil(size(seqOrder,2)/7), max(seqOrder)+2); % [datatype spec, vols, history repition, history type] datatypes = [calciumData(1), photodiode(2), rawTimes(3), durationTimes(4)] +2 because 0 can't be an indice and blanks
         traceNums = ones(1, size(historyTraces, 4));
         for eventNum = 1:size(seqOrder, 2)
             currentHistory = seqOrder(eventNum);
             if currentHistory == -2
                 historyTraces(1:size(eventTraces, 1),:, traceNums(end),end) = eventTraces(:,:,eventNum);
                 historyTraces(4, :, traceNums(end), end) = eventDurationTimes(eventNum, :);
                 traceNums(end) = traceNums(end) + 1;
             else
                 historyTraces(1:size(eventTraces, 1),:, traceNums(currentHistory + 1), currentHistory + 1) = eventTraces(:,:,eventNum);
                 historyTraces(4, :, traceNums(currentHistory+1), currentHistory+1) = eventDurationTimes(eventNum, :);
                 traceNums(currentHistory + 1) = traceNums(currentHistory + 1) + 1;
             end
             % traceNums(currentHistory) = traceNums(currentHistory) + 1;
         end

         %heaps of plots
         % if individualPlots
         %     for history = 1:(size(historyTraces, 4) -1) / 2 % (16) assumes will always be looking at 16 histories (i.e. 5-back)
         %         leftHistory = leftIndex(history);
         %         rightHistory = rightIndex(history);
         % 
         %         figure
         %         subplot(1, 2, [1])
         %         hold on
         %         plot(squeeze(historyTraces(4,:,:,leftHistory)), squeeze(historyTraces(1,:,:,leftHistory)))
         %         plot(nanmean(historyTraces(4,:,:,leftHistory), 3), nanmean(historyTraces(1,:,:,leftHistory), 3), 'LineWidth', 3, 'Color', 'k')
         %         title(dec2bin(leftHistory - 1, 5))
         % 
         %         subplot(1,2,[2])
         %         hold on
         %         plot(squeeze(historyTraces(4,:,:,rightHistory)), squeeze(historyTraces(1,:,:,rightHistory)))
         %         plot(nanmean(historyTraces(4,:,:,rightHistory), 3), nanmean(historyTraces(1,:,:,rightHistory), 3), 'LineWidth', 3, 'Color', 'k')
         %         title(dec2bin(rightHistory - 1, 5))
         % 
         %         set(gcf, 'Name', ['fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum)])
         %     end
         % end

         if sideTraces
             figure
             hold on
             plot( nanmean( historyTraces(4,:,:,leftIndex), [3,4]) , nanmean( historyTraces(1,:,:,leftIndex), [3,4] ), 'LineWidth', 3 )
             set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
             title(['fly ', num2str(currentFly(b).flyNum), ' block ', num2str(currentFly(b).blockNum), ' leftEnding histories'])

             figure
             hold on
             plot( nanmean( historyTraces(4,:,:,rightIndex), [3,4]) , nanmean( historyTraces(1,:,:,rightIndex), [3,4] ), 'LineWidth', 3)
             set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
             title(['fly ', num2str(currentFly(b).flyNum), ' block ', num2str(currentFly(b).blockNum), ' rightEnding histories'])
         end

         individualSideFeats = nan(size(historyTraces,4), 3); %[sequence history, amplitude(1)/means(2)/var(3)/max{beforeMin}(4)?] side-specific 
         for history = 1:size(individualSideFeats,1)
             currentTrace = squeeze(nanmean(historyTraces(1,:,:,history), 3));
             % peak = nanmean( maxk( currentTrace, ampSpecs ) );
             % trough = nanmean( mink( currentTrace, ampSpecs ) );

             peak = max( currentTrace );
             trough = min( currentTrace );

             if useMaxOnly
                 individualSideFeats(history, 1) = peak;
             else
                 individualSideFeats(history, 1) = peak - trough; %must be done on smoothed trace and probably only need to do single points for max/min
             end

             individualSideFeats(history, 2) = nanmean(currentTrace);

             individualSideFeats(history, 3) = nanstd(currentTrace);

         %     smoothedTrace = smooth(currentTrace, size(currentTrace, 2)/2);
         %     individualSideAmps(history, 3) = max(smoothedTrace) - min(smoothedTrace);
         end

         averageAmps = nan((size(historyTraces, 4)-1)/2, 2); %[sequence history, amplitude(1)/means(2)] side-averaged
         for history = 1:17 %hardcoded for nBack of 5
            % bothTraces = [squeeze( nanmean( historyTraces(1,:,:,history), 3)); squeeze( nanmean( historyTraces(1,:,:,size(historyTraces,4) - history), 3))]; %WRONG!!!
            if history == 17
                averageTrace = squeeze(nanmean(historyTraces(1,:,:, end), 3)); 
            else
                leftTraces = squeeze(historyTraces(1,:,:,leftIndex(history)));
                rightTraces = squeeze(historyTraces(1,:,:,rightIndex(history)));
                bothTraces = [leftTraces, rightTraces];
                averageTrace = nanmean(bothTraces);
            end

            %trimming - probably dnt need to trim
            % leftTrace(:, find(isnan(leftTrace(1,:)), 1, 'first'):end) = [];

            % averageTrace = nanmean(bothTraces);

            peak = nanmean( maxk( averageTrace, ampSpecs ) );
            trough = nanmean( mink( averageTrace, ampSpecs ) );

            if useMaxOnly
                averageAmps(history, 1) = peak;
            else
                averageAmps(history, 1) = peak - trough;
            end
            averageAmps(history, 2) = nanmean(averageTrace);
         end

         allTraceFeats = nan(size(historyTraces, 3), size(historyTraces, 4), 3); %amplitude of every indivudal trace [rep, history, amp(1)/mean(2)/var(3)]
         for hist = 1:size(historyTraces, 4) %- 1 %eveyrhting but blanks
             for rep = 1:size(historyTraces, 3)
                 if isnan(historyTraces(1, 1, rep, hist)) %skip empty reps
                     continue
                 end
                 currentAmp = max(historyTraces(1,:,rep,hist)) - min(historyTraces(1,:,rep,hist));
                 currentMean = nanmean(historyTraces(1,:,rep,hist));
                 currentVar = nanstd(historyTraces(1,:,rep,hist));

                 allTraceFeats(rep, hist, 1) = currentAmp;
                 allTraceFeats(rep, hist, 2) = currentMean;
                 allTraceFeats(rep, hist, 3) = currentVar;
             end
         end

         if loopProfiles
             %individual profiles (from allTrace features)
             allTraceSEMs = nanstd(allTraceFeats, [], 1)/sqrt(size(allTraceFeats, 1)); %little dodge in that assumes same number of reps for all histories but should be good for now
             %Amps
             %side1
             figure
             hold on
             plot(nanmean(allTraceFeats(:,leftIndex, 1), 1))
             errorbar([1:16],  nanmean(allTraceFeats(:,leftIndex, 1), 1), allTraceSEMs(1,leftIndex, 1), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
             line([1,16], [nanmean(allTraceFeats(:, 33, 1), 1), nanmean(allTraceFeats(:,33,1), 1)], 'LineWidth', 2, 'Color', 'r')
             xticks([1:16])
             xticklabels(SETitles)
             title(['fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum), ' leftEnding profile (amps)'])
             %spiffy fig sets
             set(gcf, 'Name', ['fly_', num2str(FLIES(fly).BLOCKS(b).flyNum), '_block_', num2str(FLIES(fly).BLOCKS(b).blockNum), '_left_ profile(amp)'])

             %side2
             figure
             hold on
             plot(nanmean(allTraceFeats(:,rightIndex, 1), 1))
             errorbar([1:16], nanmean(allTraceFeats(:,rightIndex, 1), 1), allTraceSEMs(:,rightIndex, 1), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
             line([1,16], [nanmean(allTraceFeats(:, 33, 1), 1), nanmean(allTraceFeats(:,33,1), 1)], 'LineWidth', 2, 'Color', 'r')
             xticks([1:16])
             xticklabels(SETitles)
             title(['fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum), ' rightEnding profile (amps)'])
             %spiffy fig sets
             set(gcf, 'Name', ['fly_', num2str(FLIES(fly).BLOCKS(b).flyNum), '_block_', num2str(FLIES(fly).BLOCKS(b).blockNum), '_right_ profile(amp)'])

             % %means
             % %side1
             % figure
             % hold on
             % plot(nanmean(allTraceFeats(:,leftIndex, 2), 1))
             % errorbar([1:16],  nanmean(allTraceFeats(:,leftIndex, 2), 1), allTraceSEMs(1,leftIndex, 2), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
             % line([1,16], [nanmean(allTraceFeats(:, 33, 2), 1), nanmean(allTraceFeats(:,33,2), 1)], 'LineWidth', 2, 'Color', 'r')
             % xticks([1:16])
             % xticklabels(SETitles)
             % title(['fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum), ' leftEnding profile (means)'])
             % %spiffy fig sets
             % set(gcf, 'Name', ['fly_', num2str(FLIES(fly).BLOCKS(b).flyNum), '_block_', num2str(FLIES(fly).BLOCKS(b).blockNum), '_left_ profile(mean)'])
             % 
             % %side2
             % figure
             % hold on
             % plot(nanmean(allTraceFeats(:,rightIndex, 2), 1))
             % errorbar([1:16], nanmean(allTraceFeats(:,rightIndex, 2), 1), allTraceSEMs(:,rightIndex, 2), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
             % line([1,16], [nanmean(allTraceFeats(:, 33, 2), 1), nanmean(allTraceFeats(:,33,2), 1)], 'LineWidth', 2, 'Color', 'r')
             % xticks([1:16])
             % xticklabels(SETitles)
             % title(['fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum), ' rightEnding profile (means)'])
             % %spiffy fig sets
             % set(gcf, 'Name', ['fly_', num2str(FLIES(fly).BLOCKS(b).flyNum), '_block_', num2str(FLIES(fly).BLOCKS(b).blockNum), '_right_ profile(mean)'])

             %variance
             %side1
             figure
             hold on
             plot(nanmean(allTraceFeats(:,leftIndex, 3), 1))
             errorbar([1:16],  nanmean(allTraceFeats(:,leftIndex, 3), 1), allTraceSEMs(1,leftIndex, 3), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
             line([1,16], [nanmean(allTraceFeats(:, 33, 3), 1), nanmean(allTraceFeats(:,33,3), 1)], 'LineWidth', 2, 'Color', 'r')
             xticks([1:16])
             xticklabels(SETitles)
             title(['fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum), ' leftEnding profile (var)'])
             %spiffy fig sets
             set(gcf, 'Name', ['fly_', num2str(FLIES(fly).BLOCKS(b).flyNum), '_block_', num2str(FLIES(fly).BLOCKS(b).blockNum), '_left_ profile(var)'])

             %side2
             figure
             hold on
             plot(nanmean(allTraceFeats(:,rightIndex, 3), 1))
             errorbar([1:16], nanmean(allTraceFeats(:,rightIndex, 3), 1), allTraceSEMs(:,rightIndex, 3), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
             line([1,16], [nanmean(allTraceFeats(:, 33, 3), 1), nanmean(allTraceFeats(:,33,3), 1)], 'LineWidth', 2, 'Color', 'r')
             xticks([1:16])
             xticklabels(SETitles)
             title(['fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum), ' rightEnding profile (var)'])
             %spiffy fig sets
             set(gcf, 'Name', ['fly_', num2str(FLIES(fly).BLOCKS(b).flyNum), '_block_', num2str(FLIES(fly).BLOCKS(b).blockNum), '_right_ profile(var)'])
         end

         %saving to struct for cross-fly analysis(es)
         SEStruct(structCount).fly = FLIES(fly).BLOCKS(b).flyNum;
         % SEStruct(structCount).block = b;
         if FLIES(fly).BLOCKS(b).flyNum == 123 %case for fly 123, where first block failed (every other block is 1 higher than reality)
             SEStruct(structCount).block = FLIES(fly).BLOCKS(b).blockNum - 1;
         else
             SEStruct(structCount).block = FLIES(fly).BLOCKS(b).blockNum;
         end
         SEStruct(structCount).historyTraces = historyTraces;
         SEStruct(structCount).sideSpecificFeats = individualSideFeats;
         SEStruct(structCount).averageAmps = averageAmps;

         disp(['data organisation complete for fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum)])
         structCount = structCount + 1;
    end
end

%% cross-fly analysis
% no
blockNums = [SEStruct.block];
% seq_eff_order

allFlyTraces = nan(size(FLIES, 2), size(historyTraces, 2), size(historyTraces, 1), size(historyTraces, 4), max(blockNums)); %[fly, vols, spec, history, block] size of these parameters for historyTraces doesn't change across loops
allFlyFeatures = nan(size(chosenFlies, 2), size(historyTraces, 4), size(SEStruct(1).sideSpecificFeats, 2), max(blockNums)); % [ fly, history, amps(1)/means(2)/var(3), block ]
allFlyAvgFeats = nan(size(allFlyFeatures, 1), size(SEStruct(1).averageAmps, 1), 1, size(allFlyFeatures, 4)); % [ fly, history, amps(1)/means(2)/var(3), block ]
currentFly = SEStruct(1).fly;
flyCount = 1;

for block = 1:size(blockNums, 2) %find(blockNums == b)
    % generating avg(crossFly) traces
    b = blockNums(block);
    if currentFly ~= SEStruct(block).fly
        currentFly = SEStruct(block).fly;
        flyCount = flyCount + 1;
    end

    %features
    allFlyFeatures(flyCount, :, :, b) = SEStruct(block).sideSpecificFeats;
    allFlyAvgFeats(flyCount, :, :, b) = SEStruct(block).averageAmps(:,1);

    %traces
    historyTraces = SEStruct(block).historyTraces;
    for dataType = 1:size(allFlyTraces, 3)
        allFlyTraces(flyCount, :, dataType, :, blockNums(block)) = nanmean( historyTraces(dataType,:,:,:), 3);
    end

end

if baselineCorrect
    for fly = 1:size(allFlyTraces, 1)
        for block = 1:size(allFlyTraces, 5)
            correctionFactor = nanmean(allFlyTraces(fly, :, 1, 33, block), 'all'); %correction factor is based on mean of blanks
            for history = 1:size(allFlyTraces, 4)
                allFlyTraces(fly, :, 1, history, block) = allFlyTraces(fly, :, 1, history, block) - correctionFactor;
            end
        end
    end
end

%generating SEM
featSEM = nanstd(allFlyFeatures, [], 1)/sqrt(size(allFlyFeatures, 1));
avgFeatSEM = nanstd(allFlyAvgFeats, [], 1)/sqrt(size(allFlyAvgFeats, 1));

%plotting all trace average for each fly/block (blanks exlcuded)
for block = 1:size(SEStruct, 2)
    currentHistories = SEStruct(block).historyTraces;
    flyName = SEStruct(block).fly;
    flyBlock = SEStruct(block).block;
    figure
    plot(nanmean(currentHistories(4,:,:,:), [3,4]), nanmean(currentHistories(1,:,:,:), [3,4]))
    title(['fly ', num2str(flyName), ' block ', num2str(flyBlock), ' all history average'])
end

%plotting all traces for each fly/block
for flyBlock = 1:size(SEStruct, 2)
    theseTraces = SEStruct(flyBlock).historyTraces;
    theseFeatures = SEStruct(flyBlock).sideSpecificFeats;
    figure
    for hist = 1:size(theseTraces, 4)
        subplot(4, 9, [hist])
        hold on
        plot(squeeze(theseTraces(4,:,:, hist)), squeeze(theseTraces(1,:,:,hist)))
        plot(nanmean(theseTraces(4,:,:,hist), 3), nanmean(theseTraces(1,:,:,hist), 3), 'LineWidth', 3, 'Color', 'k')

        if hist == 33 %redundant
            title(['blanks', newline, 'A=', num2str(theseFeatures(hist, 1)), ' M=', num2str(theseFeatures(hist, 2)), newline, ' V=', num2str(theseFeatures(hist, 3))])
        else
            title([num2str(dec2bin(hist-1, 5)), newline, 'A=', num2str(theseFeatures(hist, 1)), ' M=', num2str(theseFeatures(hist, 2)), newline, ' V=', num2str(theseFeatures(hist, 3))])
        end
    end
    set(gcf, 'Name', ['all traces fly ', num2str(SEStruct(flyBlock).fly), ' block ', num2str(SEStruct(flyBlock).block)])
end

%plotting averagetraces for each history
for block = 1:size(allFlyTraces, 5)
    figure
    for hist = 1:size(allFlyTraces, 4)
        subplot(4,9,[hist])
        hold on
        plot(squeeze(allFlyTraces(:,:,4,hist,block))', squeeze(allFlyTraces(:,:,1,hist,block))')
        plot(squeeze(nanmean(allFlyTraces(:,:,4,hist,block), 1)), squeeze(nanmean(allFlyTraces(:,:,1,hist,block), 1)), 'LineWidth', 3, 'Color', 'k' )
        set(gca, 'ylim', [-0.03, 0.03])

        if hist == 33
            title('blanks')
        else
            title([num2str(dec2bin(hist-1, 5))])
        end
    end
end

%sideSepped traces
for block = 1:size(allFlyTraces, 5)
    for side = [1,2]
        figure
        hold on
        for hist = 1:16
            if side == 1
                index = leftIndex(hist);
            else
                index = rightIndex(hist);
            end
            subplot(4,5,hist)
            hold on
            plot(squeeze(allFlyTraces(:,:,4,index,block))', squeeze(allFlyTraces(:,:,1,index,block))')
            plot(squeeze(nanmean(allFlyTraces(:,:,4,index,block), 1)), squeeze(nanmean(allFlyTraces(:,:,1,index,block), 1)), 'LineWidth', 3, 'Color', 'k' )
            set(gca, 'ylim', [-0.03, 0.03])
            title(num2str(dec2bin(leftIndex(hist) - 1, 5)))
        end
        subplot(4,5, 17)
        hold on
        plot(squeeze(allFlyTraces(:,:,4,33,block))', squeeze(allFlyTraces(:,:,1,33,block))')
        plot(squeeze(nanmean(allFlyTraces(:,:,4,33,block), 1)), squeeze(nanmean(allFlyTraces(:,:,1,33,block), 1)), 'LineWidth', 3, 'Color', 'k' )
        set(gca, 'ylim', [-0.03, 0.03])
        title('blanks')
        if side == 1
            set(gcf, 'Name', 'leftEnding')
        else
            set(gcf, 'Name', 'rightEnding')
        end
    end
end

%plotting RRRR(1) and RRRA(9) traces for each side - used for thesis
for block = 1:size(allFlyTraces, 5)
    for his = [1,9] %RRRR and RRRA index
        %leftEnding
        figure
        hold on
        plot(allFlyTraces(:,:, 4, leftIndex(his), block)', allFlyTraces(:,:,1,leftIndex(his), block)') %plot underliying individuals
        plot(squeeze(nanmean(allFlyTraces(:,:,4,leftIndex(his), block))), squeeze(nanmean(allFlyTraces(:,:,1,leftIndex(his), block))), 'LineWidth', 3, 'Color', 'k')
        %errorbar
        title(['leftEnding ', SETitles{his}, ' block ', num2str(block)])
        set(gca, 'ylim', [-0.02, 0.02])

        %rightEnding
        figure
        hold on
        plot(allFlyTraces(:,:, 4, rightIndex(his), block)', allFlyTraces(:,:,1,rightIndex(his), block)') %plot underliying individuals
        plot(squeeze(nanmean(allFlyTraces(:,:,4, rightIndex(his), block))), squeeze(nanmean(allFlyTraces(:,:,1,rightIndex(his), block))), 'LineWidth', 3, 'Color', 'k')
        %errorbar
        title(['rightEnding ', SETitles{his}, ' block ', num2str(block)])
        set(gca, 'ylim', [-0.02, 0.02])
    end
end

allTraceSEM = nanstd(allFlyTraces, [], 1)/sqrt(size(allFlyTraces, 1));
%plotting side-averaged RRRR and RRRA - for thesis(for real this time)
for block = 2:size(allFlyTraces, 5)
    for his = [1,8,9,16]
        averageTrace = nanmean(allFlyTraces(:,:,:, [leftIndex(his), rightIndex(his)], block), 4);
        % thisTraceSEM = nanstd(averageTrace, [], 1)/sqrt(size(averageTrace, 1));

        %trimming
        for f = 1:size(averageTrace, 1)
            averageTrace(f, find(averageTrace(f,:,4)>=3, 1, 'first'):end, 1) = nan;
        end
        thisTraceSEM = nanstd(averageTrace, [], 1)/sqrt(size(averageTrace, 1));

        figure
        hold on
        plot(averageTrace(:,:,4)', averageTrace(:,:,1)')
        plot(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), 'LineWidth', 3, 'Color', 'k')
        %errorbars
        errorbar(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), squeeze(thisTraceSEM(1,:,1)), 'LineStyle', 'none', 'Color', 'k')
        title(['side-averaged ', SETitles{his}, ' block ', num2str(block)])
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        set(gca, 'ylim', [-0.02, 0.02])
        fill([0 0 1 1], [-0.02, 0.02, 0.02, -0.02], 'b', 'FaceAlpha', 0.3)
    end

    % %averages on same plot
    % currentTraces = squeeze(nanmean(allFlytraces(:,:,:,[1,8,9,16],block), 1));
    % for f = 1:size(currentTraces, 1)
    %     currenTraces(f, find(currentTraces(f,:,4)>=3, 1, 'first'):end, 1) = nan;
    % end
    % figure
    % hold on
    % plot(squeeze(nanmean(allFlyTraces(:,:, 4, [1,8,9,16], block), 1)), squeeze(nanmean(allFlyTraces(:,:,1,[1,8,9,16],block), 1)) ) 
end

%plotting traces on same plots (really really for thesis this time)
SEColour = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]}; %4 colours for thesis
for block = 1:size(allFlyTraces, 5)
    %less surprising
    figure
    hold on
    colourNum = 1;
    for his = [1,16]
        averageTrace = nanmean(allFlyTraces(:,:,:, [leftIndex(his), rightIndex(his)], block), 4);

        %trimming
        for f = 1:size(averageTrace, 1)
            averageTrace(f, find(averageTrace(f,:,4)>=3, 1, 'first'):end, 1) = nan;
        end
        thisTraceSEM = nanstd(averageTrace, [], 1)/sqrt(size(averageTrace, 1));

        % plot(averageTrace(:,:,4)', averageTrace(:,:,1)')
        plot(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), 'LineWidth', 3, 'Color', SEColour{colourNum})
        %errorbars
        errorbar(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), squeeze(thisTraceSEM(1,:,1)), 'LineStyle', 'none', 'Color', SEColour{colourNum})
        % title(['side-averaged ', SETitles{his}, ' block ', num2str(block)])
        % set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        % set(gca, 'ylim', [-0.02, 0.02])
        % fill([0 0 1 1], [-0.02, 0.02, 0.02, -0.02], 'b', 'FaceAlpha', 0.3)
        colourNum = colourNum+1;
    end
     title(['side-averaged AAAA+RRRR block ', num2str(block)])
     set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
     set(gca, 'ylim', [-0.01, 0.015])
     fill([0 0 1 1], [-0.01, 0.015, 0.015, -0.01], 'b', 'FaceAlpha', 0.3)

    %most surprising
    figure
    hold on
    for his = [8,9]
        averageTrace = nanmean(allFlyTraces(:,:,:, [leftIndex(his), rightIndex(his)], block), 4);

        %trimming
        for f = 1:size(averageTrace, 1)
            averageTrace(f, find(averageTrace(f,:,4)>=3, 1, 'first'):end, 1) = nan;
        end
        thisTraceSEM = nanstd(averageTrace, [], 1)/sqrt(size(averageTrace, 1));

        % plot(averageTrace(:,:,4)', averageTrace(:,:,1)')
        plot(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), 'LineWidth', 3, 'Color', SEColour{colourNum})
        %errorbars
        errorbar(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), squeeze(thisTraceSEM(1,:,1)), 'LineStyle', 'none', 'Color', SEColour{colourNum})
        % title(['side-averaged ', SETitles{his}, ' block ', num2str(block)])
        % set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        % set(gca, 'ylim', [-0.02, 0.02])
        % fill([0 0 1 1], [-0.02, 0.02, 0.02, -0.02], 'b', 'FaceAlpha', 0.3)
        colourNum = colourNum+1;
    end
    title(['side-averaged AAAR+RRRA block ', num2str(block)])
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    set(gca, 'ylim', [-0.01, 0.015])
     fill([0 0 1 1], [-0.01, 0.015, 0.015, -0.01], 'b', 'FaceAlpha', 0.3)
end

%all on one plot (another that's really really for thesis)
for block = 2:size(allFlyTraces, 5)
    figure
    hold on
    colourNum = 1;
    for his = [1,16,8,9]
        averageTrace = nanmean(allFlyTraces(:,:,:, [leftIndex(his), rightIndex(his)], block), 4);

        %trimming
        for f = 1:size(averageTrace, 1)
            averageTrace(f, find(averageTrace(f,:,4)>=3, 1, 'first'):end, 1) = nan;
        end
        thisTraceSEM = nanstd(averageTrace, [], 1)/sqrt(size(averageTrace, 1));

        % plot(averageTrace(:,:,4)', averageTrace(:,:,1)')
        plot(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), 'LineWidth', 3, 'Color', SEColour{colourNum})
        %errorbars
        errorbar(squeeze(nanmean(averageTrace(:,:,4), 1)), squeeze(nanmean(averageTrace(:,:,1), 1)), squeeze(thisTraceSEM(1,:,1)), 'LineStyle', 'none', 'Color', SEColour{colourNum})
        colourNum = colourNum+1;
    end
    title(['side-averaged all lim traces block ', num2str(block)])
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    set(gca, 'ylim', [-0.01, 0.01])
    fill([0 0 1 1], [-0.01, 0.01, 0.01, -0.01], 'b', 'FaceAlpha', 0.3)
end

%full profile from averaged sides and features
for block = 2:size(allFlyAvgFeats, 4)
    %stats
    currentData = squeeze(allFlyAvgFeats(:,seq_eff_order(5), 1,block));
    %checking for normality
    nonNormCheck = 0;
    for g = 1:size(currentData, 2)
        normTest = lillietest(currentData(:,g));
        if ~normTest
            nonNormCheck = 1;
        end
    end

    if nonNormCheck
        disp('at least one history is not normally distributed, using nonparametric methods for stats')
        [~,~,SEStats] = kruskalwallis(currentData, [], 'off');
        SESig = multcompare(SEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
        statDescription = ['kruskal-wallis test', newline];
    else
        [~,~,SEStats] = anova1(currentData, [], 'off');
        SESig = multcompare(SEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
        statDescription = ['one-way anova', newline];
    end

    for group = 1:size(SESig, 1)
        statDescription = [statDescription, [num2str(SESig(group, 1)), 'v', num2str(SESig(group, 2)), ': ', num2str(SESig(group, 6)), newline]];
    end

    figure
    hold on
    bar(1:16, squeeze(nanmean(allFlyAvgFeats(:,seq_eff_order(5), 1,block), 1)))
    plot(1:16, squeeze(allFlyAvgFeats(:, seq_eff_order(5), 1, block)))
    % scatter(1:4, squeeze(allFlyAvgFeats(:, [1,8,9,16], 1, block)), 'filled')
    errorbar(1:16, squeeze(nanmean(allFlyAvgFeats(:,seq_eff_order(5), 1,block), 1)), avgFeatSEM(1, seq_eff_order(5), 1, block), 'LineStyle', 'none', 'Color', 'k')
    line([1,16], repmat(nanmean(allFlyAvgFeats(:,end,1,block), 1), 1, 2), 'LineWidth', 2, 'Color', 'r')
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    % figMax = get(gca, 'ylim');
    % figMax(1) = [];
    % text(barX(end), figMax/2, LRDescription)

    xticks([1:16])
    xticklabels(SETitles)

    figMax = get(gca, 'ylim');
    figMax(1) = [];
    text(17, figMax/2, statDescription)

    title(['crossFly sideAveraged, block ', num2str(block)])
    set(gcf, 'Renderer', 'painters')
end

if limitProfiles
    % for block = 1:size(allFlyFeatures, 4)
    %     for side = [1,2]
    %         if side == 1
    %             hisIndex = leftIndex([1,8,9,16]);
    %             description = 'leftEnding';
    %         else
    %             hisIndex = rightIndex([1,8,9,16]);
    %             description = 'rightEnding';
    %         end
    %         figure
    %         hold on
    %         % plot([1:4], squeeze(allFlyFeatures(:,hisIndex, 1, block)))
    %         bar([1:4], squeeze(nanmean(allFlyFeatures(:,hisIndex,1, block), 1)))
    %         plot([1:4], squeeze(allFlyFeatures(:,hisIndex, 1, block)))
    %         % plot([1:4], squeeze(nanmean(allFlyFeatures(:,hisIndex,1, block), 1)), 'LineWidth', 3, 'Color', 'k') %only plots amplitude
    %         errorbar([1:4], squeeze(nanmean(allFlyFeatures(:,hisIndex,1, block), 1)), featSEM(1, hisIndex, 1, block), 'LineStyle', 'none', 'Color', 'k')
    %         line([1,4], repmat(nanmean(allFlyFeatures(:,end,1,block), 1), 1,2), 'LineWidth', 3, 'Color', 'r')
    %         set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    %         xticks([1:4])
    %         xticklabels({SETitles{[1,8,9,16]}})
    %         title(['crossFly ', description, ' block ', num2str(block)])
    %         set(gcf, 'Renderer', 'painters')
    %     end
    % end

    % this one for thesis
    for block = 2%2:size(allFlyAvgFeats, 4)
        %stats
        currentData = squeeze(allFlyAvgFeats(:,[1,8,9,16, end], 1,block));
        %trimming currentData - abolishes fly identities, but necessary for paired statitistical analyses
        currentData(find(isnan(currentData(:,1))), :) = []; % only using first col as ref because all cols from same fly

        %checking for normality
        nonNormCheck = 0;
        for g = 1:size(currentData, 2)
            normTest = lillietest(currentData(:,g));
            if ~normTest
                nonNormCheck = 1;
            end
        end

        if nonNormCheck
            disp('at least one history is not normally distributed, using nonparametric methods for stats')
            [~,~,SEStats] = friedman(currentData, 1, 'off');
            SESig = multcompare(SEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            statDescription = ['Friedman test', newline];
        else
            [~,~,SEStats] = anova1(currentData, [], 'off');
            SESig = multcompare(SEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            statDescription = ['one-way anova', newline];
        end

        for group = 1:size(SESig, 1)
            statDescription = [statDescription, [num2str(SESig(group, 1)), 'v', num2str(SESig(group, 2)), ': ', num2str(SESig(group, 6)), newline]];
        end

        figure
        hold on
        bar(1:4, squeeze(nanmean(allFlyAvgFeats(:,[1,8,9,16], 1,block), 1)))
        if doDots
            scatter(1:4, squeeze(allFlyAvgFeats(:, [1,8,9,16], 1, block)), 10, 'k', 'filled')
        elseif doLines
            plot(1:4, squeeze(allFlyAvgFeats(:, [1,8,9,16], 1, block)))
        end
        % scatter(1:4, squeeze(allFlyAvgFeats(:, [1,8,9,16], 1, block)), 'filled')
        errorbar(1:4, squeeze(nanmean(allFlyAvgFeats(:,[1,8,9,16], 1,block), 1)), avgFeatSEM(1, [1,8,9,16], 1, block), 'LineStyle', 'none', 'Color', 'k')
        % line([1,4], repmat(nanmean(allFlyAvgFeats(:,end,1,block), 1), 1, 2), 'LineWidth', 2, 'Color', 'r')
        bar(5, squeeze(nanmean(allFlyAvgFeats(:,end,1,block), 1)))
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'ylim', [0, 0.03])
        if doDots
            scatter(5, allFlyAvgFeats(:,end,1,block), 10, 'k', 'filled')
        end
        errorbar(5, squeeze(nanmean(allFlyAvgFeats(:,end,1,block), 1)), avgFeatSEM(1,end,1,block), 'LineStyle', 'none', 'Color', 'k')
        % max = get(gca, 'ylim');
        % max(1) = [];
        % text(barX(end), max/2, LRDescription)

        % xticks([1:4])
        % xticklabels({SETitles{[1,8,9,16]}})
        xticks([1:5])
        xticklabels({SETitles{[1,8,9,16]}, 'Blanks'})
        
        figMax = get(gca, 'ylim');
        figMax(1) = [];
        text(4.5, figMax/2, statDescription)

        title(['crossFly sideAveraged, block ', num2str(block)])
        set(gcf, 'Renderer', 'painters')
    end
end

featureTitles = {'Amps', 'Means', 'Var'};
%generating SEs plot based on amplitudes (side-separated)
for flyBlock = 1:size(SEStruct, 2)
    for side = [1,2]
        for feat = [1,3] %cutting out means because doesn't work when savitsky-golayed
            figure
            % plot(SEStruct(flyBlock).averageAmps(seq_eff_order(5),1))
            if side==1
                plot(SEStruct(flyBlock).sideSpecificFeats(leftIndex, feat))
            else
                plot(SEStruct(flyBlock).sideSpecificFeats(rightIndex, feat))
            end
            line([1,16], [SEStruct(flyBlock).sideSpecificFeats(33,feat), SEStruct(flyBlock).sideSpecificFeats(33,feat)], 'Color', 'r', 'LineWidth', 2)
            % plot(SEStruct(flyBlock).averageAmps(seq_eff_order(5),2))
            % line([1,16], [SEStruct(flyBlock).sideSpecificAmps(33,2), SEStruct(flyBlock).sideSpecificAmps(33,2)], 'Color', 'r', 'LineWidth', 2)
            %xlabel = histories
            xticks([1:16])
            xticklabels(SETitles)

            if side ==1
                title(['fly ', num2str(SEStruct(flyBlock).fly), ' block ', num2str(SEStruct(flyBlock).block), ' left-start SEs (based on ', featureTitles{feat}, ')'])
            else
                title(['fly ', num2str(SEStruct(flyBlock).fly), ' block ', num2str(SEStruct(flyBlock).block), ' right-start SEs (based on ', featureTitles{feat}, ')'])
            end
            % title(['fly ', num2str(SEStruct(flyBlock).fly), ' block ', num2str(SEStruct(flyBlock).block), ' SEs (based on ', featureTitles{feat}, ')'])
        end
    end
end

%all trace plot(s)
for flyBlock = 1:size(SEStruct, 2)
    theseTraces = SEStruct(flyBlock).historyTraces;
    theseFeatures = SEStruct(flyBlock).sideSpecificFeats;
    figure
    for hist = 1:size(theseTraces, 4)
        subplot(4, 9, [hist])
        hold on
        plot(squeeze(theseTraces(4,:,:, hist)), squeeze(theseTraces(1,:,:,hist)))
        plot(nanmean(theseTraces(4,:,:,hist), 3), nanmean(theseTraces(1,:,:,hist), 3), 'LineWidth', 3, 'Color', 'k')

        if hist == 33 %redundant
            title(['blanks', newline, 'A=', num2str(theseFeatures(hist, 1)), ' M=', num2str(theseFeatures(hist, 2)), newline, ' V=', num2str(theseFeatures(hist, 3))])
        else
            title([num2str(dec2bin(hist-1, 5)), newline, 'A=', num2str(theseFeatures(hist, 1)), ' M=', num2str(theseFeatures(hist, 2)), newline, ' V=', num2str(theseFeatures(hist, 3))])
        end
    end
    set(gcf, 'Name', ['all traces fly ', num2str(SEStruct(flyBlock).fly), ' block ', num2str(SEStruct(flyBlock).block)])
end

%crossFly profiles (averaging features)
% featSEM = squeeze(nanstd(allFlyFeatures, [], 1)/sqrt(size(allFlyFeatures, 1)));

for block = 1:size(allFlyFeatures, 4)
    for side = [1,2]
        for feat = [1,3] %1:size(allFlyFeatures, 3)
            figure
            hold on
            if side==1
                plot(nanmean(allFlyFeatures(:,leftIndex, feat, block), 1))
                errorbar([1:16], nanmean(allFlyFeatures(:,leftIndex, feat, block), 1), featSEM(1, leftIndex, feat, block), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
            else
                plot(nanmean(allFlyFeatures(:,rightIndex, feat, block), 1))
                errorbar([1:16], nanmean(allFlyFeatures(:,rightIndex, feat, block), 1), featSEM(1, rightIndex, feat, block), 'LineStyle', 'none', 'Color', [0 0.4470 0.7410])
            end
            line([1,16], [nanmean(allFlyFeatures(:,33, feat, block), 1), nanmean(allFlyFeatures(:,33, feat, block), 1)], 'Color', 'r', 'LineWidth', 2)

            xticks([1:16])
            xticklabels(SETitles)

            if side ==1
                title(['cross-fly block ', num2str(block), ' leftEnding SEs (from averaged ', featureTitles{feat}, ')'])
            else
                title(['cross-fly block ', num2str(block), ' rightEnding SEs (from averaged ', featureTitles{feat}, ')'])
            end
        end
    end
end

%side averaging
% allFlyAvgFeats = nan(size(allFlyFeatures, 1), (size(allFlyFeatures, 2)-1)/2 + 1, size(allFlyFeatures, 3), size(allFlyFeatures, 4)); % [ fly, history, amps(1)/means(2)/var(3), block ] - I believe this will end up being in correct SEs hiostory order
% 
% for block = 1:size(allFlyAvgFeats, 4)
%     for his = 1:16
% 
%     end
% end

%crossFly profiles (using amplitudes of average traces for each history within a fly)
% crossFlyFeatures = nan(size(allFlyFeatures)); % [ fly, history, feature{amp(1)/max(2)}, block ]
% for block = 1:size(crossFlyFeatures, 4)
%     for fly = 1:size(crossFlyFeatures, 1)
%         for his = 1:size(allFlyTraces, 4)
%             currentTrace = allFlyTraces(fly, :, 1, his, block);
%             crossFlyFeatures(fly, his, 1, block) = max(currentTrace) - min(currentTrace);
%             crossFlyFeatures(fly, his, 2, block) = max(currentTrace);
%         end
%     end
% end
% 
% crossFlySEM = nanstd(crossFlyFeatures, [], 1)/sqrt(size(crossFlyFeatures, 1));

%calculating cross-fly amplitudes (side-averaged)
% crossFlyFeatures = nan(size(allFlyTraces, [2:5]));
% for blockType = 1:size(allFlyTraces, 5)
%     for block = find(blockNums == blockType)
%         for hist = 1:size(allFlyTraces, 4)
% 
%         end
%     end
% end

%% crossFly profiles from average (crossFly) trace for each history
%calculating amplitude (matches battery as much as possible)
SEAmps = nan(2, size(allFlyTraces, 4), size(allFlyTraces, 5)); % [ amp(1)/max(2), history, block ]
disp('amps restricted to outside imaging period and max before min')
for block = 1:size(allFlyTraces, 5)
    for hist = 1:size(allFlyTraces, 4)
        thisHistoryTrace = squeeze(nanmean(allFlyTraces(:,:,:,hist,block), 1)); % should become 40x4 matrix
        windowStart = find(thisHistoryTrace(:,4) >= 1, 1, 'first');
        windowEnd = find(~isnan(thisHistoryTrace(:,1)), 1, 'last');

        [~, minPos] = min(thisHistoryTrace(windowStart+3:windowEnd-3, 1));

        minPos = minPos + windowStart+3 - 1; %correcting for size change due to restricted window size

        [peakValues, peakIndex] = maxk(thisHistoryTrace(windowStart+3:minPos,1), ampSpecs);
        [troughValues, troughIndex] = mink(thisHistoryTrace(windowStart:windowEnd,1), ampSpecs);

        peak = nanmean(peakValues);
        trough = nanmean(troughValues);

        SEAmps(1,hist,block) = peak - trough; %amplitude
        SEAmps(2,hist,block) = peak; %max only
    end
end

%plotting side-seprated profiles for each block
for block = 1:size(SEAmps, 3)
    for side = [1,2]
        if side == 1
            index = leftIndex;
            sideName = 'leftEnding';
        else
            index = rightIndex;
            sideName = 'rightEnding';
        end
        figure
        hold on
        plot(SEAmps(1,index,block))
        xticks(1:16)
        xticklabels(SETitles)
        title([sideName, ' block ', num2str(block), ' from average multiFly traces'])
    end
end

%% cross block analysis

if size(allFlyAvgFeats, 4) < 3
    disp('not enough blocks specified')
else 
    SEOffset = [-0.222, 0, 0.222]; %hardcoded for 3 bars
    barPos = 1;
    figure 
    hold on
    for hist = [1,8, 9, 16]
        %stats
        currentData = squeeze(allFlyAvgFeats(:,hist,1,[2:4]));
        %trimming currentData - abolishes fly identities, but necessary for paired statitistical analyses
        currentData(find(isnan(currentData(:,1))), :) = []; % only using first col as ref because all cols from same fly

        nonNormCheck = 0;
        for col = 1:size(currentData, 2)
            normTest = lillietest(currentData(:,col));
            if ~normTest
                nonNormCheck = 1;
            end
        end

        if nonNormCheck
            disp(['a group was not normally distributed, Friedman test being used for stats'])
            [~,~,SEStats] = friedman(currentData, 1 , 'off');
            SESig = multcompare(SEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            SEDescription = ['Friedman test', newline];
        else
            [~,~,SEStats] = ranova(currentData, [] , 'off');
            SESig = multcompare(SEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            SEDescription = ['repeated-measures anova', newline];
        end

        for group = 1:size(SESig, 1)
            SEDescription = [SEDescription, [num2str(SESig(group,1)), 'v', num2str(SESig(group,2)), ': ', num2str(SESig(group,6)), newline]];
        end


        %creating profile
        bar(barPos, squeeze(nanmean(allFlyAvgFeats(:,hist, 1, [2:4]), 1)))
        barX = SEOffset + barPos;
        errorbar(barX, squeeze(nanmean(allFlyAvgFeats(:,hist, 1, [2:4]), 1)), squeeze(avgFeatSEM(:, hist, 1, [2:4])), 'LineStyle', 'none', 'Color', 'k')
        if doDots
            scatter(barX, squeeze(allFlyAvgFeats(:,hist,1,[2:4])), 10, 'k', 'filled')
        end
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

        figMax = get(gca, 'ylim');
        figMax(1) = [];
        text(barX(end), figMax*(3/4), SEDescription)

        title('crossBlock comparisons')
        set(gcf, 'Renderer', 'painters')

        barPos = barPos + 1;
    end
    xticks([1:4])
    xticklabels([SETitles{[1,8,9,16]}])
end

%% testing
if testing
    % plot cross-fly profile with underlying profiles
    % flyNames = cell(1, size(chosenFlies, 2));
    flies = [SEStruct(find(blockNums == 2)).fly];
    flyNames = cell(1, size(flies, 2));
    for f = 1:size(flyNames, 2)
        flyNames{f} = num2str(flies(f));
    end

    for block = 1:size(allFlyFeatures, 4)
        figure
        hold on
        plot([1:16], allFlyFeatures(:,leftIndex, 1, block), 'LineWidth', 2)
        plot([1:16], nanmean(allFlyFeatures(:,leftIndex, 1, block), 1), 'LineWidth', 3, 'Color', 'k')
        line([1,16], repmat(nanmean(allFlyFeatures(:,33,1, block), 1), 1,2), 'Color', 'r', 'LineWidth', 3)
        xticks(1:16)
        xticklabels(SETitles)
        legend([flyNames, {'mean'}, {'blank'}])
        title(['allFlies leftEnding - block ', num2str(block)])

        figure
        hold on
        plot([1:16], allFlyFeatures(:,rightIndex, 1, block), 'LineWidth', 2)
        plot([1:16], nanmean(allFlyFeatures(:,rightIndex, 1, block), 1), 'LineWidth', 3, 'Color', 'k')
        line([1,16], repmat(nanmean(allFlyFeatures(:,33,1, block), 1), 1,2), 'Color', 'r', 'LineWidth', 3)
        xticks(1:16)
        xticklabels(SETitles)
        legend([flyNames, {'mean'}, {'blank'}])
        title(['allFlies rightEnding - block ', num2str(block)])
    end

    %trying to limit the histories
    threeBackTitles = {'RRR', 'ARR', 'RAR', 'AAR', 'RRA', 'ARA', 'RAA', 'AAA'};
    twoBackTitles = {'RR', 'AR', 'RA', 'AA'};
    reorderedHisFeats = allFlyFeatures(:, [leftIndex, rightIndex, 33], :,:);

    threeBackFeats = nan(size(reorderedHisFeats, 1), 17, size(reorderedHisFeats, 3)); 
    %for block (only one block currently, so skipped for speed purposes)
    hisCount = 1;
    for his = 1:2:16
        threeBackFeats(:,hisCount, :) = nanmean(reorderedHisFeats(:, [his, his+1], :, :), 2);
        threeBackFeats(:,hisCount+8, :) = nanmean(reorderedHisFeats(:, [his+16, his+17], :, :), 2);

        hisCount = hisCount + 1;
    end

    %leftEnding 3-back
    figure
    hold on
    plot([1:8], threeBackFeats(:,1:8, 1), 'LineWidth', 2)
    plot([1:8], nanmean(threeBackFeats(:,1:8, 1), 1), 'LineWidth', 3, 'Color', 'k')
    line([1,8], repmat(nanmean(threeBackFeats(:,17, 1), 1), 1,2), 'Color', 'r', 'LineWidth', 3)
    xticks(1:8)
    xticklabels(threeBackTitles)
    legend([flyNames, {'mean'}, {'blank'}])
    title(['allFlies leftEnding - block ', num2str(block)])

    %rightEnding 3-back
    figure
    hold on
    plot([1:8], threeBackFeats(:,9:16, 1), 'LineWidth', 2)
    plot([1:8], nanmean(threeBackFeats(:,9:16, 1), 1), 'LineWidth', 3, 'Color', 'k')
    line([1,8], repmat(nanmean(threeBackFeats(:,17, 1), 1), 1,2), 'Color', 'r', 'LineWidth', 3)
    xticks(1:8)
    xticklabels(threeBackTitles)
    legend([flyNames, {'mean'}, {'blank'}])
    title(['allFlies rightEnding - block ', num2str(block)])


    %code to generate figure of i in time
    figure
    hold on
    iCount = 1;
    for l = 1:size(BTData(:, 6), 1)
        if BTData(l, 5) ~= iCount
            line([BTData(l, 6), BTData(l, 6)], [0, 1], 'Color', 'k')
            text(BTData(l,6), 0.1, ['i=', num2str(BTData(l, 5))], 'Color', 'r')
            iCount = BTData(l, 5);
        end
    end

    %time axis for full photodiode plotting
    fullPhotFittedTime = nan(1, size(FLIES(fly).BLOCKS(b).photData, 2));
    a = 1;
    for t = 1:size(volTimes, 2) - 1
        spacing = linspace(volTimes(1, t), volTimes(1, t+1), volTimes(3, t+1) - volTimes(3, t) );
        %fullPhotoTest(1, [find(isnan(fullPhotoTest), 1, 'first') : find(isnan(fullPhotoTest), 1, 'first') - 1 + size(spacing, 2)]) = spacing;
        fullPhotFittedTime( a:a+size(spacing,2)-1 ) =  spacing;
        a = a + size(spacing,2);
    end

    %recapitualying from allIVols
    first100 = nan(1, 3000);
    itest = nan(size(first100));
    lastVol = 1;
    for i = 1:100
        currentSize = find(isnan(squeeze(nanmean(allIVols(:,:,:,i), [1,2]))), 1, 'first') - 1;

        first100(1, lastVol:currentSize + lastVol - 1) = squeeze(nanmean(allIVols(:,:,1:currentSize, i), [1,2]))';
        lastVol = lastVol + currentSize + 1;
    end

    %testing filter application
    filteredEventTraces = eventTraces;
    filteredEventTraces(1,:,:) = nan;
    sgFactor = sgolayfilt(squeeze(nanmean(rawVols, [1,2]))', 1, 59);
    filtIndex = 1;

    for e = 1:size(eventTraces, 3)
        currentEvent = eventTraces(1,:,e);
        for i = 1:find(isnan(currentEvent), 1, 'first') - 1
            filteredEventTraces(1,i,e) = currentEvent(i) - sgFactor(filtIndex);
            filtIndex = filtIndex+1;
        end
    end

    %testing code to turn seq_eff_order into a title
    SETitles = cell(1,size(seq_eff_order(5), 1));
    % SETitles = [];
    axisOrder = seq_eff_order(5) - 1; %5 hardcoded here (maybe blocklength? but could be interesting to add option of varying historylength
    % for hist = 1:(size(historyTraces, 4) -1)/2 %-1 for blank, /2 for sidedness
    for hist = 1:size(axisOrder, 1)
        currentSequence = dec2bin(axisOrder(hist), 5);
        thisHistory = [];
        for i = 2:size(currentSequence, 2)
            if currentSequence(i) == currentSequence(i-1)
                thisHistory = [thisHistory, 'R'];
            else
                thisHistory = [thisHistory, 'A'];
            end
        end
        SETitles{hist} = string(thisHistory);
    end

    %creating editing seq_eff_order so ends on same side
    seq = double(dec2bin(seq_eff_order(5)-1, 5)) - 48;
    leftEnd = [];
    for s = 1:size(seq, 1)
        currentSeq = seq(s,:);
        if currentSeq(end) ~= 0
            currentSeq(currentSeq == 1) = 2;
            currentSeq(currentSeq == 0) = 1;
            currentSeq(currentSeq == 2) = 0;
        end
        leftEnd = [leftEnd; currentSeq];
    end

    rightEnd = leftEnd;
    rightEnd(rightEnd == 1) = 2; % changing rights to an arbitrary number so they're position maintained
    rightEnd(rightEnd == 0) = 1; % inverting left for right
    rightEnd(rightEnd == 2) = 0; % re-inserting rights in correct syntax ('1')

    %generating index
    leftIndex = [];
    rightIndex = [];
    for h = 1:size(leftEnd, 1)
        leftBin = num2str(leftEnd(h,:));
        leftBin(isspace(leftBin)) = [];
        leftIndex = [leftIndex, bin2dec(leftBin)];

        rightBin = num2str(rightEnd(h, :));
        rightBin(isspace(rightBin)) = [];
        rightIndex = [rightIndex, bin2dec(rightBin)];
    end
    leftIndex = leftIndex + 1;
    rightIndex = rightIndex + 1;

    %code to loop through transients based on sequence history number
    figure
    while true
        for a = 1:33
            clf
            hold on
            plot(squeeze(historyTraces(4,:,:,a)), squeeze(historyTraces(1,:,:,a)))
            plot(nanmean(historyTraces(4,:,:,a), 3), nanmean(historyTraces(1,:,:,a), 3), 'LineWidth', 3, 'Color', 'k')
            if a~=33
                title(num2str(a))
            else
                title('blanks')
            end
            drawnow
            pause(1.5)
        end
    end

    %saving to temporary datasets folder in results
    dataName = 'SEsATR';
    save(['C:\Users\s4749814\2p\2p\RESULTS\DataSets', filesep, dataName], 'allFlyAvgFeats', 'avgFeatSEM' )
end
