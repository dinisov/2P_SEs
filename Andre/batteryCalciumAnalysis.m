    % analysis script for left/right, luminosity and frequency conditions of battery paradigm

% version+mark - date decommissioned, description
% V1...Mk1 - (20/6/25), slightly muddled code, runs analysis for up to averaging across flies, with data gathering for rough timing of each condition (blanks exlusive outside of opto)
% V1...Mk2 - (7/7/25), unnecessary code (at end) removed, condition selection now based off time instead of i
% V1...Mk3 - (9/7/25), transient analysis added for other (non-opto) conditions, allIVols saving removed (potentially temporarily), time dimension added/filled for all conditions, volume shift option added, simple normalisation method added
% V1...Mk4 - (11/7/25), cross-fly ROIs added! (xDim/yDim now redundant) rescaled photodiode data added for plots, can make plots when no photodiode available
% V1...Mk5 - (11/7/25), ROI selection/saving moved to before main loop, main loop now runs based on FLIES.BLOCKS sizes (no longer fly record for nBlocks)
% V1...Mk6 - (15/7/25), ROI generation bit more informative for user, issue with normalistion for ROIs (thought to be) solved
% V1...Mk7 - (15/7/25), single-fly LR traces now plotted against time
% V1...Mk8 - (30/7/25), extras separated from allIVols to optimise and so allow high-speed recordings to be analysed, allIVols size increase for fast recordings automated
% V1...Mk9 - (1/8/25), multi-block LR trace analysis added, frequency collation and analysis added, easy access for video creation added (no saving yet)

% V2...Mk0 - (8/8/25), multiblock analysis for all conditions, amplitudes calculation avaliable for all conditions, baseline correction added (for cross-fly normalised ROI data), added option to select which graphs get generated
% V2...Mk1 - (29/8/25), LR traces split to lefts vs rights, ROI mask method added, blankVols analysis added
% V2...Mk2 - (2/9/25), additional loop added to allow cross ROI analyses within single script, extra dimension added to LRAmps to accomodate left vs right comparisons 
% V2...Mk3 - (10/9/25), more features added for amplitude calculation, cross ROI analyses added for LR amps, smoothed traces added to allIVols
% V2...Mk4 - (17/9/25), statistics applied for crossBlock analyses, statistical comparisons changed from unpaired to paired methods, statistical comparisons added for crossROI, amplitudeRestriction modified to calc max only before min, option added to not plot non-essential plots, AUC added as amplitude option
% V2...Mk5 - (), 

startTime = datetime('now');

%Analysis choices
%plotting
conditionChoice = [0, 1, 1, 1]; %[opto, LR, lum, freq] plot options
graphChoice = [1, 1, 1, 1]; %for each conditon, choice to plot: [detailed individual, individual, cross-fly, amplitude avgs]
figChoice = [0, 0, 0, 0]; %[rawData, LR, Lum, Freq]
periodChoice = 2; %1:blank period only, 2: full period(includes stimulus), 3: second 1/2 stim period + blank period

allGraphs = 0; %option to turn on plots unrelated to thesis
doDots = 1; % option to plot ugly and obtrusive but semi-informative dots.

%Analysis
doNormalisation = 1; % might expand for
fullFieldNorm = 1; % option to compare all ROIs using a consistent F0 (from the full-field)
baselineCorrect = 1; %option to shift baseline's by mean for each condition (only used in cross-fly analyses)
doSmooth = 1; %choice to smooth data (stored in allIExtras)

%amplitude
ampRestrict = 1; %choice to limit the time amps are collected from to post stimulus, also enforces that max value must occur before min during this period (currently hardcoded for 1s:end +-3vols)
windowSpecs = []; %times to use for amplitdue window calculation. ampRestrict must be on. leave empty for default (endStimulus,end)
AUC = 0; %option to represent amplitdues via area under the curve
LRBasedAmps = 0; %NOT YET FUNCTIONAL
ampSpecs = 1; %number of points to use to calculate peaks and troughs of amplitude
useMaxOnly = 1; %option to use only the max as indicator of amplitude - only valid because traces have been 0-meaned

%collection
doShift = 1; %choice to move shift timing of volumes with regards to i-timing - could potentially be automated using photodiode one day
shiftAmount = -3; % number of vols to shift by, negative to move right, positive for left. Only applicable if doShift is active.

useISequence = 0;
roughTiming = 1; %most likely temporary - used to analyse average opto stuff only. practicing cross fly averaging/getting a result quickly.
ignoreFreqSequence = 0; %choice to not use freqSequence in analysis (here primarily because I don't want to comment out code that uses blanks to calculate conditionInds)
analyseTransientOld = 0; %choice to analyse transients for each bat condition [non-functional], (vols not separated by i, transients separated roughly within loop, no saving)
analyseTransient = 1; % within each loop separates vols by i and into separate conditions which are saved for eventual cross-fly analyses
analyseAverage = 1; %choice to analyse average activity for each condition

%ROI specifications 
ROI = 1; %choice to specify region(s) of interest
% regionName = 'leftMB'; % ['EPG, EB, leftMB, rightMB, background] name of structure within ROI, used later in title and finding file
% roiList = {'bg'}; %list of regionNames, NOTE: 'm' for mask, 'f' for fullField (as additional ROI)
% roiList = {['leftMB(small)'], ['leftMB'], ['rightMB'], ['rightMB(small)'], ['EB'], ['bg']}; % list of regionNames, NOTE: 'm' for mask, 'f' for fullField (as additional ROI), recommend listing 'background' last
%roiList = {['leftMB'], ['rightMB'], ['EB'], ['bg']};
roiList = {['leftMB'], ['bg']};

if ~ROI
    % roiList = nan; %should eveluate to size 1 (single iteration through overarching loop) and will be clear that ROI not chosen
    roiList = {'full field'};
end
% doMask = 0;
%rightEPG, leftEPG, topLeftBG, bottomRightBG
%dimensions redundant since V1Mk4
% dimensionsROI = {[29, 23], [39, 32]}; %{[x,y], [x,y]}

%testing features
singleTesting = 0; %photodiode plots and calcs which only works for single block analysis

%QA
if useMaxOnly && ~baselineCorrect
    disp('WARNING: using max may not be valid as baseline correction not activated')
    crash = yes
end

if ampRestrict && (size(windowSpecs, 2) ~= 2 && ~isempty(windowSpecs))
    disp('window specification incorrect')
    crash = yes
end

overStruct = struct; %structure to allow cross-ROI analysese

saveVars = who; %saving all variables to this point
saveVars = [saveVars;'saveVars';'region']; %adding necessary variables to stored variables

for region = 1:size(roiList,2)
    clearvars('-except', saveVars{:} ); %resetting variables to what they were before loop

    regionName = roiList{region};
    if regionName == 'm' %option for mask method
        doMask = 1; %doMask being 1 should overwrite anything else and display that masking option is being used
    else
        doMask = 0;
    end

    % %ROI specifications -----!!--- will need to move if analysing multiple flies
    % ROI = 1; %choice to specify region(s) of interest
    % regionName = 'leftMB'; % ['EPG, leftMB, rightMB, background] name of structure within ROI, used later in title and finding file
    % doMask = 0;
    % %rightEPG, leftEPG, topLeftBG, bottomRightBG
    % %dimensions redundant since V1Mk4
    % dimensionsROI = {[29, 23], [39, 32]}; %{[x,y], [x,y]}

    %ROI setup
    roiFlies = struct;
    if doMask
        regionName = 'masked ROI';
    end
    % 
    % if ~ROI
    %     regionName = 'full field';

        %QA
        % if doMask
        %     disp(['warning: ROI option not activated, mask is being ignored'])
        % end
    % else
    if ROI
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
                    disp(['pre-existing ROI data loaded for fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum ), ' - ', regionName])
                else
                    if doMask %automatically generating ROI based on mask
                        disp(['Generating ROI using mask method'])
                        data = nanmean(FLIES(roiFly).BLOCKS(roiBlock).greenChannel, 3);
                        maskThreshold = nanmedian(data, 'all') + 2*nanstd(data, [], 'all'); %median + 2 standard deviations
                        roiMap = nan(64, 64);
                        roiMap(find(data > maskThreshold)) = 1;
                        roiMap(isnan(roiMap)) = 0;

                        %plotting
                        figure
                        imagesc(roiMap)
                    else %asking user to draw an ROI
                        disp('ROI data not found')
                        data = cell(1);
                        data{1} = nanmean(FLIES(roiFly).BLOCKS(roiBlock).greenChannel, 3);

                        figure
                        set(gcf, 'Name', ['ROI generation: fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum)], 'NUmberTitle', 'off' , 'Renderer', 'painters')
                        % subplot(1,2,1)
                        imagesc( nanmean( data{1} , 3) )
                        title(['Mean for fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum), ' - ', num2str(regionName)])
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
                        set(gcf, 'Name', ['ROI specs: fly ', num2str(chosenFlies(roiFly)), ' block ', num2str(FLIES(roiFly).BLOCKS(roiBlock).blockNum)], 'NumberTitle', 'off', 'Renderer', 'painters' )
                        imagesc( isIntyRes )

                        %Rename
                        roiMap = isIntyRes;
                        roiPos = [NaN]; %Standin value to simplify lower loops/etc
                    end

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

    % if ROI == 1
    %     xDim = [dimensionsROI{1}(1):dimensionsROI{2}(1)];
    %     yDim = [dimensionsROI{1}(2):dimensionsROI{2}(2)];
    % else
    %     %currently hardcoded for 64x64 pixel image size
    %     xDim = [1: 64];
    %     yDim = [1: 64];
    %     regionName = 'full field';
    % end

    %for saving, allowing cross fly analyses
    batStruct = struct; %used to store processed data
    structCount = 1; %used later to determine which row of batStruct.blocks to fill

    %QAs
    %ROIs
    % regionNameOptions = {'EPG'
    % if ROI && (regionName ~= 'EPG' || regionName ~= 'leftMB' || regionName ~= 'rightMB' || regionName ~=)
    %
    % end
    %shifting
    if doShift == 1 && shiftAmount == 0
        disp('WARNING: shift amount unspecified')
    end

    if doShift ~= 1 && shiftAmount ~= 0 %shift amount specified, but doShift not activated
        disp('WARNING: Volumes will not be shifted')
        shiftAmount = 0;
    end

    %% vol collection loop
    for fly = 1:length(chosenFlies)

        thisFly = chosenFlies(fly);
        thisFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(fly),:);

        % nBlocks = height(thisFlyBlocks);
        nBlocks = size(FLIES(fly).BLOCKS, 2);

        if nBlocks == 0
            ['-# Alert: No apparent blocks (for fly ',num2str(chosenFlies(fly)),') #-']
            %crash = yes
            continue
        end

        % for b = 1
        for b = 1:nBlocks
            currentBlock = thisFlyBlocks(b,:);

            if currentBlock.nStimuli == -2 || currentBlock.nStimuli == -5
                currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
                currentDate_2 = char(datetime(thisFlyBlocks.Date(1),'Format','ddMMyy')); %2nd date uses different format
                flyID = [currentDate_2, '_f', num2str(currentBlock.FlyOnDay), '_b', num2str(currentBlock.Block)];

                btDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btData.csv'];
                btAncillaryDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btAncillary.mat'];
                matDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_MAT.mat'];

                %Loading files
                BTData = csvread(btDataDirectory);
                load(btAncillaryDirectory);
                load(matDataDirectory);

                if size(BTData, 2) < 11
                    disp(['-_- WARNING: data pre-dates flip times from BT, analysis may be hindered or non-functional -_-'])
                    if isfield(matSave, 'batSequence') %feel like this isn't correct - find will give an array of the indices, not an integer
                        elementsLR = find(matSave.batSequence == 2);
                        elementsLum = find(matSave.batSequence == 3);
                        elemntsFreq = find(matSave.batSequence == 4);

                        LRInds = [matSave.optoElements + 1 : elementsLR];
                        lumInds = [elementsLR + 1 : elementsLum];
                        freqInds = [elementsLum + 1 : elementsFreq];
                    else
                        %calculating number of i values per stimulus condition
                        %currently assumes default battery order with no blanks
                        elementsLR = matSave.optoElements + matSave.numLR*2;
                        elementsLum = elementsLR + matSave.numLum*2;
                        elementsFreq = elementsLum + matSave.freqElements;
                        % randomSequence = matSave.randomSequence;
                        % freqSequence = matSave.freqSequence;

                        LRInds = [matSave.optoElements + 1 : elementsLR];
                        lumInds = [elementsLR + 1 : elementsLum];
                        freqInds = [elementsLum + 1 : elementsFreq];
                    end


                else
                    %generating frequency independent of freqSequence
                    flipI = [];
                    flipDurations = []; %duration of flip - derived from timeSinceLastFlip (BTData)
                    btFlipIndex = []; %last bt indice associated with the flip
                    onOff = [];
                    flipTime = 0;
                    % for i = find(BTData(:, 5)==2, 1, 'first'): size(BTData(:, 11), 1)
                    for i = 1:size(BTData(:, 11), 1)
                        if flipTime > BTData(i, 11)
                            flipDurations = [flipDurations, flipTime];
                            flipI = [flipI, BTData(i-1, 5)];
                            btFlipIndex = [btFlipIndex, i-1];
                            onOff = [onOff, BTData(i-1, 9)];
                        end
                        flipTime = BTData(i, 11);
                    end

                    %converts flipDurations to frequency (of each flip) - note: frequency of i (as reported by freqSequence) would be half these values
                    %also grabbing believed tBrightseq vals
                    flipFreq = [];
                    brightness = [];
                    for i = 1:size(flipDurations, 2)
                        freq = 1/flipDurations(i);
                        flipFreq = [flipFreq, freq];
                        brightness = [brightness, matSave.totalBrightSequence(flipI(i))];
                    end

                    timingMaster = [flipDurations; flipI; btFlipIndex; flipFreq/2; onOff; brightness];

                    %adding batSequence to timingMaster
                    batSync = [];
                    for i = 1:size(timingMaster, 2)
                        batSync = [batSync, matSave.batSequence(timingMaster(2, i))];
                    end
                    batSync(1) = 0; %disregarding first element because first blank empirically known to be at i=2
                    batSync(2) = -1; %adding blank back in
                    timingMaster = [timingMaster; batSync];

                    %                 blankLocations = timingMaster(2, find(flipDurations < 3 & flipDurations > 2)); % 5s blank spread across 2 flips
                    blankLocations = timingMaster(2, find(timingMaster(7, :) == -1));
                    %             end

                    %calculating number of i values per stimulus condition

                    elementsLR = matSave.optoElements + matSave.numLR*2;
                    elementsLum = elementsLR + matSave.numLum*2;
                    elementsFreq = elementsLum + matSave.freqElements;
                    randomSequence = matSave.randomSequence;

                    if ignoreFreqSequence
                        %QA
                        disp(['WARNING: freqSequence is being ignored, using blank locations to determine condition elements'])
                        %currently doesn't work correctly - reports LRInds[7201:7203], lumInds[7204:7604], freqInds seems to work

                        %finding battery order - only useful once I code in the modifiable order
                        conditionValue = -1;
                        batOrder = [];
                        for i = 1:size(matSave.batSequence, 2)
                            if matSave.batSequence(i) ~= conditionValue && matSave.batSequence(i) ~= -1
                                batOrder = [batOrder, matSave.batSequence(i)];
                                conditionValue = matSave.batSequence(i);
                            end
                        end
                        %QA
                        if size(batOrder, 2) < 4
                            disp(['-_- WARNING: not all battery conditions found for f:', num2str(currentBlock.Fly), ' b:', num2str(currentBlock.Block), ' -_-'])
                        end

                        %using blanks to create matrices containing indices of each condition
                        optoInds = [];
                        LRInds = [];
                        lumInds = [];
                        freqInds = [];
                        %blanks again hardcoded ugh
                        conditionCount = 1;
                        for i = 1:2:size(blankLocations, 2)
                            if i == 1
                                conditionStart = 3; %start at i=3 because initial blank
                            else
                                conditionStart = blankLocations(i - 1) + 1;
                            end

                            conditionEnd = blankLocations(i) - 1;
                            thisCondition = [conditionStart: conditionEnd];

                            if conditionCount <= size(batOrder, 2)
                                if batOrder(conditionCount) == 1
                                    optoInds = thisCondition;
                                elseif batOrder(conditionCount) == 2
                                    LRInds = thisCondition;
                                elseif batOrder(conditionCount) == 3
                                    lumInds = thisCondition;
                                else
                                    freqInds = thisCondition;
                                end
                            end
                            conditionCount = conditionCount + 1;
                        end
                    else
                        %                     optoInds = find(matSave.batSequence == 1);
                        %                     LRInds = find(matSave.batSequence == 2);
                        %                     lumInds = find(matSave.batSequence == 3);
                        %                     freqInds = find(matSave.batSequence == 4);
                        %                     blankInds = find(matSave.batSequence == -1);
                        %                     blankInds(end) = []; %remove random blank appended to end of batSequence
                        %%% ^^^ old news, assumes perfect timing which is wrong, below will use better but assumed timing

                        %                     optoInds = find(timingMaster(7, :) == 1);
                        %                     LRInds = find(timingMaster(7, :) == 2);
                        %                     lumInds = find(timingMaster(7, :) == 3);
                        %                     freqInds = find(timingMaster(7, :) == 4);
                        %                     blankInds = find(timingMaster(7, :) == -1);
                        %%% ^^ useful for sub-i timing, provides timingMaster index

                        optoInds = timingMaster(2, find(timingMaster(7, :) == 1));
                        LRInds = timingMaster(2, find(timingMaster(7, :) == 2));
                        lumInds = timingMaster(2, find(timingMaster(7, :) == 3));
                        freqInds = timingMaster(2, find(timingMaster(7, :) == 4));
                        blankInds = timingMaster(2, find(timingMaster(7, :) == -1));
                        blankInds = unique(blankInds); %removing double-ups of i caused by sub-i nature of timingMaster

                        freqSequence = matSave.freqSequence;
                    end
                end

                %             fudge

                %raw data configuration stuff
                iSequence = FLIES(fly).BLOCKS(b).iSequence;
                volTimes = FLIES(fly).BLOCKS(b).volTimes;
                volTimes(1, :) = volTimes(1,:) - volTimes(1,1); %timing correction for volTimes (validated by photodiode synchronization)
                rawVols = FLIES(fly).BLOCKS(b).greenChannel;
                volI = [];
                iValRaw = 1;

                % ROI Implementation
                if ROI
                    roiMap = roiStruct(roiCount).roiMap;
                    for i = 1:size(rawVols, 3)
                        thisIVols = squeeze(rawVols(:,:,i));
                        thisIVols(~roiMap) = nan;
                        rawVols(:,:,i) = thisIVols;
                    end
                end

                %using time rather than iSequence
                iTimes = []; %will contain time when each i-element changed (from BTData)
                iCount = 0;
                for l = 1:size(BTData(:, 6), 1)
                    if BTData(l, 5) ~= iCount
                        iTimes = [iTimes, BTData(l,6)];
                        iCount = BTData(l, 5);
                    end
                end

                % %making ROI if necessary
                % if ROI
                %     roiDirectory = [dataDirectory, filesep, currentDate, filesep, 'ROI'];
                %     roiStructureDirectory = [roiDirectory, filesep, regionName];
                %     roiFileName = [FLIES(fly).BLOCKS(b).flyID, '_roiShortCut.mat'];
                %     try
                %         load([roiDirectory, filesep, roiFileName]);
                %         roiMap = roiShortCut.roiMap;
                %         disp('pre-existing ROI data loaded')
                %     catch
                %         disp('ROI data not found')
                %         data = cell(1);
                %         data{1} = nanmean(FLIES(fly).BLOCKS(b).greenChannel, 3);
                %
                %         figure
                %         % subplot(1,2,1)
                %         imagesc( nanmean( data{1} , 3) )
                %         title(['Mean'])
                %         %figure
                %         % subplot(1,2,2)
                %         % imagesc( nanstd( data{1} , [], 3) )
                %         % title(['SD'])
                %         %Make ROI
                %         disp(['Please draw polygon ROI on plot'])
                %         roy = drawpolygon;
                %
                %         %Find coordinates
                %         allX = repmat( [1:size(data{1},2)]', size(data{1},2), 1 );
                %         allY = repmat( [1:size(data{1},2)], size(data{1},2), 1 );
                %         allY = reshape( allY, size(allY,1)*size(allY,2), 1 );
                %
                %         %Find what inside
                %         [isInty] = isinterior( polyshape(roy.Position), allY, allX );
                %         isIntyRes = reshape( isInty, size(data{1},1) , size(data{1},2) );
                %
                %         %Plot
                %         figure
                %         imagesc( isIntyRes )
                %
                %         %Rename
                %         roiMap = isIntyRes;
                %         roiPos = [NaN]; %Standin value to simplify lower loops/etc
                %
                %         roiShortCut = struct;
                %         roiShortCut.roiMap = roiMap;
                %
                %         if exist(roiDirectory) ~= 7
                %             mkdir(roiDirectory)
                %             disp('ROI file made at ', roiDirectory)
                %         end
                %         save([roiDirectory, filesep, roiFileName], 'roiShortCut')
                %         disp('ROI specs saved')
                %
                %     end
                % end

                %Normalising data [(F - F0)/F0] - F0 currently full-field + full recording avg
                %changes rawVols - anywhere greenChannel is used direclty will not be normalised!!
                if doNormalisation
                    % f0 = nanmean(nanmean(rawVols, [1,2]));
                    if fullFieldNorm
                        f0 = nanmean(FLIES(fly).BLOCKS(b).greenChannel, 'all');
                    else
                       f0 = nanmean(nanmean(rawVols, [1,2]));
                    end
                    rawVols = (rawVols - f0)/f0;
                end

                %smoothing data
                if doSmooth
                    smoothVols = smooth(squeeze(nanmean(rawVols, [1,2])), 8); %8 currently subjective because using GCaMP6S at 6 steps(128x128)
                end

                if isfield(FLIES(fly).BLOCKS(b), 'photData')
                    %photodiode
                    photData = FLIES(fly).BLOCKS(b).photData;
                    % smoothed photodiode
                    smoothPhot = smooth( FLIES(fly).BLOCKS(b).photData , 0.032*30000 );
                    %smoothPhot for sepcific times
                    timedPhoto = [];
                    for t = 1:size(volTimes, 2)
                        timedPhoto = [timedPhoto, smoothPhot(volTimes(3, t))];
                    end

                    %full photodiode data (all 5mill or so of it!)
                    fullPhotFittedTime = nan(1, size(FLIES(fly).BLOCKS(b).photData, 2));
                    a = 1;
                    for t = 1:size(volTimes, 2) - 1
                        spacing = linspace(volTimes(1, t), volTimes(1, t+1), volTimes(3, t+1) - volTimes(3, t) );
                        %fullPhotoTest(1, [find(isnan(fullPhotoTest), 1, 'first') : find(isnan(fullPhotoTest), 1, 'first') - 1 + size(spacing, 2)]) = spacing;
                        fullPhotFittedTime( a:a+size(spacing,2)-1 ) =  spacing;
                        a = a + size(spacing,2);
                    end

                    %photodiode for each i element, with 2nd element in 3rd dimension representing interpolated/linspaced time
                    % use iTimes
                    % % photSep = nan(matSave.sequenceLength, 100000, 2); %rows = i, cols = photodiode index/time, 3rd dim = determines what stored in cols
                    % % for i = 1:matSave.sequenceLength
                    % %     if i~= matSave.sequenceLength
                    % %         iStart = find( fullPhotFittedTime >= iTimes(i), 1, 'first');
                    % %         iEnd = find(fullPhotFittedTime <= iTimes(i+1), 1, 'last');
                    % %         photSep(i, :, 1) = smoothPhot()';
                    % %     end
                    % % end
                    %slow ^^^
                    %
                    % % tic
                    % % photSep = nan(matSave.sequenceLength + 1, 100000, 2); %rows = i, cols = photodiode index/time, 3rd dim = determines what stored in cols
                    % % splitPhot = nan(2, 100000);
                    % % iCount = 1;
                    % % lastPhotIndex = 0;
                    % % for i = 1:size(fullPhotFittedTime, 2)
                    % %     splitPhot(1, i - lastPhotIndex) = smoothPhot(i);
                    % %     splitPhot(2, i - lastPhotIndex) = fullPhotFittedTime(i);
                    % %
                    % %     if iCount == matSave.sequenceLength + 1
                    % %         if i == size(fullPhotFittedTime)
                    % %            photSep(iCount, 1:size(splitPhot, 2), 1) = splitPhot(1, :);
                    % %            photSep(iCount, 1:size(splitPhot, 2), 2) = splitPhot(2, :);
                    % %         end
                    % %     elseif fullPhotFittedTime(i) >= iTimes(iCount)
                    % %         photSep(iCount, 1:size(splitPhot, 2), 1) = splitPhot(1, :);
                    % %         photSep(iCount, 1:size(splitPhot, 2), 2) = splitPhot(2, :);
                    % %
                    % %         %resetting
                    % %         splitPhot = nan(2, 100000);
                    % %         iCount = iCount + 1;
                    % %         lastPhotIndex = i;
                    % %     end
                    % % end
                    % % toc
                    % photSep(1, :, :) = []; % remove first
                else
                    smoothPhot = nan(48000000,1);
                end


                %             for vol = 1:size(iSequence, 2)
                %                 if iSequence(vol) ~= iValRaw
                %                     volI = [volI, vol];
                %                     iValRaw = iSequence(vol);
                %                 end
                %             end
                %             beef


                %organising volumes
                if size(iSequence, 2) ~= size(volTimes, 2)
                    %                 volMaster = [iSequence', [volTimes(1, 1):volTimes(1, end - (size(volTimes, 2) - size(iSequence, 2)))]']; %assumes volTimes will always be longer than iSequence
                    %                 volMaster = [iSequence', [volTimes(1, 1):volTimes(1, end - 1)]'];
                    volTimes(:, size(iSequence, 2) + 1 : size(volTimes, 2)) = [];
                end
                volMaster = [iSequence', volTimes'];

                %code for adding approximate volume durations
                firstTest = [];
                lastTest = [];
                itest = [];

                iVolDurs = [];
                volDurs = [];
                volICount = 1;
                %             firstVolTime = volMaster(1, 2);
                %             lastVolTime = 0;

                for i = unique(volMaster(:, 1))'
                    firstVolTime = volMaster(find(volMaster(:, 1) == i, 1, 'first'), 2);
                    lastVolTime = volMaster(find(volMaster(:, 1) == i, 1, 'last'), 2);

                    volInfo = [(lastVolTime - firstVolTime); i];

                    iVolDurs = [iVolDurs, volInfo];
                end
                %
                %             iTime = [];
                %             for i = 1:size(iVolDurs, 2)
                %                 scriptedTime = sum( timingMaster(1, find( timingMaster( 2, :) == iVolDurs(2, i) ) ) );
                %                 iTime = [iTime, scriptedTime];
                %             end

                %             for i = 1:size(volMaster, 1)
                %                 if i ~= size(volMaster, 1)
                %                     volD = volMaster(i+1) - volMaster(i);
                %                     volDurs = [volDurs, volD];
                %                     if volMaster(i, 1) ~= volICount
                %                         iVolDurs = [iVolDurs, (lastVolTime - firstVolTime)];
                %                         volIcount = volMaster(i, 1);
                %                         firstVolTime = volMaster(i, 2);
                %                         firstTest = [firstTest, firstVolTime];
                %                         itest = [itest, volMaster(i, 1)];
                %                     else
                %                         lastVolTime = volMaster(i, 2);
                %                     end
                %                 else
                %                     volDurs = [volDurs, nan];
                %                     lastVolTime = volMaster(i, 2);
                %                     iVolDurs = [iVolDurs, (lastVolTime - firstVolTime)];
                %                 end
                %                 lastTest = [lastTest, lastVolTime];
                %             end

                %             volMaster = [volMaster; volDurs'];
                %^^not very useful^^ I want approximte durations for i NOT each volume

                %             transientsLR = nan(size(LRInds, 2) ,10 ); %10 arbitrary large number

                % %% add ROIs loop here?
                % if ROI
                %     roiMap = roiStruct(roiCount).roiMap;
                %     for i = 1:size(rawVols, 3)
                %         thisIVols = squeeze(rawVols(:,:,i));
                %         thisIVols(~roiMap) = nan;
                %         rawVols(:,:,i) = thisIVols;
                %     end
                % end

                %%%%%% POST TIMING STUFF %%%%%%%%%%%%%%%%%%%%%
                if roughTiming %will likely be removed, added so I can quickly access opto data exclusively doesn't factor in any analysis variables except for ROI, will likely go under if analyseAverage bool when ready to see other conditions
                    %% Opto, Blanks, allIVols
                    %generating indices used to grab blanks from iSequence
                    %                 blankIndex = [find(iSequence == blankInds(1), 1, 'first'): find(iSequence == blankInds(1), 1, 'last')]; %for opto condition, eventually maybe (probs not)
                    blankIndex = nan(5, 1000); %1000 cols arbitrary, first row will probably not be used, but would store opto, last row is for final blank (useful maybe for data pre-dating blanks)
                    conditionIndex = [];
                    for condition = 2:size(blankInds, 2) - 1 % -1 for now because not currently looking at post stimuli
                        if useISequence
                            conditionIndex = [find(iSequence == blankInds(condition), 1, 'first'): find(iSequence == blankInds(condition) + 1, 1, 'last')]; %uses iSequence
                        else
                            conditionIndex = [find(volTimes(1,:) >= iTimes(blankInds(condition)), 1, 'first'): find(volTimes(1,:) <= iTimes(blankInds(condition) + 2), 1, 'last')]; %using time to grab vols
                        end
                        blankIndex(condition, [1:size(conditionIndex, 2)]) = conditionIndex;
                    end

                    % total average for normalisation/comparison - similar to Kai methods
                    totalFullAvg = nanmean(rawVols, 3);

                    %first blank
                    %                 optoBlankVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, [find(iSequence == 2, 1, 'first'): find(iSequence == 2, 1, 'last')]); %hardcoded for now
                    if useISequence
                        optoBlankVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, [find(iSequence == blankInds(1), 1, 'first'): find(iSequence == blankInds(1), 1, 'last')]); %non-hardcoded, uses iSequence
                    else
                        % optoBlankVols = FLIES(fly).BLOCKS(b).greenChannel(yDim, xDim, [find(volTimes(1,:) >= iTimes(blankInds(1)), 1, 'first'): find(volTimes(1,:) <= iTimes(blankInds(1)+1), 1, 'last')]); %uses time
                        optoBlankVols = rawVols(:, :, [find(volTimes(1,:) >= iTimes(blankInds(1)), 1, 'first'): find(volTimes(1,:) <= iTimes(blankInds(1)+1), 1, 'last')]); %uses time
                    end

                    %                 optoBlankVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, blankIndices(1, :)); % non-hardcoded
                    %                 disp(['-_- optoblank indices has been hardcoded -_-'])
                    %             optoBlankFullAvg = nanmean( FLIES(fly).BLOCKS(b).greenChannel(:, :, [optoBlankVols(1): optoBlankVols(end)]), 3);
                    optoBlankFullAvg = nanmean(optoBlankVols, 3);
                    optoBlankFullSD = nanstd(optoBlankVols, [], 3);

                    optoBlankAvgVar = squeeze(nanmean(optoBlankFullSD, [1,2]));
                    optoBlankAvgVarSEM = nanstd(optoBlankFullSD, [], [1,2])/sqrt(size(optoBlankFullSD, 1)*size(optoBlankFullSD, 1));

                    optoBlankAvgs = [];
                    for vol = 1:size(optoBlankVols, 3)
                        volAvg = squeeze(nanmean(optoBlankVols(:, :, vol), [1,2]));
                        optoBlankAvgs = [optoBlankAvgs, volAvg];
                    end




                    % %                 figure
                    % %                 plot(FLIES(fly).BLOCKS(b).volTimes(1, [find(iSequence == 2, 1, 'first'): find(iSequence == 2, 1, 'last')]), optoBlankAvgs)
                    % %                 title('first blank transient')

                    %first 5 sec opto
                    %                 firstOptoVol = find(iSequence > 2, 1, 'first'); %assumes opto is first condition
                    if useISequence
                        firstOptoVol = find(iSequence >= optoInds(1), 1, 'first'); %non-hardcoded version
                    else
                        firstOptoVol = find( volTimes(1,:) >= iTimes(optoInds(1)), 1, 'first');
                    end
                    % first5OptoVols = FLIES(fly).BLOCKS(b).greenChannel(yDim, xDim, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 5)]);
                    first5OptoVols = rawVols(:, :, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 5)]);
                    first5OptoFullAvg = nanmean(first5OptoVols, 3);
                    first5OptoFullSD = nanstd(first5OptoVols, [], 3);

                    first5OptoAvgVar = squeeze(nanmean(first5OptoFullSD, [1,2]));
                    first5OptoAvgVarSEM = nanstd(first5OptoFullSD, [], [1,2]) / sqrt(size(first5OptoFullSD, 1)*size(first5OptoFullSD, 2));

                    first5OptoAvgs = [];
                    for vol = 1:size(first5OptoVols, 3)
                        volAvg = squeeze(nanmean(first5OptoVols(:, :, vol), [1,2]));
                        first5OptoAvgs = [first5OptoAvgs, volAvg];
                    end

                    %full opto
                    fullOptoVols = rawVols(:, :, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 120)]);
                    % fullOptoVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 120)]);
                    fullOptoAvg = [];
                    for vol = 1:size(fullOptoVols, 3)
                        volAvg = squeeze(nanmean(fullOptoVols(:, :, vol), [1,2]));
                        fullOptoAvg = [fullOptoAvg, volAvg];
                    end

                    fullOptoSD = nanstd(fullOptoVols, [], 3);
                    fullOptoVar = squeeze(nanmean(fullOptoSD, [1,2]));

                    % %                 figure
                    % %                 plot(FLIES(fly).BLOCKS(b).volTimes(1, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 5)]), first5OptoAvgs)
                    % %                 title('first 5 seconds optomotor transient')

                    %boxplotting opto things
                    SEMOptoBlank = nanstd(optoBlankAvgs)/sqrt(size(optoBlankAvgs, 2));
                    SEMFirst5Opto = nanstd(first5OptoAvgs)/sqrt(size(first5OptoAvgs, 2));

                    % %                 %optoBox = [optoBlankAvgs', first5OptoAvgs'];
                    % %                 figure
                    % %                 bar([1,2], [nanmean(optoBlankAvgs), nanmean(first5OptoAvgs)])
                    % %                 line([1, 1], [nanmean(optoBlankAvgs) - SEMOptoBlank, nanmean(optoBlankAvgs) + SEMOptoBlank], 'LineWidth', 2, 'Color', 'k')
                    % %                 line([2,2], [nanmean(first5OptoAvgs) - SEMFirst5Opto, nanmean(first5OptoAvgs) + SEMFirst5Opto], 'LineWidth', 2, 'Color', 'k')
                    % %
                    % %
                    % %                 figure
                    % %                 subplot( 1, 5, [1,2] )
                    % %                 image( optoBlankFullAvg*0.1 )
                    % %                 subplot( 1, 5, [3] )
                    % %                 bar([1,2], [nanmean(optoBlankAvgs), nanmean(first5OptoAvgs)])
                    % %                 line([1, 1], [nanmean(optoBlankAvgs) - SEMOptoBlank, nanmean(optoBlankAvgs) + SEMOptoBlank], 'LineWidth', 2, 'Color', 'k')
                    % %                 line([2,2], [nanmean(first5OptoAvgs) - SEMFirst5Opto, nanmean(first5OptoAvgs) + SEMFirst5Opto], 'LineWidth', 2, 'Color', 'k')
                    % %                 %             ylim([347,353])
                    % %                 ylim([min( [nanmean(optoBlankAvgs) - SEMOptoBlank, nanmean(first5OptoAvgs) - SEMFirst5Opto] ) - 1, max( [nanmean(optoBlankAvgs) + SEMOptoBlank, nanmean(first5OptoAvgs) + SEMFirst5Opto] + 1 )])
                    % %                 subplot( 1, 5, [4,5] )
                    % %                 image( first5OptoFullAvg*0.1 )
                    % %
                    % %                 %comparing SD images
                    % %                 figure
                    % %                 subplot(1,5,[1,2])
                    % %                 image( optoBlankFullSD )
                    % %                 subplot(1, 5, [3])
                    % %                 bar([1,2], [optoBlankAvgVar, first5OptoAvgVar])
                    % %                 line([1,1], [optoBlankAvgVar - optoBlankAvgVarSEM, optoBlankAvgVar + optoBlankAvgVarSEM], 'LineWidth', 2, 'Color', 'k')
                    % %                 line([2,2], [first5OptoAvgVar - first5OptoAvgVarSEM, first5OptoAvgVar + first5OptoAvgVarSEM], 'LineWidth', 2, 'Color', 'k')
                    % %                 subplot(1,5,[4,5])
                    % %                 image( first5OptoFullSD )
                    % %                 sgtitle(['average variance of blank and 5s opto'])

                    %stats for SD comparison
                    % % optoBlankVar = reshape(optoBlankFullSD, 4096, 1);
                    % % first5OptoVar = reshape(first5OptoFullSD, 4096, 1);
                    % % [h, pSDOpto] = ttest2(optoBlankVar, first5OptoVar);

                    %Grabbing all vols (using time)
                    %estimating volIndex associated with each i
                    % iTimesComp = repmat(iTimes', 1, size(volTimes, 2));
                    % volTimesComp = repmat(volTimes(1,:), size(iTimes, 2), 1);
                    % % iVsVolTimesComp = volTimesComp - iTimesComp;
                    % iVsVolTimesComp = iTimesComp - volTimesComp;

                    % % iVsVolTimesComp = repmat(volTimes(1,:), size(iTimes, 2), 1) - repmat(iTimes', 1, size(volTimes, 2));
                    % % [timeFromI, iVols] = min( abs(iVsVolTimesComp), [], 1);
                    % %
                    % % volIndex = nan(matSave.sequenceLength, 200);
                    % % for i = 1:size(volIndex, 1)
                    % %     vols = find( iVols == i);
                    % %     % vols = [vols(1)-1,vols,vols(end)+1];
                    % %
                    % %     %shifting based on percentage
                    % %     % vols = ceil(vols(1) + size(vols,2)*0.15) : floor( vols(end) + size(vols,2)*0.15 );
                    % %     % vols( vols < 1 ) = [];
                    % %     % vols( vols > size(rawVols,3) ) = [];
                    % %
                    % %     if ~isempty(vols)
                    % %         vols = ceil(vols(1) + size(vols,2)*0.4) : floor( vols(end) + size(vols,2)*0.4 );
                    % %         vols( vols < 1 ) = [];
                    % %         vols( vols > size(rawVols,3) ) = [];
                    % %
                    % %         volIndex(i, 1:size(vols, 2)) = vols;
                    % %     end
                    % % end

                    % % tic
                    % % allIVols = nan(64, 64, 200, 2, matSave.sequenceLength); %[pixelX, pixelY, volumes, imaging/other necessary info (time+photoData), i]
                    % % for i = 1:matSave.sequenceLength
                    % %     % vol = [volIndex(i, :) : find(~isnan(volIndex(i, :)), 1, 'last')];
                    % %     for vol = 1:size(volIndex, 2)
                    % %         if ~isnan(volIndex(i, vol))
                    % %             allIVols(:, :, vol, 1, i) = rawVols(:, :, volIndex(i, vol) );
                    % %             allIVols(1, 1, vol, 2, i) = volTimes(1, volIndex(i, vol));
                    % %             allIVols(1, 2, vol, 2, i) = smoothPhot( volTimes(3, volIndex(i, vol) ) );
                    % %         else
                    % %             break %speeding up loop
                    % %         end
                    % %     end
                    % % end
                    % % toc

                    disp('generating allIVols') %change eventually as only makes sense to me
                    tic
                    iCount = 1;
                    lastVol = 0;
                    volLength = 200; %standard length for 6-step 128x128 recordings (my default for battery/imaging)
                    if size(FLIES(fly).BLOCKS(b).greenChannel, 3) > 100000
                        volLength = 800; %could be a smarter way to get this number
                        disp(['large volume number detected. allIVols size increased for fly ', num2str(FLIES(fly).BLOCKS(b).flyNum), ' block ', num2str(FLIES(fly).BLOCKS(b).blockNum)])
                        %currently assumes all large battery recordings will be due to faster imaging
                    end

                    % allIVols = nan(64, 64, 200, 2, matSave.sequenceLength); %[pixelX, pixelY, volumes, imaging/other necessary info (time+photoData), i]
                    allIVols = nan(64, 64, volLength, matSave.sequenceLength); %[pixelX, pixelY, volumes, i]
                    allIExtras = nan(1, 2, volLength, matSave.sequenceLength); %same as allIVols but for time and photIndex (separated by dimension 2)
                    if doSmooth
                        allIExtras = nan(1, 3, volLength, matSave.sequenceLength);
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
                            thisIVols(1,3, vol - lastVol, 2) = smoothVols(vol + shiftAmount); % smoothed vols
                        end

                        if iCount == matSave.sequenceLength
                            %something to collect last few vols
                            if vol == size(volTimes, 2)
                                allIVols( :, :, :, iCount) = thisIVols(:,:,:,1);
                                allIExtras( :, :, :, iCount) = thisIVols(1,[1:size(allIExtras, 2)],:,2);
                            end
                        elseif volTimes(1, vol) >= iTimes(iCount+1) % could i+1 be changed to i+volAdjust?? <- Nope
                            % allIVols( :, :, :, 1, iCount ) = thisIVols(:, :, :, 1);
                            % allIVols( :, :, :, 2, iCount ) = thisIVols(:, :, :, 2);
                            allIVols(:, :, :, iCount ) = thisIVols(:,:,:,1); %maybe simpler/faster than above but unsure of functionality
                            allIExtras( :, :, :, iCount) = thisIVols(1, [1:size(allIExtras, 2)], :, 2);
                            lastVol = vol;

                            % iCount = iCount + 1;
                            iCount = find(iTimes<=volTimes(1,vol), 1, 'last');
                            % iCount = find(iTimes>=volTimes(1,i), 1, 'first'); % FORBIDDEN! - don't use this line, it corrupts timing
                            thisIVols = nan(size(thisIVols));
                        end
                    end
                    toc

                    % boom
                    %test for number of vols put into each i-element of allIVols
                    % volNumTest = nan(1, size(allIVols, 5));
                    % for i = 1:size(allIVols, 5)
                    %     for j = 1:size(allIVols, 3)
                    %         if isnan(allIVols(1,1,j,1,i) )
                    %             volNumTest(i) = j;
                    %             break
                    %         end
                    %     end
                    % end

                    %non-opto things
                    %LR
                    LRBlankVols = rawVols(:, :, blankIndex(2, [1: find(~isnan(blankIndex(2, :)), 1, 'last')]) );
                    % LRBlankVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, blankIndex(2, [1: find(~isnan(blankIndex(2, :)), 1, 'last')]) );
                    LRBlankAvg = squeeze(nanmean(LRBlankVols, [1,2]));


                    %Lum
                    lumBlankVols = rawVols(:, :, blankIndex(3, [1: find(~isnan(blankIndex(3, :)), 1, 'last')]) );
                    % lumBlankVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, blankIndex(3, [1: find(~isnan(blankIndex(3, :)), 1, 'last')]) );
                    lumBlankAvg = squeeze(nanmean(lumBlankVols, [1,2]));

                    %Freq
                    freqBlankVols = rawVols(:, :, blankIndex(4, [1: find(~isnan(blankIndex(4, :)), 1, 'last')]) );
                    % freqBlankVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, blankIndex(4, [1: find(~isnan(blankIndex(4, :)), 1, 'last')]) );
                    freqBlankAvg = squeeze(nanmean(freqBlankVols, [1,2]));

                    %saving for cross fly analyses
                    % disp(['Analysis complete for fly ', num2str(chosenFlies(fly)), ' block', num2str(b)])
                    batStruct.blocks(structCount).fly = thisFly;
                    % batStruct.blocks(structCount).block = b; %may not be correct if not all blocks analysed
                    if thisFly == 101 %case for fly 101, where first block failed
                        batStruct.blocks(structCount).block = FLIES(fly).BLOCKS(b).blockNum - 1; 
                    else
                        batStruct.blocks(structCount).block = FLIES(fly).BLOCKS(b).blockNum;
                    end
                    batStruct.blocks(structCount).wholeExpAvg = totalFullAvg;

                    %opto
                    batStruct.blocks(structCount).blankRaw = optoBlankVols;
                    if useISequence
                        batStruct.blocks(structCount).blankRawTimes = volTimes(1, [find(iSequence == blankInds(1), 1, 'first'): find(iSequence == blankInds(1), 1, 'last')]);
                    else
                        batStruct.blocks(structCount).blankRawTimes = volTimes(1, [find(volTimes(1,:) >= iTimes(blankInds(1)), 1, 'first'): find(volTimes(1,:) <= iTimes(blankInds(1)+1), 1, 'last')]);
                    end

                    batStruct.blocks(structCount).optoRaw = first5OptoVols;
                    batStruct.blocks(structCount).optoRawTimes = volTimes(1, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 5)]);
                    batStruct.blocks(structCount).fullOptoRaw = fullOptoVols;
                    batStruct.blocks(structCount).fullOptoRawTimes = volTimes(1, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 120)]);

                    batStruct.blocks(structCount).blankAvg = nanmean(optoBlankAvgs);
                    batStruct.blocks(structCount).blankVar = optoBlankAvgVar;

                    %                 batStruct.blocks(structCount).optoRaw = first5OptoVols;

                    %                 batStruct.blocks(structCount).optoAvg = nanmean(fullOptoAvg);
                    %                 batStruct.blocks(structCount).optoAvgFiveS = nanmean(first5OptoAvgs);
                    batStruct.blocks(structCount).optoAvg = nanmean(first5OptoAvgs);
                    batStruct.blocks(structCount).optoVar = first5OptoAvgVar;

                    batStruct.blocks(structCount).fullOptoAvg = nanmean(fullOptoAvg);
                    batStruct.blocks(structCount).fullOptoVar = fullOptoVar;

                    %LR
                    batStruct.blocks(structCount).LRBlankRaw = LRBlankVols;
                    batStruct.blocks(structCount).LRBlankRawTimes = volTimes(1, blankIndex(2, [1: find(~isnan(blankIndex(2, :)), 1, 'last')]) );
                    batStruct.blocks(structCount).LRBlankAvg = nanmean(LRBlankAvg);

                    %Lum
                    batStruct.blocks(structCount).lumBlankRaw = lumBlankVols;
                    batStruct.blocks(structCount).lumBlankRawTimes = volTimes(1, blankIndex(3, [1: find(~isnan(blankIndex(3, :)), 1, 'last')]) );
                    batStruct.blocks(structCount).lumBlankAvg = nanmean(lumBlankAvg);

                    %everything, everywhere, all at once
                    % % batStruct.blocks(structCount).iSortedVols = allIVols;

                    %Freq
                    batStruct.blocks(structCount).freqBlankRaw = freqBlankVols;
                    batStruct.blocks(structCount).freqBlankRawTimes = volTimes(1, blankIndex(4, [1: find(~isnan(blankIndex(4, :)), 1, 'last')]) );
                    batStruct.blocks(structCount).freqBlankAvg = nanmean(freqBlankAvg);

                    %code for grabbing photoDiode data from a specific time window
                    % for data = [blankRawTimes, optoTimes, ... ] %example of how multiple dataTypes can be looped
                    % for x = 1:size(batStruct.blocks(1).blankRawTimes, 2)
                    %     photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).blankRawTimes(x) ) );
                    %     photo = FLIES.BLOCKS(b).photData(photoIndex);
                    %     blankPhoto = [blankPhoto, photo];
                    % end
                    % end

                    %% LR traces
                    disp('LR vols being collected')
                    LRVols = nan(size(unique(LRInds), 2)/2, size(allIVols, 3), 3); %maybe 3rd dim for time and photo?
                    LRImages = nan(64, 64, size(allIVols, 3), size(unique(LRInds), 2)/2);
                    LRCount = 1;

                    for i = LRInds(1) + 1:2:LRInds(end)+1 %7204:2:7603 %LRInds hardcoded [83]
                        % for i = 7204:2:7243 % [87]
                        % stimVol = allIVols(yDim, xDim, :, 1, i);
                        stimVolImages = allIVols(:, :, :, i);
                        %disp([num2str(LRCount),' - ',num2str(stimVol(1,1,1))])
                        stimPhot = allIExtras(1, 2, :, i);
                        stimTime = squeeze(allIExtras(1,1,:,i));
                        % stimVol(isnan(stimVol(:,))) = [];
                        % nonStimVol = allIVols(yDim, xDim, :, 1, i+1);
                        nonStimVolImages = allIVols(:, :, :, i+1);
                        nonStimPhot = allIExtras(1, 2, :, i+1);
                        nonStimTime = squeeze(allIExtras(1,1,:,i+1));
                        % nonstimVol(isnan(nonstimVol)) = [];

                        %averaging to get trace
                        if doSmooth
                            stimVol = squeeze(allIExtras(1,3,:,i));
                            nonStimVol = squeeze(allIExtras(1,3,:,i+1));
                        else
                            stimVol = squeeze(nanmean(stimVolImages, [1,2]));
                            nonStimVol = squeeze(nanmean(nonStimVolImages, [1,2]));
                        end

                        stimPhot = squeeze(nanmean(stimPhot, [1,2]));
                        nonStimPhot = squeeze(nanmean(nonStimPhot, [1,2]));

                        %correcting after squeeze
                        stimVol = stimVol';
                        nonStimVol = nonStimVol';

                        stimPhot = stimPhot';
                        nonStimPhot = nonStimPhot';

                        stimTime = stimTime';
                        nonStimTime = nonStimTime';

                        %trimming
                        stimVol(isnan(stimVol)) = [];
                        nonStimVol(isnan(nonStimVol)) = [];

                        stimPhot(isnan(stimPhot)) = [];
                        nonStimPhot(isnan(nonStimPhot)) = [];

                        stimTime(isnan(stimTime)) = [];
                        nonStimTime(isnan(nonStimTime)) = [];

                        %reducing last nonStimVol because empirically runs into next blank
                        if i+1 == 7603 %LRInds(end)
                            nonStimVol = nonStimVol(1, 1:ceil( size( nonStimVol, 2 )/2 ) );
                            nonStimPhot = nonStimPhot(1, 1:ceil( size(nonStimPhot, 2 )/2 ) );
                            nonStimTime = nonStimTime(1, 1:ceil( size(nonStimTime, 2 )/2 ) );
                        end

                        % stimCheck(LRCount, 1:size(stimVol, 2)) = stimVol;

                        %images (64x64 avg)
                        LRImages(:, :, [1: size(stimVol, 2)], LRCount) = stimVolImages(:, :, 1: size(stimVol, 2));
                        LRImages(:, :, [size(stimVol, 2) + 1 : size(stimVol, 2) + size(nonStimVol, 2)], LRCount) = nonStimVolImages(:, :, 1:size(nonStimVol, 2));

                        %vols
                        LRVols(LRCount, [1: size(stimVol, 2)], 1) = stimVol;
                        LRVols(LRCount, [size(stimVol, 2) + 1 : size(stimVol, 2) + size(nonStimVol, 2)], 1) = nonStimVol;
                        %photodiode
                        LRVols(LRCount, [1: size(stimPhot, 2)], 2) = stimPhot;
                        LRVols(LRCount, [size(stimPhot, 2) + 1 : size(stimPhot, 2) + size(nonStimPhot, 2)], 2) = nonStimPhot;
                        %time
                        LRVols(LRCount, 1: size(stimTime, 2), 3) = stimTime;
                        LRVols(LRCount, size(stimTime, 2) + 1 : size(stimTime, 2) + size(nonStimTime, 2), 3) = nonStimTime;
                        LRCount = LRCount + 1;
                    end

                    %currently not saved and unsure of where to put
                    photoCorrected = squeeze(LRVols(:, :, 2));
                    for rep = 1:size(photoCorrected, 1)
                        photoCorrected(rep, :) = photoCorrected(rep, :) - photoCorrected(rep, 1);
                    end

                    %trimming
                    LRSizes = nan(1, size(LRVols, 1));
                    for rep = 1:size(LRVols, 1)
                        thisLRSize = find(~isnan(LRVols(rep, :, 1)), 1, 'last');
                        LRSizes(rep) = thisLRSize;
                    end

                    LRSizes = mode(LRSizes);
                    for rep = 1:size(LRVols, 1)
                        if ~isnan(LRVols(rep, LRSizes + 1, 1))
                            LRVols(rep, LRSizes + 1: end, :) = nan;
                        end
                    end

                    %making time consistent for each rep
                    segmentedLRTimes = nan(size(LRVols, [1,2]));
                    for rep = 1:size(LRVols, 1)
                        segmentedLRTimes(rep, :) = LRVols(rep, :, 3) - LRVols(rep,1,3);
                    end

                    %code for avg image
                    % avgLRImage = nanmean(LRImages, 4);
                    % figure
                    % for v = 1:LRSizes
                    %     image(avgLRImage(:,:,v)*200)
                    %     title([num2str(v)])
                    %     drawnow
                    %     pause(0.1)
                    % end

                    batStruct.blocks(structCount).LRVols = LRVols;
                    batStruct.blocks(structCount).segmentedLRTimes = segmentedLRTimes;
                    batStruct.blocks(structCount).LRFrames = nanmean(LRImages, 4);
                    %maybe save image here

                    %% Lum traces
                    disp('lum vols being collected')
                    % lumVols = nan(size(unique(lumInds), 2)/2, 200, 3, 7); %7 hardcoded but pretty safe bet given how stimulus script coded
                    lumVols = nan(size(unique(lumInds), 2)/(2*7), size(allIVols, 3), 3, 7); %7 hardcoded but pretty safe bet given how stimulus script coded
                    lumCount = ones(1, 7);
                    tic

                    for i = lumInds(1) + 1:2:lumInds(end) + 1 %7605:2:7884 %lumInds hardcoded [83]
                        lumIndex = matSave.totalBrightSequence(i-1) - 1; % -1 (within loop) because i out by ~1, -1 (end) because 7 luminosities but 8 options in totalBrightSequence

                        % for i = 7245:2:7524 % [87]
                        % stimVol = allIVols(yDim, xDim, :, 1, i);
                        stimVol = allIVols(:, :, :, i);
                        %disp([num2str(LRCount),' - ',num2str(stimVol(1,1,1))])
                        stimPhot = allIExtras(1, 2, :, i);
                        stimTime = squeeze(allIExtras(1, 1, :, i));
                        % stimVol(isnan(stimVol(:,))) = [];
                        % nonStimVol = allIVols(yDim, xDim, :, 1, i+1);
                        nonStimVol = allIVols(:, :, :, i+1);
                        nonStimPhot = allIExtras(1, 2, :, i+1);
                        nonStimTime = squeeze(allIExtras(1, 1, :, i+1));
                        % nonstimVol(isnan(nonstimVol)) = [];

                        %averaging to get trace
                        if doSmooth
                            stimVol = squeeze(allIExtras(1,3,:,i));
                            nonStimVol = squeeze(allIExtras(1,3,:,i+1));
                        else
                            stimVol = squeeze(nanmean(stimVol, [1,2]));
                            nonStimVol = squeeze(nanmean(nonStimVol, [1,2]));
                        end

                        stimPhot = squeeze(nanmean(stimPhot, [1,2]));
                        nonStimPhot = squeeze(nanmean(nonStimPhot, [1,2]));

                        %correcting after squeeze
                        stimVol = stimVol';
                        nonStimVol = nonStimVol';

                        stimPhot = stimPhot';
                        nonStimPhot = nonStimPhot';

                        stimTime = stimTime';
                        nonStimTime = nonStimTime';

                        %trimming
                        stimVol(isnan(stimVol)) = [];
                        nonStimVol(isnan(nonStimVol)) = [];

                        stimPhot(isnan(stimPhot)) = [];
                        nonStimPhot(isnan(nonStimPhot)) = [];

                        stimTime(isnan(stimTime)) = [];
                        nonStimTime(isnan(nonStimTime)) = [];

                        %reducing last nonStimVol because empirically runs into next blank
                        if i+1 == lumInds(end) + 1
                            nonStimVol = nonStimVol(1, 1:ceil( size( nonStimVol, 2 )/2 ) );
                            nonStimPhot = nonStimPhot(1, 1:ceil( size(nonStimPhot, 2 )/2 ) );
                            nonStimTime = nonStimTime(1, 1:ceil( size(nonStimTime, 2 )/2 ) );
                        end

                        % stimCheck(LRCount, 1:size(stimVol, 2)) = stimVol;

                        %vols
                        lumVols(lumCount(lumIndex), 1: size(stimVol, 2), 1, lumIndex) = stimVol;
                        lumVols(lumCount(lumIndex), size(stimVol, 2) + 1 : size(stimVol, 2) + size(nonStimVol, 2), 1, lumIndex) = nonStimVol;
                        %photodiode
                        lumVols(lumCount(lumIndex), 1: size(stimPhot, 2), 2, lumIndex) = stimPhot;
                        lumVols(lumCount(lumIndex), size(stimPhot, 2) + 1 : size(stimPhot, 2) + size(nonStimPhot, 2), 2, lumIndex) = nonStimPhot;
                        %time
                        lumVols(lumCount(lumIndex), 1: size(stimTime, 2), 3, lumIndex) = stimTime;
                        lumVols(lumCount(lumIndex), size(stimTime, 2) + 1 : size(stimTime, 2) + size(nonStimTime, 2), 3, lumIndex) = nonStimTime;
                        lumCount(lumIndex) = lumCount(lumIndex) + 1;
                    end
                    toc

                    %
                    lumPhotoCorrected = squeeze(lumVols(:, :, 2));
                    for rep = 1:size(lumPhotoCorrected, 1)
                        lumPhotoCorrected(rep, :) = lumPhotoCorrected(rep, :) - lumPhotoCorrected(rep, 1);
                    end

                    %trimming (so that cross repetition analyses can be done)
                    lumSizes = nan(size(lumVols, [4, 1]));
                    for lum = 1:size(lumVols, 4)
                        for rep = 1:size(lumVols, 1)
                            thisLumSize = find(~isnan(lumVols(rep, :, 1, lum)), 1, 'last');
                            lumSizes(lum, rep) = thisLumSize;
                        end
                    end

                    lumSizes = mode(lumSizes, 2);
                    for lum = 1:size(lumVols, 4)
                        for rep = 1:size(lumVols, 1)
                            if ~isnan(lumVols(rep, lumSizes(lum) + 1, 1, lum) )
                                lumVols(rep, lumSizes(lum)+1:end, :, lum) = nan;
                            end
                        end
                    end


                    %making time consistent for each luminosity type
                    segmentedLumTimes = nan(size(lumVols, [1,2,4]));
                    for lum = 1:size(lumVols, 4)
                        for rep = 1:size(lumVols, 1)
                            segmentedLumTimes(rep, :, lum) = lumVols(rep, :, 3, lum) - lumVols(rep, 1, 3, lum);
                        end
                    end

                    batStruct.blocks(structCount).lumVols = lumVols;
                    batStruct.blocks(structCount).segmentedLumTimes = segmentedLumTimes;

                    %% Freq
                    disp('freq vols being collected')
                    sepFreqInds = []; %an array containing indices where freqSequence changes
                    for i = freqInds(1):freqInds(end) + 1 %7885:8925 %hardcoded freqInds
                        if freqSequence(i-1) ~= freqSequence(i)
                            % if freqSequence(i) == 0.5
                            freqInfo = [i;freqSequence(i)];
                            sepFreqInds = [sepFreqInds, freqInfo];
                        end
                    end

                    % sepFreqInds(2,:) = sepFreqInds(2,:) + 1; %correction for freqSequence because i is empirically out by ~1 i - faulty should be applied to row 1, but wont be because freqSequence is being used without tampering

                    %where the magic happens
                    freqVols = nan(matSave.numFreq*2, size(allIVols, 3), 3, size(matSave.batFrequencies, 2));
                    % freqVols = nan(matSave.numFreq*size(matSave.batFrequencies, 2)*2, 200, 3, size(matSave.batFrequencies, 2));
                    % freqVols = nan(matSave.numFreq*size(matSave.batFrequencies, 2)*2, 200, 3); % , size(matSave.batFrequencies, 2)...
                    freqCount = ones(1, size(freqVols, 4));
                    freqOrder = matSave.freqOrder;
                    freqOrdering = [];
                    stimVol = [];
                    nonStimVol = []; %prob add more for phot and time

                    tic
                    for i = sepFreqInds(1,:)%7885:8925 %hardcoded freqInds - may want this to be first inds (presentFereqInds
                        if ismember(freqSequence(i), matSave.batFrequencies) %ismember(freqSequence(i-1), matSave.batFrequencies)
                            freqIndice = find(matSave.batFrequencies == freqSequence(i));
                            stimVol = []; %may have to pre-allocate if slow
                            stimPhot = [];
                            stimTime = [];

                            if freqSequence(i) == 1
                                freqI = i+1;
                                if doSmooth
                                    currentVol = squeeze(allIExtras(1,3,:,freqI));
                                else
                                    currentVol = squeeze(nanmean(allIVols(:, :, :, freqI), [1,2])); %should probably make yDim, xDim here
                                end
                                currentVol(isnan(currentVol)) = []; %cleaning up current vol

                                currentPhot = squeeze(allIExtras(1,2,:,freqI));
                                currentPhot(isnan(currentPhot)) = [];

                                currentTime = squeeze(allIExtras(1,1,:,freqI));
                                currentTime(isnan(currentTime)) = [];

                                stimVol = [stimVol, currentVol'];
                                stimPhot = [stimPhot, currentPhot'];
                                stimTime = [stimTime, currentTime'];

                                fullI = 0;
                            end

                            %number of i's to colllect is freqsequence(i-1)
                            for freqI = i+1: i + freqSequence(i) - 1
                                % if freqI == i
                                %     %soemthing to remove blank at start
                                % elseif freqI == i + freqSequence(i) - 1
                                %     % something to remove blank at end
                                %     fullI = 0;
                                % else
                                if doSmooth
                                    currentVol = squeeze(allIExtras(1,3,:,freqI));
                                else
                                    currentVol = squeeze(nanmean(allIVols(:, :, :, freqI), [1,2])); %should probably make yDim, xDim here
                                end
                                currentVol(isnan(currentVol)) = []; %cleaning up current vol

                                currentPhot = squeeze(allIExtras(1,2,:,freqI));
                                currentPhot(isnan(currentPhot)) = [];

                                currentTime = squeeze(allIExtras(1,1,:,freqI));
                                currentTime(isnan(currentTime)) = [];

                                stimVol = [stimVol, currentVol'];
                                stimPhot = [stimPhot, currentPhot'];
                                stimTime = [stimTime, currentTime'];

                                fullI = 0;
                                % end
                                % stimVol = [stimVol, currentVol'];
                                % stimPhot = [stimPhot, currentPhot'];

                            end
                        else % some form of blank
                            nonStimVol = [];
                            nonStimPhot = [];
                            nonStimTime = [];
                            for imagingI = [i, i+1]
                                if freqSequence(i-1) ~= 1 || freqSequence(i-1) == 1 && imagingI ~= i %empirical adjustment for 1Hz
                                    if doSmooth
                                        currentNonStimVol = squeeze(allIExtras(1,3,:,imagingI));
                                    else
                                        currentNonStimVol = squeeze(nanmean(allIVols(:,:,:,imagingI), [1,2]));
                                    end
                                    currentNonStimVol(isnan(currentNonStimVol)) = [];

                                    currentNonStimPhot = squeeze(allIExtras(1,2,:,imagingI));
                                    currentNonStimPhot(isnan(currentNonStimPhot)) = [];

                                    currentNonStimTime = squeeze(allIExtras(1,1,:,imagingI));
                                    currentNonStimTime(isnan(currentNonStimTime)) = [];

                                    nonStimVol = [nonStimVol; currentNonStimVol];
                                    nonStimPhot = [nonStimPhot; currentNonStimPhot];
                                    nonStimTime = [nonStimTime; currentNonStimTime];
                                end

                            end

                            fullI = 1;
                        end

                        if fullI
                            freqVols(freqCount(freqIndice), 1:size(stimVol, 2), 1, freqIndice) = stimVol;
                            freqVols(freqCount(freqIndice), size(stimVol, 2) + 1: size(stimVol, 2) + size(nonStimVol, 1), 1, freqIndice) = nonStimVol';

                            freqVols(freqCount(freqIndice), 1:size(stimPhot, 2), 2, freqIndice) = stimPhot;
                            freqVols(freqCount(freqIndice), size(stimPhot, 2) + 1 : size(stimPhot, 2) + size(nonStimPhot, 1), 2, freqIndice) = nonStimPhot';

                            freqVols(freqCount(freqIndice), 1: size(stimTime, 2), 3, freqIndice) = stimTime;
                            freqVols(freqCount(freqIndice), size(stimTime, 2) + 1 : size(stimTime, 2) + size(nonStimTime, 1), 3, freqIndice) = nonStimTime';

                            freqCount(freqIndice) = freqCount(freqIndice) + 1;
                        end
                    end
                    toc


                    % tic
                    % for i = sepFreqInds(1,:)%7885:8925 %hardcoded freqInds - may want this to be first inds (presentFereqInds
                    %     if ismember(freqSequence(i), matSave.batFrequencies) %ismember(freqSequence(i-1), matSave.batFrequencies)
                    %         freqIndice = find(matSave.batFrequencies == freqSequence(i));
                    %         stimVol = []; %may have to pre-allocate if slow
                    %         stimPhot = [];
                    %         stimTime = [];
                    %
                    %         %number of i's to colllect is freqsequence(i-1)
                    %         for freqI = i+1: i + freqSequence(i) - 1
                    %             % if freqI == i
                    %             %     %soemthing to remove blank at start
                    %             % elseif freqI == i + freqSequence(i) - 1
                    %             %     % something to remove blank at end
                    %             %     fullI = 0;
                    %             % else
                    %             currentVol = squeeze(nanmean(allIVols(:, :, :, freqI), [1,2])); %should probably make yDim, xDim here
                    %             currentVol(isnan(currentVol)) = []; %cleaning up current vol
                    %
                    %             currentPhot = squeeze(allIExtras(1,2,:,freqI));
                    %             currentPhot(isnan(currentPhot)) = [];
                    %
                    %             currentTime = squeeze(allIExtras(1,1,:,freqI));
                    %             currentTime(isnan(currentTime)) = [];
                    %
                    %             stimVol = [stimVol, currentVol'];
                    %             stimPhot = [stimPhot, currentPhot'];
                    %             stimTime = [stimTime, currentTime'];
                    %
                    %             fullI = 0;
                    %             % end
                    %             % stimVol = [stimVol, currentVol'];
                    %             % stimPhot = [stimPhot, currentPhot'];
                    %
                    %         end
                    %     else % some form of blank
                    %         nonStimVol = [];
                    %         nonStimPhot = [];
                    %         nonStimTime = [];
                    %         for imagingI = [i, i+1]
                    %             currentNonStimVol = squeeze(nanmean(allIVols(:,:,:,imagingI), [1,2]));
                    %             currentNonStimVol(isnan(currentNonStimVol)) = [];
                    %
                    %             currentNonStimPhot = squeeze(allIExtras(1,2,:,imagingI));
                    %             currentNonStimPhot(isnan(currentNonStimPhot)) = [];
                    %
                    %             currentNonStimTime = squeeze(allIExtras(1,1,:,imagingI));
                    %             currentNonStimTime(isnan(currentNonStimTime)) = [];
                    %
                    %             nonStimVol = [nonStimVol; currentNonStimVol];
                    %             nonStimPhot = [nonStimPhot; currentNonStimPhot];
                    %             nonStimTime = [nonStimTime, currentNonStimTime];
                    %
                    %         end
                    %         % nonStimVol = squeeze(nanmean(allIVols(:,:,:,i), [1,2]));
                    %         % nonStimVol(isnan(nonStimVol)) = [];
                    %         %
                    %         % nonStimPhot = squeeze(allIExtras(1,2,:,i));
                    %         % nonStimPhot(isnan(nonStimPhot)) = [];
                    %         %
                    %         % nonStimTime = squeeze(allIExtras(1,1,:,i));
                    %         % nonStimTime(isnan(nonStimTime)) = [];
                    %
                    %         fullI = 1;
                    %     end
                    %
                    %     if fullI
                    %         freqVols(freqCount(freqIndice), 1:size(stimVol, 2), 1, freqIndice) = stimVol;
                    %         freqVols(freqCount(freqIndice), size(stimVol, 2) + 1: size(stimVol, 2) + size(nonStimVol, 1), 1, freqIndice) = nonStimVol';
                    %
                    %         freqVols(freqCount(freqIndice), 1:size(stimPhot, 2), 2, freqIndice) = stimPhot;
                    %         freqVols(freqCount(freqIndice), size(stimPhot, 2) + 1 : size(stimPhot, 2) + size(nonStimPhot, 1), 2, freqIndice) = nonStimPhot';
                    %
                    %         freqVols(freqCount(freqIndice), 1: size(stimTime, 2), 3, freqIndice) = stimTime;
                    %         freqVols(freqCount(freqIndice), size(stimTime, 2) + 1 : size(stimTime, 2) + size(nonStimTime, 1), 3, freqIndice) = nonStimTime';
                    %
                    %         freqCount(freqIndice) = freqCount(freqIndice) + 1;
                    %     end
                    % end
                    % toc

                    % %trimming
                    freqSizes = nan(size(freqVols, [4, 1]));
                    for freq = 1:size(freqVols, 4)
                        for rep = 1:size(freqVols, 1)
                            thisFreqSize = find(~isnan(freqVols(rep, :, 1, freq)), 1, 'last');
                            freqSizes(freq, rep) = thisFreqSize;
                        end
                    end

                    freqSizes = mode(freqSizes, 2);
                    for freq = 1:size(freqVols, 4)
                        for rep = 1:size(freqVols, 1)
                            if ~isnan(freqVols(rep, freqSizes(freq) + 1, 1, freq) )
                                freqVols(rep, freqSizes(freq)+1:end, :, freq) = nan;
                            end
                        end
                    end


                    %making time consistent for each luminosity type
                    segmentedFreqTimes = nan(size(freqVols, [1,2,4]));
                    for freq = 1:size(freqVols, 4)
                        for rep = 1:size(freqVols, 1)
                            segmentedFreqTimes(rep, :, freq) = freqVols(rep, :, 3, freq) - freqVols(rep, 1, 3, freq);
                        end
                    end

                    batStruct.blocks(structCount).freqVols = freqVols;
                    batStruct.blocks(structCount).segmentedFreqTimes = segmentedFreqTimes;

                    disp(['Analysis complete for fly ', num2str(chosenFlies(fly)), ' block', num2str(b)])

                    %% blanks
                    blankVols = nan(size(blankInds, 2), size(allIVols, 3), 3);
                    blankCount = 1;
                    for i = blankInds
                        if i == blankInds(1) || i == blankInds(end) % opto blank and end of recording- empirically only 1 i in length
                            blankVols(blankCount, :, 1) = squeeze(nanmean(allIVols(:,:,:,i), [1,2]));
                            blankVols(blankCount, :, 2) = squeeze(allIExtras(1,2,:,i));
                            blankVols(blankCount, :, 3) = squeeze(allIExtras(1,1,:,i));
                        else
                            firstBlankVol = squeeze(nanmean(allIVols(:,:, :, i), [1,2]));
                            secondBlankVol = squeeze(nanmean(allIVols(:,:, :, i+1), [1,2]));

                            firstBlankPhot = squeeze(allIExtras(1,2,:,i));
                            secondBlankPhot = squeeze(allIExtras(1,2,:,i+1));

                            firstBlankTime = squeeze(allIExtras(1,1,:,i));
                            secondBlankTime = squeeze(allIExtras(1,1,:,i+1));

                            %trimming
                            firstBlankVol(isnan(firstBlankVol)) = [];
                            secondBlankVol(isnan(secondBlankVol)) = [];

                            firstBlankPhot(isnan(firstBlankPhot)) = [];
                            secondBlankPhot(isnan(secondBlankPhot)) = [];

                            firstBlankTime(isnan(firstBlankTime)) = [];
                            secondBlankTime(isnan(secondBlankTime)) = [];

                            %adding to blankVols
                            %vols
                            blankVols(blankCount, 1:size(firstBlankVol, 1), 1) = firstBlankVol';
                            blankVols(blankCount, size(firstBlankVol, 1) + 1: size(firstBlankVol, 1) + size(secondBlankVol, 1), 1) = secondBlankVol';
                            %photodiode
                            blankVols(blankCount, 1:size(firstBlankPhot, 1), 2) = firstBlankPhot';
                            blankVols(blankCount, size(firstBlankPhot, 1) + 1: size(firstBlankPhot, 1) + size(secondBlankPhot, 1), 2) = secondBlankPhot';
                            %time
                            blankVols(blankCount, 1:size(firstBlankTime, 1), 3) = firstBlankTime';
                            blankVols(blankCount, size(firstBlankTime, 1) + 1: size(firstBlankTime, 1) + size(secondBlankTime, 1), 3) = secondBlankTime';
                        end
                        blankCount = blankCount + 1;
                    end

                    segmentedBlankTimes = nan(size(blankVols, [1,2]));
                    for rep = 1:size(blankVols, 1)
                        segmentedBlankTimes(rep, :) = blankVols(rep, :, 3) - blankVols(rep,1,3);
                    end

                    batStruct.blocks(structCount).blankVols = blankVols;
                    batStruct.blocks(structCount).segmentedBlankTimes = segmentedBlankTimes;

                else
                    %%% AVERAGES %%%
                    if analyseAverage
                        %                     batStruct = struct; %used to store processed data

                        %             optoBlankVols = [find(iSequence == 2, 1, 'first'): find(iSequence == 2, 1, 'last')]; %hardcoded for now
                        %             optoBlankFullAvg = nanmean(FLIES(fly).BLOCKS(b).greenChannel(:, :, [optoBlankVols(1): optoBlankVols(end)]), 3);


                        %raw data collection for each condition
                        %total
                        %                 transientsRaw = []; %used later when plotting too
                        %                 for vol = 1:size(rawVols, 3)
                        %                     %                 average = squeeze(nanmean(transientsRaw([17:25], [38:44], vol), [1,2])); %hardCoded ROI
                        %                     %                 average = squeeze(nanmean(transientsRaw(:, :, vol), [1,2]));
                        %                     average = squeeze(nanmean(rawVols(yDim, xDim, vol), [1,2]));
                        %                     transientsRaw = [transientsRaw, average];
                        %                 end

                        % total average for normalisation/comparison - similar to Kai methods
                        totalFullAvg = nanmean(rawVols, 3);
                        %             totalSD =

                        %blanks
                        if isfield(matSave, 'blank')
                            %first blank is full second i
                            firstBlank = find(iSequence == 2);
                            rawBlanks = nan*ones(size(blankInds, 2), size(firstBlank, 2) + 5);
                            %                     for blank = 2:size(blankInds, 2)
                            % %                         currentBlank =
                            %                     end

                        end

                        %             clunk
                        %first blank
                        optoBlankVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, [find(iSequence == 2, 1, 'first'): find(iSequence == 2, 1, 'last')]); %hardcoded for now
                        disp(['-_- optoblank indices has been hardcoded -_-'])
                        %             optoBlankFullAvg = nanmean( FLIES(fly).BLOCKS(b).greenChannel(:, :, [optoBlankVols(1): optoBlankVols(end)]), 3);
                        optoBlankFullAvg = nanmean(optoBlankVols, 3);
                        optoBlankFullSD = nanstd(optoBlankVols, [], 3);

                        optoBlankAvgVar = squeeze(nanmean(optoBlankFullSD, [1,2]));
                        optoBlankAvgVarSEM = nanstd(optoBlankFullSD, [], [1,2])/sqrt(size(optoBlankFullSD, 1)*size(optoBlankFullSD, 1));

                        %QA for hardcoding
                        if ROI == 1
                            disp(['WARNING: code currently uses full size for SEM calc, not ROI specs'])
                        end


                        optoBlankAvgs = [];
                        for vol = 1:size(optoBlankVols, 3)
                            volAvg = squeeze(nanmean(optoBlankVols(:, :, vol), [1,2]));
                            optoBlankAvgs = [optoBlankAvgs, volAvg];
                        end

                        batStruct.opto.images.blankMean = optoBlankFullAvg;
                        batStruct.opto.images.blankSD = optoBlankFullSD;

                        %                 batStruct.FLIES.


                        figure
                        plot(FLIES(fly).BLOCKS(b).volTimes(1, [find(iSequence == 2, 1, 'first'): find(iSequence == 2, 1, 'last')]), optoBlankAvgs)
                        title('first blank transient')

                        %first 5 sec opto
                        firstOptoVol = find(iSequence > 2, 1, 'first'); %assumes opto is first condition
                        first5OptoVols = FLIES(fly).BLOCKS(b).greenChannel(:, :, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 5)]);
                        first5OptoFullAvg = nanmean(first5OptoVols, 3);
                        first5OptoFullSD = nanstd(first5OptoVols, [], 3);

                        first5OptoAvgVar = squeeze(nanmean(first5OptoFullSD, [1,2]));
                        first5OptoAvgVarSEM = nanstd(first5OptoFullSD, [], [1,2]) / sqrt(size(first5OptoFullSD, 1)*size(first5OptoFullSD, 2));

                        first5OptoAvgs = [];
                        for vol = 1:size(first5OptoVols, 3)
                            volAvg = squeeze(nanmean(first5OptoVols(:, :, vol), [1,2]));
                            first5OptoAvgs = [first5OptoAvgs, volAvg];
                        end

                        figure
                        plot(FLIES(fly).BLOCKS(b).volTimes(1, [firstOptoVol : find((volTimes(1, :) - volTimes(1, firstOptoVol)) >= 5)]), first5OptoAvgs)
                        title('first 5 seconds optomotor transient')

                        %boxplotting opto things
                        SEMOptoBlank = nanstd(optoBlankAvgs)/sqrt(size(optoBlankAvgs, 2));
                        SEMFirst5Opto = nanstd(first5OptoAvgs)/sqrt(size(first5OptoAvgs, 2));

                        %optoBox = [optoBlankAvgs', first5OptoAvgs'];
                        figure
                        bar([1,2], [nanmean(optoBlankAvgs), nanmean(first5OptoAvgs)])
                        line([1, 1], [nanmean(optoBlankAvgs) - SEMOptoBlank, nanmean(optoBlankAvgs) + SEMOptoBlank], 'LineWidth', 2, 'Color', 'k')
                        line([2,2], [nanmean(first5OptoAvgs) - SEMFirst5Opto, nanmean(first5OptoAvgs) + SEMFirst5Opto], 'LineWidth', 2, 'Color', 'k')


                        figure
                        subplot( 1, 5, [1,2] )
                        image( optoBlankFullAvg*0.1 )
                        subplot( 1, 5, [3] )
                        bar([1,2], [nanmean(optoBlankAvgs), nanmean(first5OptoAvgs)])
                        line([1, 1], [nanmean(optoBlankAvgs) - SEMOptoBlank, nanmean(optoBlankAvgs) + SEMOptoBlank], 'LineWidth', 2, 'Color', 'k')
                        line([2,2], [nanmean(first5OptoAvgs) - SEMFirst5Opto, nanmean(first5OptoAvgs) + SEMFirst5Opto], 'LineWidth', 2, 'Color', 'k')
                        %             ylim([347,353])
                        ylim([min( [nanmean(optoBlankAvgs) - SEMOptoBlank, nanmean(first5OptoAvgs) - SEMFirst5Opto] ) - 1, max( [nanmean(optoBlankAvgs) + SEMOptoBlank, nanmean(first5OptoAvgs) + SEMFirst5Opto] + 1 )])
                        subplot( 1, 5, [4,5] )
                        image( first5OptoFullAvg*0.1 )

                        %comparing SD images
                        figure
                        subplot(1,5,[1,2])
                        image( optoBlankFullSD )
                        subplot(1, 5, [3])
                        bar([1,2], [optoBlankAvgVar, first5OptoAvgVar])
                        line([1,1], [optoBlankAvgVar - optoBlankAvgVarSEM, optoBlankAvgVar + optoBlankAvgVarSEM], 'LineWidth', 2, 'Color', 'k')
                        line([2,2], [first5OptoAvgVar - first5OptoAvgVarSEM, first5OptoAvgVar + first5OptoAvgVarSEM], 'LineWidth', 2, 'Color', 'k')
                        subplot(1,5,[4,5])
                        image( first5OptoFullSD )
                        sgtitle(['average variance of blank and 5s opto'])

                        %stats for SD comparison
                        optoBlankVar = reshape(optoBlankFullSD, 4096, 1);
                        first5OptoVar = reshape(first5OptoFullSD, 4096, 1);
                        [h, pSDOpto] = ttest2(optoBlankVar, first5OptoVar);

                    end



                    %%%   TRANSIENTS   %%%
                    if analyseTransientOld
                        %Raw
                        transientsRaw = []; %used later when plotting too
                        for vol = 1:size(rawVols, 3)
                            %                 average = squeeze(nanmean(transientsRaw([17:25], [38:44], vol), [1,2])); %hardCoded ROI
                            %                 average = squeeze(nanmean(transientsRaw(:, :, vol), [1,2]));
                            % average = squeeze(nanmean(rawVols(yDim, xDim, vol), [1,2]));
                            average = squeeze(nanmean(rawVols(:, :, vol), [1,2]));
                            transientsRaw = [transientsRaw, average];
                        end


                        a = 1; %used for anything?
                        %getting frames for each respective i associated with battery conditions
                        %LR
                        transientsLR = nan(64, 64, 20, matSave.numLR/2, 2); %5-D matrix; (x, y, time/volume#, repitition, side)
                        photoLR = nan( 1, matSave.numLR/2, 2 ); % -D matrix containing photodiode data (photodiode points/TSIndex?, repetition, side)
                        for inds = LRInds(2): 2: LRInds(end)
                            %                 if freqSequence(inds) == 1
                            %                     noStimInds = [find(FLIES(fly).BLOCKS(b).iSequence == inds)];
                            %                     stimInds = [find(FLIES.BLOCKS(b).iSequence == inds - 1)]; % all volumes associated with stimulus period
                            %                     stimOffInds = [stimInds(ceil(size(stimInds, 2)/2)): max(stimInds)]; % volumes in stimulus period where stimulus is off (currently assumed as second half)

                            %surrounding i's are being accessed because freqSequence is currently out of phase
                            propPrevI = ( ( 1/freqSequence(inds-1) )/2 ) / ( ( (1/freqSequence(inds-2) ) ) +  ( ( 1/freqSequence(inds-1) )/2 ) ); %proportion of volumes to take from previous i
                            propNextI = ( ( 1/freqSequence(inds) )/2 ) / ( ( (1/freqSequence(inds+1) )/2 ) + ( (1/freqSequence(inds) )/2 ) ) ; %proportion of volumes to take from following i

                            allPrevVols = [find(FLIES.BLOCKS(b).iSequence == inds - 1)];
                            allNextVols = [find(FLIES.BLOCKS(b).iSequence == inds + 1)];

                            firstIVols = allPrevVols(end+1 - ceil(size(allPrevVols, 2)*propPrevI) : end); % end+1 because how indexing works (could've omitted if chosen floor instead of ceil, but ceil shows logic my vol selection more)
                            currentIVols = [find(FLIES.BLOCKS(b).iSequence == inds)];
                            secondIVols = allNextVols(end);

                            if periodChoice == 1
                                blankInds = noStimInds;
                            elseif periodChoice == 2
                                blankInds = [stimInds, noStimInds];
                            elseif periodChoice == 3
                                blankInds = [stimOffInds, noStimInds];
                            end


                            volumes = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);

                            repLR = 1;

                            if randomSequence(inds) == 1
                                transientsLR(:, :, [1:size(blankInds, 2)], repLR, 1) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                                repLR = repLR + 1;
                            elseif randomSequence(inds) == 2
                                transientsLR(:, :, [1:size(blankInds, 2)], repLR, 2) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                                repLR = repLR + 1;
                            else
                                disp(['|-_-| unknown side requested in side condition |-_-|'])
                            end

                            if repLR > size(transientsLR, 4)
                                repLR = 1;
                            end

                        end

                        %%old - uses freqSequence%%
                        %             for inds = LRInds
                        %                 if freqSequence(inds) == 1
                        %                     noStimInds = [find(FLIES.BLOCKS(b).iSequence == inds)];
                        %                     stimInds = [find(FLIES.BLOCKS(b).iSequence == inds - 1)]; % all volumes associated with stimulus period
                        %                     stimOffInds = [stimInds(ceil(size(stimInds, 2)/2)): max(stimInds)]; % volumes in stimulus period where stimulus is off (currently assumed as second half)
                        %
                        %                     if plottingChoice == 1
                        %                         blankInds = noStimInds;
                        %                     elseif plottingChoice == 2
                        %                         blankInds = [stimInds, noStimInds];
                        %                     elseif plottingChoice == 3
                        %                         blankInds = [stimOffInds, noStimInds];
                        %                     end
                        %
                        %
                        %                     volumes = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                        %
                        %                     repLR = 1;
                        %
                        %                     if randomSequence(inds) == 1
                        %                         transientsLR(:, :, [1:size(blankInds, 2)], repLR, 1) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                        %                         repLR = repLR + 1;
                        %                     elseif randomSequence(inds) == 2
                        %                         transientsLR(:, :, [1:size(blankInds, 2)], repLR, 2) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                        %                         repLR = repLR + 1;
                        %                     else
                        %                         disp(['|-_-| unknown side requested in side condition |-_-|'])
                        %                     end
                        %
                        %                     if repLR > size(transientsLR, 4)
                        %                         repLR = 1;
                        %                     end
                        %                 end
                        %             end

                        %             for inds = LRInds
                        %                 thisInds = [find(FLIES.BLOCKS(b).iSequence == inds)];
                        % %                 disp( [num2str(inds), ' - ', num2str(thisInds)] )
                        %                 volumes = FLIES.BLOCKS(b).greenChannel(:, :, thisInds);
                        %                 volAvg = squeeze(nanmean(volumes, [1,2]));
                        %                 transientsLR( a, [1:size(volAvg,1)] ) = volAvg;
                        %                 a = a + 1;
                        %             end

                        %Lum
                        transientsLum = nan(64, 64, 20, matSave.numLum/(7*2),  8, 2); %6-D matrix; (x, y, time/volume#, repition, brightness, side
                        repLum = ones(2, size(matSave.brightIndex, 2));
                        batLuminosities = [0: size(matSave.brightIndex, 2)];

                        for inds = lumInds
                            if freqSequence(inds) == 1
                                %                     blankInds = [find(FLIES.BLOCKS(b).iSequence == inds)];
                                noStimInds = [find(FLIES.BLOCKS(b).iSequence == inds)];
                                stimInds = [find(FLIES.BLOCKS(b).iSequence == inds - 1)]; % all volumes associated with stimulus period
                                stimOffInds = [stimInds(ceil(size(stimInds, 2)/2)): max(stimInds)]; % volumes in stimulus period where stimulus is off (currently assumed as second half)

                                if periodChoice == 1
                                    blankInds = noStimInds;
                                elseif periodChoice == 2
                                    blankInds = [stimInds, noStimInds];
                                elseif periodChoice == 3
                                    blankInds = [stimOffInds, noStimInds];
                                end

                                volumes = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);

                                lumIndex = matSave.totalBrightSequence(inds - 1);

                                if randomSequence(inds) == 1
                                    transientsLum(:, :, [1:size(blankInds, 2)], repLum(1, lumIndex), lumIndex, 1) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                                    repLum(1, lumIndex) = repLum(1, lumIndex) + 1;
                                elseif randomSequence(inds) == 2
                                    transientsLum(:, :, [1:size(blankInds, 2)], repLum(2, lumIndex), lumIndex, 2) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                                    repLum(2, lumIndex) = repLum(2, lumIndex) + 1;
                                else
                                    disp(['|-_-| unknown side requested in luminosity condition |-_-|'])
                                end
                            end
                        end


                        %             transientsFreq = nan(size(freqInds, 2) , 10 , size(matSave.batFrequencies, 2) ); %first col will have frequency tagline, 2-: will have volumes for transient
                        transientsFreq = nan(64, 64, 10, matSave.numFreq, size(matSave.batFrequencies, 2), 2); %6-D matrix; (x, y, time/volume#, repitition, frequency, side)
                        repFreq = ones(2, size(matSave.batFrequencies, 2)); %A matrix used to track repitions for each frequency

                        %Freq
                        %currently set to ignore stimulus period - could be inaccurate given timing issue
                        for inds = freqInds
                            if freqSequence(inds) == 0.5 && freqSequence(inds - 1) ~= 0.5
                                blankInds = [find(FLIES.BLOCKS(b).iSequence == inds)];
                                %                     disp( [num2str(inds), ' - ', num2str(blankInds)] )
                                volumes = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);

                                freqTypeIndex = find(matSave.batFrequencies == freqSequence(inds - 1));

                                if randomSequence(inds) == 1
                                    transientsFreq(:, :, [1:size(blankInds, 2)], repFreq(1, freqTypeIndex), freqTypeIndex, 1) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                                    repFreq(1, freqTypeIndex) = repFreq(1, freqTypeIndex) + 1;
                                elseif randomSequence(inds) == 2
                                    transientsFreq(:, :, [1:size(blankInds, 2)], repFreq(2, freqTypeIndex), freqTypeIndex, 2) = FLIES.BLOCKS(b).greenChannel(:, :, blankInds);
                                    repFreq(2, freqTypeIndex) = repFreq(2, freqTypeIndex) + 1;
                                else
                                    disp(['|-_-| weird things happening in freq condition - unknown side called for |-_-|'])
                                end

                                %                     repFreq = repFreq + 1;

                                %                     volAvg = squeeze(nanmean(volumes, [1,2]));

                                %                     transientsFreq(a, 1) = freqSequence(inds - 1);
                                %                     transientOrder
                                %                     transientsFreq( a, [1:size(volAvg, 1)] ) = volAvg;
                                %                     a = a + 1;
                            end
                        end
                    end

                    %%%% Plotting Stuff %%%%

                    if figChoice(1) == 1
                        figure
                        hold on
                        plot(FLIES.BLOCKS(b).volTimes(1, :), transientsRaw)
                        % plot(FLIES.BLOCKS(b).volTimes(1, :), highpass( transientsRaw, 0.25, 3.33 ) + 1300)
                        iValRaw = 1;
                        for vol = 1:size(iSequence, 2)
                            if iSequence(vol) ~= iValRaw
                                line([volTimes(1, vol), volTimes(1, vol)], [min(transientsRaw) - 10, max(transientsRaw) + 10], 'Color', 'k')
                                %                     text(volTimes(vol), min(rawAvgs) - 2, ['i=', num2str(iValRaw)], 'Color', 'r')
                                text(volTimes(1, vol), min(transientsRaw) - 2, ['i=', num2str(iSequence(vol))], 'Color', 'r')
                                iValRaw = iSequence(1, vol);
                            end
                        end
                        if ROI == 1
                            title(['(', num2str(dimensionsROI{1}(1)), ',', num2str(dimensionsROI{1}(2)), ')', ' x ', '(', num2str(dimensionsROI{2}(1)), ',', num2str(dimensionsROI{2}(2)), ')', ' - ', regionName])
                        else
                            title([regionName])
                        end
                    end

                    %             nugget

                    %LR
                    if figChoice(2) == 1
                        figure
                        hold on
                        for side = [1,2]
                            LRFullAvg = [];
                            for rep = 1:size(transientsLR, 4)
                                LRRepAvg = squeeze(nanmean(transientsLR(:, :, :, rep, side), [1,2]));
                                LRFullAvg = [LRFullAvg, LRRepAvg];
                            end
                            LRFullAvg = nanmean(LRFullAvg, 2);
                            plot(LRFullAvg)

                        end
                        title(['Single Bar Transients'])
                        legend([{'Left'}, {'Right'}])
                    end

                    %Lum
                    if figChoice(3) == 1
                        for side = [1,2]
                            lumLegend = [];

                            figure
                            hold on
                            for lum = 1:size(transientsLum, 5)
                                %writing legend
                                lumLevel = matSave.brightIndex(lum);
                                lumLegend = [lumLegend, {[num2str(floor(lumLevel*100)), '%']}];

                                lumFullAvg = [];
                                for rep = 1:size(transientsLum, 4)
                                    lumRepAvg = squeeze(nanmean(transientsLum(:, :, :, rep, lum, side), [1,2]));
                                    lumFullAvg = [lumFullAvg, lumRepAvg]; %collecting reps per frequency to be averaged outside loop
                                end
                                lumFullAvg = nanmean(lumFullAvg, 2);
                                plot(lumFullAvg)
                            end
                            legend(lumLegend)
                            if side == 1
                                title(['f:', num2str(thisFly), ' b:', num2str(b), ' Luminosity Transients Left'])
                            elseif side == 2
                                %                     title(['Luminosity Transients Right'])
                                title(['f:', num2str(thisFly), ' b:', num2str(b), ' Luminosity Transients Right'])
                            end
                        end
                    end

                    %Freq
                    if figChoice(4) == 1
                        %             fullFieldAvgFreq = [];
                        for side = [1,2]
                            freqLegend = [];

                            figure
                            hold on
                            for freq = 1:size(transientsFreq, 5)
                                %writing legend
                                freqType = matSave.batFrequencies(freq);
                                freqLegend = [freqLegend, {num2str(freqType)}];

                                freqFullAvg = [];
                                for rep = 1:size(transientsFreq, 4)
                                    freqRepAvg = squeeze(nanmean(transientsFreq(:, :, :, rep, freq, side), [1,2]));
                                    freqFullAvg = [freqFullAvg, freqRepAvg]; %collecting reps per frequency to be averaged outside loop
                                    %                     freqType = matSave.batFrequencies(freq);
                                    %                     freqLeg = [freqLeg, {num2str(freqType)}];
                                end
                                freqFullAvg = nanmean(freqFullAvg, 2);
                                plot(freqFullAvg)
                            end
                            legend(freqLegend)
                            if side == 1
                                title(['f:', num2str(thisFly), ' b:', num2str(b), ' Frequency Transients Left'])
                            elseif side == 2
                                title(['f:', num2str(thisFly), ' b:', num2str(b), ' Frequency Transients Right'])
                            end
                        end
                    end

                    %             timingLR = [];
                    %             %calculating time between each i for LR
                    %             for i = LRInds
                    %                 timeDiff = BTData(find(BTData(:, 5) == i, 1,'last'), 6) - BTData(find(BTData(:, 5) == i, 1, 'first'), 6);
                    %                 timingLR = [timingLR, timeDiff];
                    %             end
                    %
                    %             timingLum = [];
                    %             for i = lumInds
                    %                 timeDiff = BTData(find(BTData(:, 5) == i, 1,'last'), 6) - BTData(find(BTData(:, 5) == i, 1, 'first'), 6);
                    %                 timingLum = [timingLum, timeDiff];
                    %             end
                    %
                    %             timingFreq = [];
                    %             for i = freqInds
                    %                 timeDiff = BTData(find(BTData(:, 5) == i, 1,'last'), 6) - BTData(find(BTData(:, 5) == i, 1, 'first'), 6);
                    %                 timingFreq = [timingFreq, timeDiff];
                    %             end
                    %
                    %             %sorting transient frames into blanks vs stim on
                    % %             for i = 1: size(transientsLR, 1)
                    % %                 if ghh
                    % %                 end
                    % %             end

                end
            end
            structCount = structCount + 1;
            if ROI
                roiCount = roiCount + 1;
            end
        end
    end

    % woopah
    %% Cross fly Analysis
    blockNums = [batStruct.blocks.block];
    if conditionChoice(1)
    %opto

    %averaging over blocks
    % optoBlankData = []; %maybe rename to firstBlankData
    % optoData = [];
    % flyChanges = []; %locations where fly switches in batStruct - used to index specific blocks later
    % previousFly = 0;

    totalExpAvg = nan(max([batStruct.blocks.block]), size(chosenFlies, 2));
    optoData = nan(max([batStruct.blocks.block]), size(chosenFlies, 2), 10); %currently 4 for 3rd dim, (this dim will contain data type e.g. opto, blank, raw + normalised) - 4 dims + respective normalised versions = 8
    blockNums = [batStruct.blocks.block];

    for row = 1:size(batStruct.blocks, 2)
        totalExpAvg(blockNums(row), find(chosenFlies == batStruct.blocks(row).fly)) = squeeze(nanmean(batStruct.blocks(row).wholeExpAvg, [1,2]));

        optoData(blockNums(row), find(chosenFlies == batStruct.blocks(row).fly), 1) = batStruct.blocks(row).optoAvg;
        optoData(blockNums(row), find(chosenFlies == batStruct.blocks(row).fly), 2) = batStruct.blocks(row).blankAvg;
        optoData(blockNums(row), find(chosenFlies == batStruct.blocks(row).fly), 5) = batStruct.blocks(row).optoVar;
        optoData(blockNums(row), find(chosenFlies == batStruct.blocks(row).fly), 6) = batStruct.blocks(row).blankVar;
        optoData(blockNums(row), find(chosenFlies == batStruct.blocks(row).fly), 9) = batStruct.blocks(row).fullOptoAvg;
        optoData(blockNums(row), find(chosenFlies == batStruct.blocks(row).fly), 10) = batStruct.blocks(row).fullOptoVar;
        %may have to reconfig to have structure of [avgs, var] but probably best to just turn it into a structure instead

    end


    % for i = 1:size(batStruct.blocks, 2)
    % %     %index help
    %     if batStruct.blocks(i).fly ~= previousFly
    %         flyChanges = [flyChanges, i];
    %         previousFly = batStruct.blocks(i).fly;
    %     end
    %     %data retrievement
    %     structRow = batStruct.blocks(i);
    %     optoBlankData = [optoBlankData, [structRow.blankAvg ; structRow.blankVar] ];
    %     optoData = [optoData, [structRow.optoAvg ; structRow.optoVar] ];
    % end

    % Normalising
    for dataType = [1,2 , 5,6] %unmodified data always stored in these locations
        for block = 1:size(optoData, 1)
            correctionFactor = nanmean( optoData(block, :, dataType) );
            optoData(block, :, dataType + 2) =  optoData(block, :, dataType) - correctionFactor; %stored in dataType +2 because currently alway 2 conditions (blank and batCondition), therefore next logical location for corrected is 2 along
        end
    end
    %CURRENTLY WRONG CALC ^^^

    % optoBlankCorrected = [];
    % optoCorrected = [];
    %
    % for rep = 1:size(optoData, 2)
    %     correctFactorAvg = nanmean([optoBlankData(1, rep), optoData(1, rep)]);
    %     correctFactorVar = nanmean([optoBlankData(2, rep), optoData(2, rep)]);
    %     optoBlankCorrected = [optoBlankCorrected, [optoBlankData(1, rep) - correctFactorAvg ; optoBlankData(2, rep) - correctFactorVar]];
    %     optoCorrected = [optoCorrected, [optoData(1, rep) - correctFactorAvg ; optoData(2, rep) - correctFactorVar]];
    % end

    %Stats
    %whole Average
    totalSEM = nan(size(totalExpAvg, 1), 1);
    for block = 1:size(totalExpAvg, 1)
        totalSEM(block, 1) = nanstd(totalExpAvg(block, :), [], 2) / sqrt( size( find( ~isnan( totalExpAvg(block, :) ) ), 2 ) );
    end

    %opto Data
    optoSEM = nan(size(optoData, 1), size(optoData, 3)); %storage of SEM for opto condition, rows=block, cols=dataCondition
    for dataType = 1:size(optoData, 3)
        for block = 1:size(optoData, 1)
            optoSEM(block, dataType) = nanstd(optoData(block, :, dataType), [], 2)/ sqrt(size( find( ~isnan( optoData(block, :, dataType) ) ), 2 ) );
        end
    end


    % %which ttest to use?
    % %ttest = One-sample and paired-sample t-test.
    % %ttest2 = Two-sample t-test with pooled or unpooled variance estimate. <- unpaired t-test
    % [optoH, optoP] = ttest2(optoBlankCorrected(1, :), optoCorrected(1, :));
    % % [optoVarH, optoVarP] = [];

    %%% PLOTTING STUFF %%%

    % %code for creating subplot of 4 flies
    % figure
    % for i = 1:size(blockNums, 2)
    %     if batStruct.blocks(i).block == 1
    %         subplot(2, 2, [count])
    %         imagesc(batStruct.blocks(i).wholeExpAvg)
    %         title(['fly ', num2str(batStruct.blocks(i).fly)])
    %         count = count + 1;
    %     end
    % end

    %opto
    %quick plot of opto blank against photo using allIVols(isortedVols;iVols) for b1
    %put in a memory-saving variable option which determines whether allIVols is saved
    % % for currentFly = find(blockNums == 1)
    % %     iVols = batStruct.blocks(currentFly).iSortedVols(:,:,:,:,2); %2 hardcoded blank Ind for opto
    % %     figure
    % %     hold on
    % %     plot(squeeze(iVols(1,1,:,2) ), squeeze(nanmean(iVols(:, :, :, 1), [1,2])))
    % %
    % %     % photoRange = range(squeeze(iVols(1,2,:,2)));
    % %     % volRange = range(squeeze(nanmean(iVols(:, :, :, 1), [1,2])));
    % %
    % %     photoScale = (range(squeeze(nanmean(iVols(:, :, :, 1), [1,2]))) / range(squeeze(iVols(1,2,:,2))) ) /2 ; % half of calcium range / photodiode range
    % %
    % %     plot(squeeze(iVols(1,1,:,2) ), squeeze(iVols(1,2,:,2) )*photoScale , 'Color', 'b')
    % %     xlabel('time (s)')
    % %     ylabel('mean flurorescence')
    % %     set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    % %
    % %     title(['fly ', num2str(batStruct.blocks(batFly).fly), ' block 1 calcium trace opto blank (from allIVols)'])
    % % end

    %getting avg trace
    % % % optoTraces = nan(size(chosenFlies, 2), 100, 2, 2); %fly, vols, time or data, blank or opto
    % % % for block = blockNums == 1
    % % %     optoTraces(find(chosenFlies == batStruct.blocks(block).fly), :, 1, 1) = squeeze(nanmean(batStruct.blocks(block).blankRaw, [1,2]));
    % % %     optoTraces(find(chosenFlies == batStruct.blocks(block).fly), :, 2, 1) = batStruct.blocks(block).blankRawTimes;
    % % %     optoTraces(find(chosenFlies == batStruct.blocks(block).fly), :, 1, 2) = squeeze(nanmean(batStruct.blocks(block).optoRaw, [1,2]));
    % % %     optoTraces(find(chosenFlies == batStruct.blocks(block).fly), :, 2, 2) = batStruct.blocks(block).optoRawTimes;
    % % % end

    % if conditionChoice(1)
        %traces
        %generating lots of figures
        for i = find( blockNums == 1) %1:size(blockNums, 2)
            figure
            hold on
            % subplot(1, 2, [2])
            plot([batStruct.blocks(i).optoRawTimes], [squeeze(nanmean(batStruct.blocks(i).optoRaw, [1,2]))])
            plot([batStruct.blocks(i).blankRawTimes], [squeeze(nanmean(batStruct.blocks(i).blankRaw, [1,2]))], 'Color', [1, .5, 0])
            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
            ylabel('mean Flurorescence')
            xlabel('time(s)')
            title(['fly ',num2str(batStruct.blocks(i).fly), ' - b',num2str(batStruct.blocks(i).block), ' first 5s opto']) %temp title

            % sgtitle(['fly ',num2str(batStruct.blocks(i).fly), ' - b',num2str(batStruct.blocks(i).block)]) %temp title
        end

        %little tidbit to add empty space to blank
        % thisFly = batStruct.blocks(7); %fly choice
        % plot([thisFly.blankRawTimes(end), thisFly.optoRawTimes(1)], [squeeze(nanmean(thisFly.blankRaw(:,:,end), [1,2])), squeeze( nanmean( thisFly.optoRaw(:,:, 1), [1,2]))], 'Color', [1, .5, 0])

        %full opto traces
        for i = 1:size(blockNums, 2)
            figure
            hold on
            plot([batStruct.blocks(i).fullOptoRawTimes], [squeeze(nanmean(batStruct.blocks(i).fullOptoRaw, [1,2]))])
            plot([batStruct.blocks(i).blankRawTimes], [squeeze(nanmean(batStruct.blocks(i).blankRaw, [1,2]))], 'Color', [1, .5, 0])
            set(gca, 'LineWidth', 1.25, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
            ylabel('mean Flurorescence')
            xlabel('time(s)')
            title(['fly ',num2str(batStruct.blocks(i).fly), ' - b',num2str(batStruct.blocks(i).block), ' full opto']) %temp title

            yl = ylim; %assuming matlab knows what's best for me
            for x = linspace(batStruct.blocks(i).fullOptoRawTimes(1), batStruct.blocks(i).fullOptoRawTimes(end), 20) %assumes numOpto of 10
                line([x, x], [yl(1), yl(2)], 'Color', 'k')
            end
        end

        % subplot(max(blockNums), size(chosenFlies, 2), [blockNum + flyNum])

        %uncorrected
        %5s
        figure
        hold on
        % can either
        optoBlankFig = bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))], 'FaceColor', 'flat', 'LineWidth', 2);
        optoBlankFig.CData(1, :) = [1, 0.5, 0];
        % bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))])
        errorbar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))], [optoSEM(1,2), optoSEM(1,1)], 'LineStyle', 'none', 'LineWidth', 2.5, 'Color', 'k') % I like width of 1.25

        for i = 1:size(optoData, 2)
            plot([1,2], [optoData(1, i, 2), optoData(1, i, 1)], '-o', 'Color', 'r', 'LineWidth', 2)
        end

        xticks([1,2])
        xticklabels({'Blank period', 'First 5s Optomotor'})
        xtickangle(30)
        xlabel([newline, 'Condition'])
        ylabel(['Mean Fluorescence'])
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

        title(['Mean fluorescence blank vs first'])

        %full
        figure
        hold on
        % can either
        fullOptoFig = bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 9))], 'FaceColor', 'flat', 'LineWidth', 2);
        fullOptoFig.CData(1, :) = [1, 0.5, 0];
        % bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))])
        errorbar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 9))], [optoSEM(1,2), optoSEM(1,9)], 'LineStyle', 'none', 'LineWidth', 2.5, 'Color', 'k') % I like width of 1.25

        for i = 1:size(optoData, 2)
            plot([1,2], [optoData(1, i, 2), optoData(1, i, 9)], '-o', 'Color', 'r', 'LineWidth', 2)
        end

        xticks([1,2])
        xticklabels({'Blank period', 'Full Optomotor'})
        xtickangle(30)
        xlabel([newline, 'Condition'])
        ylabel(['Mean Fluorescence'])
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

        title(['Mean fluorescence blank vs full'])

        %code to plot against total average instead of blank
        figure
        optoAvgFig = bar([1,2], [nanmean(totalExpAvg(1, :)), nanmean(optoData(1, :, 1))], 'FaceColor', 'flat', 'LineWidth', 2);
        optoAvgFig.CData(1, :) = [1, 0.5, 0];
        xticklabels({'totalExpAvg', 'first 5s opto'})
        line([1,1], [nanmean(totalExpAvg(1,:)) - totalSEM(1), nanmean(totalExpAvg(1, :)) + totalSEM(1)], 'LineWidth', 2, 'Color', 'k')
        line([2,2], [nanmean(optoData(1, :, 1)) - optoSEM(1,1), nanmean(optoData(1, :, 1)) + optoSEM(1,1)], 'LineWidth', 2, 'Color', 'k')
        title('first 5s opto vs avg')
        set(gca, 'box', 'off', 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

        % code to generate very rough imaging traces
        % for i = 1:size(batStruct.blocks,2)
        %     figure
        %     hold on
        %     plot([batStruct.blocks(i).optoRawTimes], [squeeze(nanmean(batStruct.blocks(i).optoRaw, [1,2]))])
        %     plot([batStruct.blocks(i).blankRawTimes], [squeeze(nanmean(batStruct.blocks(i).blankRaw, [1,2]))], 'Color', [1, .5, 0])
        %     title(['fly ',num2str(batStruct.blocks(i).fly), ' - b',num2str(batStruct.blocks(i).block)])
        % end

        %opto averages (blank vs first5) for first 3 blocks
        figure
        hold on
        subplot(1, 3, [1])
        optoFig = bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))], 'FaceColor', 'flat', 'LineWidth', 2);
        optoFig.CData(1, :) = [1, 0.5, 0];
        line([1,1], [nanmean(optoData(1, :, 2)) - optoSEM(1,2), nanmean(optoData(1, :, 2)) + optoSEM(1,2)], 'LineWidth', 2, 'Color', 'k')
        line([2,2], [nanmean(optoData(1, :, 1)) - optoSEM(1,1), nanmean(optoData(1, :, 1)) + optoSEM(1,1)], 'LineWidth', 2, 'Color', 'k')
        xticks([1,2])
        xticklabels({'Blank period', 'First 5s Optomotor'})
        title('Block 1')
        set(gca, 'box', 'off', 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

        % ylim([min( optoData(1, :, [1,2]) ) - max(optoSEM(1, [1,2]) ) , max( optoData(1, :, [1,2]) + max(optoSEM(1, [1,2]) ) )])

        subplot(1, 3, [2])
        optoFig = bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))], 'FaceColor', 'flat', 'LineWidth', 2);
        optoFig.CData(1, :) = [1, 0.5, 0];
        line([1,1], [nanmean(optoData(1, :, 2)) - optoSEM(1,2), nanmean(optoData(1, :, 2)) + optoSEM(1,2)], 'LineWidth', 2, 'Color', 'k')
        line([2,2], [nanmean(optoData(1, :, 1)) - optoSEM(1,1), nanmean(optoData(1, :, 1)) + optoSEM(1,1)], 'LineWidth', 2, 'Color', 'k')
        xticks([1,2])
        xticklabels({'Blank period', 'First 5s Optomotor'})
        title('Block 2')
        set(gca, 'box', 'off', 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

        subplot(1, 3, [3])
        optoFig = bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))], 'FaceColor', 'flat', 'LineWidth', 2);
        optoFig.CData(1, :) = [1, 0.5, 0];
        line([1,1], [nanmean(optoData(1, :, 2)) - optoSEM(1,2), nanmean(optoData(1, :, 2)) + optoSEM(1,2)], 'LineWidth', 2, 'Color', 'k')
        line([2,2], [nanmean(optoData(1, :, 1)) - optoSEM(1,1), nanmean(optoData(1, :, 1)) + optoSEM(1,1)], 'LineWidth', 2, 'Color', 'k')
        xticks([1,2])
        xticklabels({'Blank period', 'First 5s Optomotor'})
        title('Block 3')
        set(gca, 'box', 'off', 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

        %^^This will probably need to be made into a loop because number of blocks is variable^^

        % subplot(1, 3, [2])
        % bar

        % ylim([min( optoData(1, :, [1,2]) ) - max(optoSEM(1, [1,2]) ) , max( optoData(1, :, [1,2]) + max(optoSEM(1, [1,2]) ) )])

        %variance
        figure
        hold on
        % can either
        optoVarFig = bar([1,2], [nanmean(optoData(1, :, 6)), nanmean(optoData(1, :, 5))], 'FaceColor', 'flat', 'LineWidth', 2);
        optoVarFig.CData(1, :) = [1, 0.5, 0];
        % bar([1,2], [nanmean(optoData(1, :, 2)), nanmean(optoData(1, :, 1))])
        errorbar([1,2], [nanmean(optoData(1, :, 6)), nanmean(optoData(1, :, 5))], [optoSEM(1,6), optoSEM(1,5)], 'LineStyle', 'none', 'LineWidth', 2.5, 'Color', 'k') % I like width of 1.25

        for i = 1:size(optoData, 2)
            plot([1,2], [optoData(1, i, 6), optoData(1, i, 5)], '-o', 'Color', 'r', 'LineWidth', 2)
        end

        xticks([1,2])
        xticklabels({'Blank period', 'First 5s Opto Var'})
        xtickangle(30)
        % xlabel([newline, 'Condition'])
        ylabel(['SD of Mean Fluorescence'])
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        set(gcf, 'Name', ['cross-fly var blank v first5sec - ', regionName], 'NumberTitle', 'off', 'Renderer', 'painters')

        title(['Variance blank vs first'])
    end

    %% LR

    if conditionChoice(2)
        %block 1 analysis

        %individual raw traces with mean for each fly - LR recollection?
        tic
        % allFlyLR = nan( size(batStruct.blocks, 2), size(LRVols, 2), size(LRVols, 3) );
        % allFlyLR = nan( size(chosenFlies, 2), size(LRVols, 2), size(LRVols, 3) );
        allFlyLR = nan( size(chosenFlies, 2), size(LRVols, 2), size(LRVols, 3), max(unique(blockNums)) ); %[fly, volume, datatype, block]
        % allFlyLeft = nan(size(chosenFlies, 2), size(LRVols, 2)/2, size(LRVols, 3), max(unique(blockNums)));
        allFlyLeft = nan(size(allFlyLR));
        allFlyRight = nan(size(allFlyLeft));
        flyCount = 1;
        currentFly = batStruct.blocks(1).fly;
        for batFly = 1:size(blockNums, 2)  %find(blockNums == 2)
            if batStruct.blocks(batFly).fly ~= currentFly
                flyCount = flyCount + 1;
                currentFly = batStruct.blocks(batFly).fly;
            end

            LRVols = batStruct.blocks(batFly).LRVols;
            segmentedLRTimes = batStruct.blocks(batFly).segmentedLRTimes;
            allFlyLR(flyCount, :, 1, blockNums(batFly) ) = nanmean(LRVols(:,:,1));
            allFlyLR(flyCount, :, 2, blockNums(batFly)) = nanmean(LRVols(:,:, 2));
            allFlyLR(flyCount, :, 3, blockNums(batFly)) = nanmean(segmentedLRTimes(:, :));
            %potentially add LRVols(:,:,3) here and error bars for each average trace

            sideNum = size(LRVols, 1)/2;
            %left
            allFlyLeft(flyCount, :, 1, blockNums(batFly) ) = nanmean(LRVols(1:sideNum,:,1));
            allFlyLeft(flyCount, :, 2, blockNums(batFly)) = nanmean(LRVols(1:sideNum,:, 2));
            allFlyLeft(flyCount, :, 3, blockNums(batFly)) = nanmean(segmentedLRTimes(1:sideNum, :));

            %right
            allFlyRight(flyCount, :, 1, blockNums(batFly) ) = nanmean(LRVols(sideNum+1:end,:,1));
            allFlyRight(flyCount, :, 2, blockNums(batFly)) = nanmean(LRVols(sideNum+1:end,:, 2));
            allFlyRight(flyCount, :, 3, blockNums(batFly)) = nanmean(segmentedLRTimes(sideNum+1:end, :));

            %generating SEM for indiividual fly plots
            currentFlyLRSEM = nanstd(LRVols(:,:,1), [], 1)/sqrt(size(LRVols, 1));

            if allGraphs 
                %plotting each individual fly
                % % figure
                % % hold on
                % % plot(segmentedLRTimes', LRVols(:,:,1)')
                % % plot(nanmean(segmentedLRTimes), nanmean(LRVols(:,:,1)), 'LineWidth', 3, 'Color', 'k')
                % % errorbar(nanmean(segmentedLRTimes), nanmean(LRVols(:,:,1)), currentFlyLRSEM, '.', 'Color', 'k')
                % % % plot(nanmean(LRVols(:,:,2))*500 + (nanmean(LRVols(:,1,1))/1.75), 'lineWidth', 3, 'Color', 'b') %redo when decided how to plot photodiode against all
                % % yLims = get(gca, 'ylim');
                % % fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
                % % %photodiode
                % % % plot(segmentedLRTimes', normalize(LRVols(:,:,2)', 'range', [0, max(range(LRVols(:,:,1)))/2]) + max(LRVols(:,:,1), [], 'all') + max(range(LRVols(:,:,1)))/4 ) % raw photos
                % % % plot(nanmean(segmentedLRTimes), nanmean((normalize(LRVols(:,:,2)', 'range', [0, max(range(LRVols(:,:,1)))/2]) + max(LRVols(:,:,1), [], 'all') + max(range(LRVols(:,:,1)))/4)'), 'LineWidth', 3, 'Color', 'b' ) %avg photo
                % %
                % % xlabel('time (s)') %... for now
                % % ylabel('mean fluorescence')
                % % set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                % % set(gcf, 'Name', ['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(blockNums(batFly)), ' calcium trace - ', regionName], 'NumberTitle', 'off', 'Renderer', 'painters')
                % %
                % % title(['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(blockNums(batFly)), ' calcium trace - ', regionName])

                %plot side-separated
                leftSEM = nanstd(LRVols(1:size(LRVols,1)/2, :, 1))/sqrt(size(LRVols, 1)/2);
                rightSEM = nanstd(LRVols(size(LRVols,1)/2 + 1:end, :, 1))/sqrt(size(LRVols, 1)/2);
                %
                % %consecutive
                % p1 = plot(LRVols( [1:size(LRVols,1)/2], :, 3),  LRVols([1:size(LRVols,1)/2], :, 1), 'LineWidth', 3, 'Color', 'b');
                % p2 = plot(LRVols(size(LRVols,1)/2 + 1:end, :, 3), LRVols(size(LRVols,1)/2 + 1:end, :, 1), 'LineWidth', 3, 'Color', 'k');
                % legend([p1,p2], {'left', 'right'})
                %
                % errorbar(LRVols( [1:size(LRVols,1)/2], :, 3),  LRVols([1:size(LRVols,1)/2], :, 1), leftSEM, '.', 'Color','b')
                % errorbar(LRVols(size(LRVols,1)/2 + 1:end, :, 3), LRVols(size(LRVols,1)/2 + 1:end, :, 1), rightSEM, '.', 'Color', 'k')

                %size-sep per stim
                %left
                figure
                hold on
                plot(segmentedLRTimes(1:size(LRVols,1)/2, :)', LRVols(1:size(LRVols,1)/2, :, 1)')
                p1 = plot(nanmean(segmentedLRTimes(1:size(LRVols,1)/2, :)), nanmean(LRVols(1:size(LRVols,1)/2, :, 1)), 'LineWidth', 3, 'Color', 'k');
                % p2 = plot(nanmean(segmentedLRTimes(size(LRVols,1)/2 + 1:end, :)), nanmean(LRVols(size(LRVols,1)/2 + 1:end, :, 1)), 'LineWidth', 3, 'Color', 'k');

                errorbar(nanmean(segmentedLRTimes(1:size(LRVols,1)/2, :)), nanmean(LRVols(1:size(LRVols,1)/2, :, 1)), leftSEM, '.', 'Color','k')
                % errorbar(nanmean(segmentedLRTimes(size(LRVols,1)/2 + 1:end, :)), nanmean(LRVols(size(LRVols,1)/2 + 1:end, :, 1)), rightSEM, '.', 'Color', 'k')
                % legend([p1,p2], {'left', 'right'})

                yLims = get(gca, 'ylim');
                fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
                % legend([p1,p2], {'left', 'right'}) <- this is correct legend

                xlabel('time (s)')
                ylabel('mean fluorescence')
                set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                set(gcf, 'Name', ['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(blockNums(batFly)), ' lefts - ', regionName], 'NumberTitle', 'off', 'Renderer', 'painters')
                title(['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(blockNums(batFly)), ' lefts - ', regionName])

                % title(['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(blockNums(batFly)), ' side-separated calcium trace - ', regionName])

                %right
                figure
                hold on
                p2 = plot(nanmean(segmentedLRTimes(size(LRVols,1)/2 + 1:end, :)), nanmean(LRVols(size(LRVols,1)/2 + 1:end, :, 1)), 'LineWidth', 3, 'Color', 'k');
                errorbar(nanmean(segmentedLRTimes(size(LRVols,1)/2 + 1:end, :)), nanmean(LRVols(size(LRVols,1)/2 + 1:end, :, 1)), rightSEM, '.', 'Color', 'k')
                yLims = get(gca, 'ylim');
                fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
                % legend([p1,p2], {'left', 'right'}) <- this is correct legend

                xlabel('time (s)')
                ylabel('mean fluorescence')
                set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                set(gcf, 'Name', ['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(blockNums(batFly)), ' Rights - ', regionName], 'NumberTitle', 'off', 'Renderer', 'painters')
                title(['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(blockNums(batFly)), ' Rights - ', regionName])
            end
        end
        toc

        %0-mean method
        if baselineCorrect
            for fly = 1:size(allFlyLR, 1)
                for block = 1:size(allFlyLR, 4)
                    %combined
                    correctionFactor = nanmean(allFlyLR(fly, :, 1, block)); %correction factor is just the mean because trying to move baseline to 0 for all flies/blocks
                    allFlyLR(fly, :, 1, block) = allFlyLR(fly, :, 1, block) - correctionFactor;

                    %lefts
                    leftCorrect = nanmean(allFlyLeft(fly, :, 1, block));
                    allFlyLeft(fly, :, 1, block) = allFlyLeft(fly, :, 1, block) - leftCorrect;

                    %rights
                    rightCorrect = nanmean(allFlyRight(fly, :, 1, block));
                    allFlyRight(fly, :, 1, block) = allFlyRight(fly, :, 1, block) - rightCorrect;
                end
            end
        end

        %cross-fly SEM
        LRSEM = nanstd(allFlyLR(:, :, 1, :), [], 1)/sqrt(size(allFlyLR, 1));
        leftSEM = nanstd(allFlyLeft(:,:,1,:), [], 1)/sqrt(size(allFlyLeft, 1));
        rightSEM = nanstd(allFlyRight(:, :, 1, :), [], 1)/sqrt(size(allFlyRight, 1));

        %cross-fly plots
        for batBlock = 1:size(allFlyLR, 4) %find(blockNums == blockType)
            %combined
            figure
            hold on
            plot(allFlyLR(:,:,3, batBlock)', allFlyLR(:,:,1, batBlock)')
            plot(nanmean(allFlyLR(:,:,3, batBlock)), nanmean(allFlyLR(:,:,1, batBlock)), 'LineWidth', 3, 'Color', 'k')
            errorbar(nanmean(allFlyLR(:,:,3, batBlock)), nanmean(allFlyLR(:, :, 1, batBlock)), LRSEM(:,:,batBlock), 'Color', 'k')

            % plotReadyPhot = normalize(allFlyLR(:,:,2, batBlock)', 'range', [0, max(range(allFlyLR(:,:,1, batBlock)))/2]) + max(allFlyLR(:,:,1, batBlock), [], 'all') + max(range(allFlyLR(:,:,1,batBlock)))/4;
            % plot(allFlyLR(:,:,3,batBlock)', plotReadyPhot)
            % plot(nanmean(allFlyLR(:,:,3,batBlock)), nanmean(plotReadyPhot, 2), 'LineWidth', 3, 'Color', 'b')

            xlabel('time (s)')
            ylabel('mean fluorescence')
            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
            % set(gcf, 'Name', ['cross-fly LR traces - ', regionName])
            set(gcf, 'Name', ['cross-fly LR traces - ', regionName, ' block ', num2str(batBlock)], 'NumberTitle', 'off', 'Renderer', 'painters')

            yLims = get(gca, 'ylim');
            fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)

            %lefts
            figure;
            hold on
            plot(allFlyLeft(:,:,3, batBlock)', allFlyLeft(:,:,1, batBlock)')
            plot(nanmean(allFlyLeft(:,:,3, batBlock)), nanmean(allFlyLeft(:,:,1, batBlock)), 'LineWidth', 3, 'Color', 'k')
            errorbar(nanmean(allFlyLeft(:,:,3, batBlock)), nanmean(allFlyLeft(:, :, 1, batBlock)), leftSEM(:,:,batBlock), 'Color', 'k')

            % plotReadyPhot = normalize(allFlyLR(:,:,2, batBlock)', 'range', [0, max(range(allFlyLR(:,:,1, batBlock)))/2]) + max(allFlyLR(:,:,1, batBlock), [], 'all') + max(range(allFlyLR(:,:,1,batBlock)))/4;
            % plot(allFlyLR(:,:,3,batBlock)', plotReadyPhot)
            % plot(nanmean(allFlyLR(:,:,3,batBlock)), nanmean(plotReadyPhot, 2), 'LineWidth', 3, 'Color', 'b')

            xlabel('time (s)')
            ylabel('mean fluorescence')
            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
            % set(gcf, 'Name', ['cross-fly LR traces - ', regionName])
            set(gcf, 'Name', ['cross-fly Left traces - ', regionName, ' block ', num2str(batBlock)], 'NumberTitle', 'off', 'Renderer', 'painters')

            yLims = get(gca, 'ylim');
            fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)

            %right
            figure
            hold on
            plot(allFlyRight(:,:,3, batBlock)', allFlyRight(:,:,1, batBlock)')
            plot(nanmean(allFlyRight(:,:,3, batBlock)), nanmean(allFlyRight(:,:,1, batBlock)), 'LineWidth', 3, 'Color', 'k')
            errorbar(nanmean(allFlyRight(:,:,3, batBlock)), nanmean(allFlyRight(:, :, 1, batBlock)), rightSEM(:,:,batBlock), 'Color', 'k')

            % plotReadyPhot = normalize(allFlyLR(:,:,2, batBlock)', 'range', [0, max(range(allFlyLR(:,:,1, batBlock)))/2]) + max(allFlyLR(:,:,1, batBlock), [], 'all') + max(range(allFlyLR(:,:,1,batBlock)))/4;
            % plot(allFlyLR(:,:,3,batBlock)', plotReadyPhot)
            % plot(nanmean(allFlyLR(:,:,3,batBlock)), nanmean(plotReadyPhot, 2), 'LineWidth', 3, 'Color', 'b')

            xlabel('time (s)')
            ylabel('mean fluorescence')
            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
            % set(gcf, 'Name', ['cross-fly LR traces - ', regionName])
            set(gcf, 'Name', ['cross-fly Right traces - ', regionName, ' block ', num2str(batBlock)], 'NumberTitle', 'off', 'Renderer', 'painters')

            yLims = get(gca, 'ylim');
            fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)

            %photodiode
            % for block = 1:size(blockNums, 2)
            %     figure
            %     hold on
            %     plot(LRVols(:,:,2)') %eventually do against time rather than vols
            %     plot(nanmean(LRVols(:,:,2)), 'LineWidth', 3, 'Color', 'k')
            %     % plot(nanmean(LRVols(:,:,2))*500 + (nanmean(LRVols(:,1,1))/1.75), 'lineWidth', 3, 'Color', 'b') %redo when decided how to plot photodiode against all
            %
            %     xlabel('volumes') %... for now
            %     ylabel('mean flurorescence')
            %     set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
            %
            %     title(['fly ', num2str(batStruct.blocks(batFly).fly), ' block ', num2str(block), ' photodiode'])
            % end
        end

        %amp calcs
        LRAmps = nan(size(allFlyLR, 1), size(allFlyLR, 4), 3); % [flyAmplitude, block, combined(1)/left(2)/right(3), amp(1)/AUC(2)]
        LRAUC = nan(size(LRAmps));
        LRAmpPositions = nan(size(allFlyLR, 1), ampSpecs, size(allFlyLR, 4), 2, 3); % [fly, position, block, peak/trough, combined(1)/left(2)/right(3)] <-- could potentially save peak and trough nums here by extended dim2

        for batBlock = 1:size(allFlyLR, 4) %find(blockNums == blockType)
            for fly = 1:size(allFlyLR, 1)
                for sideType = 1:size(allFlyLR, 3) %maybe move to outer loop to save time if needed
                    %determining which crossFly variable to use for amp calc
                    if sideType == 1
                        currentFlyLR = allFlyLR;
                    elseif sideType == 2
                        currentFlyLR = allFlyLeft;
                    elseif sideType == 3
                        currentFlyLR = allFlyRight;
                    end

                    if isnan(currentFlyLR(fly, 1, 3, batBlock))
                        continue
                    end

                    if ampRestrict
                        if isempty(windowSpecs)
                            windowStart = find(currentFlyLR(fly, :, 3, batBlock) >= 1, 1, 'first');
                            windowEnd = find(~isnan(currentFlyLR(fly, :, 3, batBlock)), 1, 'last'); 
                        else
                            windowStart = find(currentFlyLR(fly, :, 3, batBlock) >= windowSpecs(1), 1, 'first'); %eventually add windowSpecs here
                            windowEnd = find(currentFlyLR(fly, :, 3, batBlock) >= windowSpecs(2), 1, 'first');
                        end
                        
                        [~,minPos] = min(currentFlyLR(fly, windowStart+3:windowEnd-3, 1, batBlock));

                        minPos = minPos + windowStart+3 - 1; %correcting for size change due to restricted window size

                        % [peakValues, peakIndex] = maxk(currentFlyLR(fly, windowStart+3:windowEnd-3, 1, batBlock), ampSpecs);
                        [peakValues, peakIndex] = maxk(currentFlyLR(fly, windowStart+3:minPos, 1, batBlock), ampSpecs);
                        [troughValues, troughIndex] = mink(currentFlyLR(fly, windowStart+3:windowEnd-3, 1, batBlock), ampSpecs);

                        if AUC
                            LRAUC(fly, batBlock, sideType) = trapz(squeeze(currentFlyLR(fly, windowStart+3:windowEnd-3, 3, batBlock)), squeeze(currentFlyLR(fly, windowStart+3:windowEnd-3, 1, batBlock)), 2);
                        end
                    else
                        [peakValues, peakIndex] = maxk(currentFlyLR(fly, :, 1, batBlock), ampSpecs);
                        [troughValues, troughIndex] = mink(currentFlyLR(fly, :, 1, batBlock), ampSpecs);
                    end

                    % [peakValues, peakIndex] = maxk(currentFlyLR(fly, :, 1, batBlock), ampSpecs);
                    % [troughValues, troughIndex] = mink(currentFlyLR(fly, :, 1, batBlock), ampSpecs);

                    LRAmpPositions(fly, :, batBlock, 1) = currentFlyLR(fly, peakIndex, 3, batBlock);
                    LRAmpPositions(fly, :, batBlock, 2) = currentFlyLR(fly, troughIndex, 3, batBlock);

                    currentPeak = nanmean(peakValues);
                    currentTrough = nanmean(troughValues);

                    % LRAmps(fly, [1,2], batBlock) = [currentPeak, currentTrough];
                    if useMaxOnly
                        LRAmps(fly, batBlock, sideType) = currentPeak;
                    else
                        LRAmps(fly, batBlock, sideType) = currentPeak - currentTrough;
                    end
                end
            end
        end
        %rudimentary plot of amplitudes for LR
        % figure
        % bar([1,2,3], [[nanmean(LRAmps(:,1))], [nanmean(LRAmps(:,2))], [nanmean(LRAmps(:,3))]])

    end
    %% lum
    if conditionChoice(3)
        lumColours = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]}; %list of colours so that individual lum plots match collated plots
        %traces for each fly
        allFlyLum = nan(size(chosenFlies, 2), size(lumVols, 2), size(lumVols, 3), size(lumVols, 4), max(blockNums));
        individualLumSEM = nan(size(chosenFlies, 2), size(lumVols, 2), size(lumVols, 4), max(blockNums));
        flyCount = 1;
        currentFly = batStruct.blocks(1).fly;
        for block = 1:size(blockNums, 2) %find(blockNums == 1)
            if batStruct.blocks(block).fly ~= currentFly
                currentFly = batStruct.blocks(block).fly;
                flyCount = flyCount + 1;
            end

            lumVols = batStruct.blocks(block).lumVols;
            segmentedLumTimes = batStruct.blocks(block).segmentedLumTimes;

            if allGraphs
                figure
                set(gcf, 'Name', ['Fly ', num2str(batStruct.blocks(block).fly), 'block ', num2str(blockNums(block)), ': Luminosities'], 'NumberTitle', 'off', 'Renderer', 'painters')
            end
            for lum = 1:size(lumVols, 4)
                allFlyLum(flyCount, :, 1, lum, blockNums(block)) = nanmean(lumVols(:, :, 1, lum));
                allFlyLum(flyCount, :, 2, lum, blockNums(block)) = nanmean(lumVols(:, :, 2, lum));
                allFlyLum(flyCount, :, 3, lum, blockNums(block)) = nanmean(segmentedLumTimes(:, :, lum));

                individualLumSEM(flyCount, :, lum, blockNums(block)) = nanstd(lumVols(:, :, 1, lum), [], 1)/sqrt(size(lumVols, 1));

                if allGraphs
                    %trace
                    subplot(2, size(lumVols, 4), lum)
                    hold on
                    plot(segmentedLumTimes(:,:,lum)', lumVols(:,:,1,lum)') %probably going to need to plot against a normalised time eventually
                    plot(nanmean(segmentedLumTimes(:,:,lum)), nanmean(lumVols(:,:,1,lum)), 'LineWidth', 3, 'Color', 'k')
                    yLims = get(gca, 'ylim');
                    % fill([0, 0, 14, 14], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3) % based on vols for fly 85
                    fill([0, 0, 1, 1], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3) %based on assumed duration in seconds - is wrong but may just be out of phase
                    title(['Luminosity=', num2str(round(lum/size(lumVols, 4), 2))])

                    %photodiode
                    subplot(2,size(lumVols, 4), lum + size(lumVols, 4))
                    hold on
                    plot(segmentedLumTimes(:,:,lum)', lumVols(:,:,2,lum)')
                    plot(nanmean(segmentedLumTimes(:,:,lum)), nanmean(lumVols(:,:,2,lum)), 'LineWidth', 3, 'Color', 'k')
                    yLims = get(gca, 'ylim');
                    fill([0, 0, 1, 1], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
                    % title(num2str(round(1/lum, 2)))
                end
            end
        end

        %0-mean - DONT DO HERE (OBLITERATES LUMINESS)
        if baselineCorrect
            for fly = 1:size(allFlyLum, 1)
                for block = 1:size(allFlyLum, 5)
                    for lum = 1:size(allFlyLum, 4)
                        correctionFactor = nanmean(allFlyLum(fly, :, 1, lum, block));
                        allFlyLum(fly, :, 1, lum, block) = allFlyLum(fly, :, 1, lum, block) - correctionFactor;
                    end
                end
            end
        end


        % correctionFactor = nanmean(allFlyLR(fly, :, 1, block)); %correction factor is just the mean because trying to move baseline to 0 for all flies/blocks
        % allFlyLR(fly, :, 1, block) = allFlyLR(fly, :, 1, block) - correctionFactor;


        if allGraphs
            %plotting individuals
            % lumColours = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
            lumTitles = {};
            for f = 1:size(allFlyLum, 1)
                for block = 1% 1:size(allFlyLum, 5)
                    if isnan(allFlyLum(f, 1, 2, lum, block))
                        continue %ensuring only necessary plots are made
                    end
                    figure
                    hold on
                    for lum = 1:7
                        plot(squeeze(allFlyLum(f, :, 3, lum, block)), squeeze(allFlyLum(f, :, 1, lum, block)), 'LineWidth', 3)
                        % errorbar( squeeze( allFlyLum(f, :, 3, lum, block)), squeeze(allFlyLum(f, :, 1, lum, block)), individualLumSEM(f, :, lum, block), '.')

                        ylabel('mean fluorescence (\DeltaF/F)')
                        xlabel('time (s)')
                        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                        set(gcf, 'Name', ['fly ', num2str(chosenFlies(f)), ' block ', num2str(block), ' - all luminosity traces'], 'Renderer', 'painters')
                        legend()
                    end
                end
            end

            ylims = get(gca, 'ylim'); %pretty janky atm, just sets limits (of fill - for stimulus periods) based on last fly's collective plot
            ylims(2) = ylims(2) + ylims(2)/10; %arbitrary addition just to get code to function - endgame might involve a ylims variable that is size of flyNum or blockNum

            %plot each luminosity on individual plot with errorbars - fairly hardcoded for fly 86 rn
            %maybe move above into collectionm loop - because ylims will no longer work (given multi-block functionality)
            % ^^ otherwise as a for fly loop
            for block = 1:size(allFlyLum, 5)
                for lum = 1:7
                    figure
                    hold on
                    plot(squeeze(allFlyLum(1, :, 3, lum, block)), squeeze(allFlyLum(1, :, 1, lum, block)), 'LineWidth', 3, 'Color', lumColours{lum})
                    errorbar(squeeze(allFlyLum(1, :, 3, lum, block)), squeeze(allFlyLum(1, :, 1, lum, block)), individualLumSEM(1, :, lum, block), '.', 'Color', lumColours{lum})
                    fill([0, 0, 1.08938, 1.08938], [ylims(1), ylims(2), ylims(2), ylims(1)], 'b', 'FaceAlpha', 0.3) %ylims was hardcoded from collective plot

                    ylabel('mean fluorescence (\DeltaF/F)')
                    xlabel('time (s)')
                    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                    % title( ['fly ', num2str() num2str(round(lum/7, 2)), ' - EPG'] )
                end
            end
        end

        %cross-fly lum analysis
        lumSEM = squeeze(nanstd(allFlyLum(:,:,1,:,:), [], 1)/sqrt(size(allFlyLum, 1)));
        %specific lum plots and amplitude calculation
        lumAmps = nan(size(allFlyLum, [1,4,5])); % [flyAmplitude, lum, block]
        lumAUC = nan(size(lumAmps));
        lumAmpPositions = nan(size(allFlyLum, 1), ampSpecs, size(allFlyLum, 4), size(allFlyLum, 5), 2); %[fly, position, lum, block, peak/trough]
        for block = 1:size(allFlyLum, 5)
            for lum = 1:size(allFlyLum, 4)
                %maybe plot if statement here
                % figure
                % hold on
                % plot(allFlyLum(:,:,3,lum,block)', allFlyLum(:,:,1,lum,block)', 'Color', lumColours{lum})
                % plot(nanmean(allFlyLum(:, :, 3, lum, block)), nanmean(allFlyLum(:,:,1,lum,block)), 'LineWidth', 3, 'Color', lumColours{lum})
                % errorbar(nanmean(allFlyLum(:, :, 3, lum, block)), nanmean(allFlyLum(:,:,1,lum,block)), lumSEM(:,lum,block), '.', 'Color', lumColours{lum})

                % insert spiffy plot specs here

                for fly = 1:size(allFlyLum, 1)
                    % if LRBasedAmps
                    %     windowStart = find(allFlyLum(fly, :, 3, lum, block) >= 1, 1, 'first');
                    %     [peakValues, peakIndex] = maxk(allFlyLum(fly, windowStart+3:end-3, 1, lum, block), ampSpecs);
                    %     [troughValues, troughIndex] = mink(allFlyLum(fly, windowStart+3:end-3, 1, lum, block), ampSpecs);

                    if isnan(allFlyLum(fly, 1, 3, lum, block))
                        continue
                    end

                    if ampRestrict
                        if isempty(windowSpecs)
                            windowStart = find(allFlyLum(fly, :, 3, lum, block) >= 1, 1, 'first');
                            windowEnd = find(~isnan(allFlyLum(fly, :, 3, lum, block)), 1, 'last');
                        else
                            windowStart = find(allFlyLum(fly, :, 3, lum, block) >= windowSpecs(1), 1, 'first');
                            windowEnd = find(allFlyLum(fly, :, 3, lum, block) >= windowSpecs(2), 1, 'first');
                        end

                        [~,minPos] = min(allFlyLum(fly, windowStart+3:windowEnd-3, 1, lum, block));

                        minPos = minPos + windowStart+3 - 1; %correcting for size change due to restricted window size

                        % [peakValues, peakIndex] = maxk(allFlyLum(fly, windowStart+3:windowEnd-3, 1, lum, block), ampSpecs);
                        [peakValues, peakIndex] = maxk(allFlyLum(fly, windowStart+3:minPos, 1, lum, block), ampSpecs);

                        [troughValues, troughIndex] = mink(allFlyLum(fly, windowStart+3:windowEnd-3, 1, lum, block), ampSpecs);
                        if AUC
                            lumAUC(fly, lum, block) = trapz(squeeze(allFlyLum(fly, windowStart+3:windowEnd-3, 3, lum, block)), squeeze(allFlyLum(fly, windowStart+3:end-3, 1, lum, block)), 2);
                        end
                    else
                        [peakValues, peakIndex] = maxk(allFlyLum(fly, :, 1, lum, block), ampSpecs);
                        [troughValues, troughIndex] = mink(allFlyLum(fly, :, 1, lum, block), ampSpecs);
                    end

                    lumAmpPositions(fly, :, lum, block, 1) = allFlyLum(fly, peakIndex, 3, lum, block);
                    lumAmpPositions(fly, :, lum, block, 2) = allFlyLum(fly, troughIndex, 3, lum, block);

                    currentPeak = nanmean(peakValues);
                    currentTrough = nanmean(troughValues);

                    % LRAmps(fly, [1,2], batBlock) = [currentPeak, currentTrough];
                    if useMaxOnly
                        lumAmps(fly, lum, block) = currentPeak;
                    else
                        lumAmps(fly, lum, block) = currentPeak - currentTrough;
                    end
                end
            end
        end

        for block = 1:size(allFlyLum, 5)
            figure
            hold on
            for lum = 1:size(allFlyLum, 4)
                plot(nanmean(allFlyLum(:,:,3,lum,block)), nanmean(allFlyLum(:,:,1,lum,block)), 'LineWidth', 3, 'Color', lumColours{lum})
                set(gcf, 'Name', ['cross-fly luminosities, block ', num2str(block)], 'Renderer', 'painters')
            end
        end
    end
    %% freq
    if conditionChoice(4)
        allFlyFreq = nan(size(chosenFlies, 2), size(freqVols, 2), size(freqVols, 3), size(freqVols, 4), max(blockNums));
        individualFreqSEM = nan(size(chosenFlies, 2), size(freqVols, 2), size(freqVols, 4), max(blockNums));
        flyCount = 1;
        currentFly = batStruct.blocks(1).fly;
        extensivePlots = 1; % option to plot all bg traces and photodiode (put earlier in code and extend it to other conditions)

        for block = 1:size(blockNums,2)%find(blockNums == 2)
            if batStruct.blocks(block).fly ~= currentFly
                currentFly = batStruct.blocks(block).fly;
                flyCount = flyCount + 1;
            end
            freqVols = batStruct.blocks(block).freqVols;
            segmentedFreqTimes = batStruct.blocks(block).segmentedFreqTimes;

            for freq = 1:size(freqVols, 4)
                %gathering mean for each block
                allFlyFreq(flyCount, :, 1, freq, blockNums(block)) = nanmean(freqVols(:,:,1,freq));
                allFlyFreq(flyCount, :, 2, freq, blockNums(block)) = nanmean(freqVols(:,:,2,freq));
                allFlyFreq(flyCount, :, 3, freq, blockNums(block)) = nanmean(segmentedFreqTimes(:,:,freq));

                individualFreqSEM(flyCount, :, freq, blockNums(block)) = nanstd(freqVols(:, :, 1, freq), [], 1)/sqrt(size(freqVols, 1));

                %add extensive plots here
            end
        end

        %matching length of crossFly frequencies
        %gathering sizes
        flyFreqSizes = nan(size(allFlyFreq, 1), size(allFlyFreq, 4), size(allFlyFreq, 5)); %[ fly, freq, block ]
        for block = 1:size(flyFreqSizes, 3)
            for f = 1:size(allFlyFreq, 1)
                for freq = 1:size(flyFreqSizes, 2)
                    if isnan(allFlyFreq(f, 1, 3, freq, block))
                        continue
                    end
                    flyFreqSizes(f, freq, block) = find(~isnan(allFlyFreq(f,:,1,freq,block)), 1, 'last');
                end
            end
        end

        %trimming based on most consistent size
        flyFreqSizes = mode(flyFreqSizes, 1);
        for block = 1:size(flyFreqSizes, 3)
            for f = 1:size(allFlyFreq, 1)
                for freq = 1:size(allFlyFreq, 4)
                    if isnan(allFlyFreq(f, 1, 3, freq, block))
                        continue
                    end
                    allFlyFreq(f, flyFreqSizes(1, freq, block) + 1:end, :, freq, block) = nan;
                end
            end
        end

        % 0-mean - potentially move between crossFly collective plot and individual plots
        %^^ is 0-mean best - would setting start to 0 be better???
        if baselineCorrect
            for fly = 1:size(allFlyFreq, 1)
                for block = 1:size(allFlyFreq, 5)
                    for freq = 1:size(allFlyFreq, 4)
                        correctionFactor = nanmean(allFlyFreq(fly, :, 1, freq, block));
                        allFlyFreq(fly, :, 1, freq, block) = allFlyFreq(fly, :, 1, freq, block) - correctionFactor;
                    end
                end
            end
        end

        freqStimEstimate = [0.454, 0.8171, 0.9079, 1.089, 1.725]; %empirical stimulus cut-off times for each frequency

        if allGraphs
            %individual fly plots
            for f = 1:size(allFlyFreq, 1) %:size(chosenFlies, 2)
                for block = 1:size(allFlyFreq,5)
                    figure
                    hold on
                    for freq = 1:size(allFlyFreq, 4)
                        plot(squeeze(allFlyFreq(f, :, 3, freq, block)), squeeze(allFlyFreq(f, :, 1, freq, block)), 'LineWidth', 3)
                        % errorbar( squeeze( allFlyFreq(f, :, 3, freq)), squeeze(allFlyFreq(f, :, 1, freq)), individualFreqSEM(f, :, freq), '.')

                        ylabel('mean fluorescence (\DeltaF/F)')
                        xlabel('time (s)')
                        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                        set(gcf, 'Name', ['fly ', num2str(chosenFlies(f)), ' block ', num2str(block), ' - all frequencies (', regionName, ')'], 'Renderer', 'painters')

                        legend({'1Hz', '3Hz', '6Hz', '12Hz', '25Hz'})
                    end
                end
            end

            ylims = get(gca, 'ylim');

            %individual frequency plots
            for fly = 1:size(allFlyFreq, 1)
                for block = 1:size(allFlyFreq, 5)
                    for freq = 1:size(allFlyFreq, 4)
                        figure
                        hold on
                        plot(squeeze(allFlyFreq(fly, :, 3, freq, block)), squeeze(allFlyFreq(fly, :, 1, freq, block)), 'LineWidth', 3, 'Color', lumColours{freq})
                        errorbar(squeeze(allFlyFreq(fly, :, 3, freq, block)), squeeze(allFlyFreq(fly, :, 1, freq, block)), individualFreqSEM(1, :, freq, block), '.', 'Color', lumColours{freq})
                        fill([0, 0, freqStimEstimate(freq), freqStimEstimate(freq)], [ylims(1), ylims(2), ylims(2), ylims(1)], 'b', 'FaceAlpha', 0.3) %ylims was hardcoded from collective plot

                        ylabel('mean fluorescence (\DeltaF/F)')
                        xlabel('time (s)')
                        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                        set(gcf, 'Name', ['fly ', num2str(chosenFlies(fly)), ' block ', num2str(block), ' - ', num2str(matSave.batFrequencies(freq)), ' Hz (', regionName, ')'], 'Renderer', 'painters')
                        title( ['fly ', num2str(chosenFlies(fly)), ' block ', num2str(block), ' - ', num2str(matSave.batFrequencies(freq)), ' Hz'])
                    end
                end
            end
        end

        % cross fly frequency
        freqSEM = squeeze(nanstd(allFlyFreq(:,:,1,:,:), [], 1)/sqrt(size(allFlyFreq, 1)));
        freqAmps = nan(size(allFlyFreq, [1,4,5])); % [flyAmplitude, freq, block]
        freqAUC = nan(size(freqAmps));
        freqAmpPositions = nan(size(allFlyFreq, 1), ampSpecs, size(allFlyFreq, 4), size(allFlyFreq, 5), 2); %[fly, position, freq, block, peak/trough]

        for block = 1:size(allFlyFreq, 5)
            for freq = 1:size(allFlyFreq, 4)
                if allGraphs
                    figure
                    hold on
                    plot(squeeze(allFlyFreq(:, :, 3, freq, block))', squeeze(allFlyFreq(:, :, 1, freq, block))', 'Color', lumColours{freq})
                    plot(nanmean(allFlyFreq(:,:,3,freq,block)), nanmean(allFlyFreq(:,:,1,freq,block)), 'LineWidth', 3, 'Color', lumColours{freq})
                    errorbar(nanmean(allFlyFreq(:,:,3,freq,block)), nanmean(allFlyFreq(:,:,1,freq,block)), freqSEM(:, freq, block), '.', 'Color', lumColours{freq})
                    % fill([0, 0, freqStimEstimate(freq), freqStimEstimate(freq)], [ylims(1), ylims(2), ylims(2), ylims(1)], 'b', 'FaceAlpha', 0.3)

                    ylabel('mean fluorescence (\DeltaF/F)')
                    xlabel('time (s)')
                    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
                    set(gcf, 'Name', ['cross-fly block ', num2str(block), ' - ', num2str(matSave.batFrequencies(freq)), ' Hz (', regionName, ')'], 'Renderer', 'painters')
                    title( ['cross-fly ', ' block ', num2str(block), ' - ', num2str(matSave.batFrequencies(freq)), ' Hz'])
                end

                for fly = 1:size(allFlyFreq, 1)
                    if isnan(allFlyFreq(fly, 1, 3, freq, block))
                        continue
                    end

                    if ampRestrict
                        if isempty(windowSpecs)
                            windowStart = find(allFlyFreq(fly, :, 3, freq, block)>=1.27, 1, 'first');
                            windowEnd = find(~isnan(allFlyFreq(fly, :, 3, freq, block)), 1, 'last');
                        else
                            windowStart = find(allFlyFreq(fly, :, 3, freq, block)>=windowSpecs(1), 1, 'first');
                            windowEnd = find(allFlyFreq(fly, :, 3, freq, block)>=windowSpecs(2), 1, 'first');
                        end

                        [~,minPos] = min(allFlyFreq(fly, windowStart+3:windowEnd-3, 1, freq, block));

                        minPos = minPos + windowStart+3 - 1; %correcting for size change due to restricted window size

                        % [peakValues, peakIndex] = maxk(allFlyFreq(fly, windowStart+3:windowEnd-3, 1, freq, block), ampSpecs);
                        [peakValues, peakIndex] = maxk(allFlyFreq(fly, windowStart+3:minPos, 1, freq, block), ampSpecs);

                        [troughValues, troughIndex] = mink(allFlyFreq(fly, windowStart+3:windowEnd-3, 1, freq, block), ampSpecs);
                        if AUC
                            freqAUC(fly, freq, block) = trapz(squeeze(allFlyFreq(fly, windowStart+3:windowEnd-3, 3, freq, block)), squeeze(allFlyFreq(fly, windowStart+3:windowEnd-3, 1, freq, block)), 2);
                        end
                    else
                        [peakValues, peakIndex] = maxk(allFlyFreq(fly, :, 1, freq, block), ampSpecs);
                        [troughValues, troughIndex] = mink(allFlyFreq(fly, :, 1, freq, block), ampSpecs);
                    end

                    freqAmpPositions(fly, :, freq, block, 1) = allFlyFreq(fly, peakIndex, 3, freq, block);
                    freqAmpPositions(fly, :, freq, block, 2) = allFlyFreq(fly, troughIndex, 3, freq, block);

                    currentPeak = nanmean(peakValues);
                    currentTrough = nanmean(troughValues);

                    % LRAmps(fly, [1,2], batBlock) = [currentPeak, currentTrough];
                    if useMaxOnly
                        freqAmps(fly, freq, block) = currentPeak;
                    else
                        freqAmps(fly, freq, block) = currentPeak - currentTrough;
                    end
                end
            end
        end

        for block = 1:size(allFlyFreq, 5)
            figure
            hold on
            for freq = [2:4] % 1:size(allFlyFreq, 4)
                % plot(nanmean(allFlyFreq(:,:,3,freq,block)), nanmean(allFlyFreq(:,:,1,freq,block)), 'Color', lumColours{freq})
                plot(nanmean(allFlyFreq(:,:,3,freq,block)), nanmean(allFlyFreq(:,:,1,freq,block)), 'LineWidth', 3)
            end
            %figure sets
            xlabel('time(s)')
            ylabel('mean fluorescence')
            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
            set(gcf, 'Name', ['cross-fly block ', num2str(block), ' - ', ' - all frequencies (', regionName, ')'], 'Renderer', 'painters')
            title('Cross-fly frequencies - block ', num2str(block))
        end
    end

    %% blanks
    %no conditionChoice because blanks likely needed for all conditions
    allFlyBlanks = nan(size(chosenFlies, 2), size(blankVols, 2), size(blankVols, 3), max(blockNums));
    individualBlankSEM = nan(size(chosenFlies, 2), size(blankVols, 2), max(blockNums));
    flyCount = 1;
    currentFly = batStruct.blocks(1).fly;

    for block = 1:size(blockNums, 2)
        if batStruct.blocks(block).fly ~= currentFly
            currentFly = batStruct.blocks(block).fly;
            flyCount = flyCount + 1;
        end
        blankVols = batStruct.blocks(block).blankVols;
        segmentedBlankTimes = batStruct.blocks(block).segmentedBlankTimes;

        allFlyBlanks(flyCount, :, 1, blockNums(block)) = nanmean(blankVols(:,:,1));
        allFlyBlanks(flyCount, :, 2, blockNums(block)) = nanmean(blankVols(:,:,2));
        allFlyBlanks(flyCount, :, 3, blockNums(block)) = nanmean(segmentedBlankTimes);

        individualBlankSEM(flyCount, :, blockNums(block)) = nanstd(blankVols(:, :, 1), [], 1)/sqrt(size(blankVols, 1));

        %individual fly plots
    end

    %0-mean method
    if baselineCorrect
        for fly = 1:size(allFlyBlanks, 1)
            for block = 1:size(allFlyBlanks, 4)
                if isnan(allFlyBlanks(fly, 1, 3, block))
                    continue
                end
                correctionFactor = nanmean(allFlyBlanks(fly, :, 1, block)); %correction factor is just the mean because trying to move baseline to 0 for all flies/blocks
                allFlyBlanks(fly, :, 1, block) = allFlyBlanks(fly, :, 1, block) - correctionFactor;
            end
        end
    end

    %cross-fly SEM
    blankSEM = squeeze( nanstd(allFlyBlanks(:, :, 1, :), [], 1)/sqrt(size(allFlyBlanks, 1)));

    %true cross-fly LR Analysis
    blankAmps = nan(size(allFlyBlanks, 1), size(allFlyBlanks, 4)); % [flyAmplitude, block]
    blankAmpPositions = nan(size(allFlyBlanks, 1), ampSpecs, size(allFlyBlanks, 4), 2); % [fly, position, block, peak/trough] <-- could potentially save peak and trough nums here by extended dim2

    for batBlock = 1:size(allFlyBlanks, 4) %find(blockNums == blockType)
        figure
        hold on
        plot(allFlyBlanks(:,:,3, batBlock)', allFlyBlanks(:,:,1, batBlock)')
        plot(nanmean(allFlyBlanks(:,:,3, batBlock)), nanmean(allFlyBlanks(:,:,1, batBlock)), 'LineWidth', 3, 'Color', 'k')
        errorbar(nanmean(allFlyBlanks(:,:,3, batBlock)), nanmean(allFlyBlanks(:, :, 1, batBlock)), LRSEM(:,batBlock), 'Color', 'k')

        plotReadyPhot = normalize(allFlyBlanks(:,:,2, batBlock)', 'range', [0, max(range(allFlyBlanks(:,:,1, batBlock)))/2]) + max(allFlyBlanks(:,:,1, batBlock), [], 'all') + max(range(allFlyBlanks(:,:,1,batBlock)))/4;
        plot(allFlyBlanks(:,:,3,batBlock)', plotReadyPhot)
        plot(nanmean(allFlyBlanks(:,:,3,batBlock)), nanmean(plotReadyPhot, 2), 'LineWidth', 3, 'Color', 'b')

        xlabel('time (s)')
        ylabel('mean fluorescence')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        % set(gcf, 'Name', ['cross-fly LR traces - ', regionName])
        set(gcf, 'Name', ['cross-fly blank traces - ', regionName, 'block ', num2str(batBlock)], 'NumberTitle', 'off', 'Renderer', 'painters')

        for fly = 1:size(allFlyBlanks, 1)
            if isnan(allFlyBlanks(fly, 1, 3, batBlock))
                continue
            end

            [peakValues, peakIndex] = maxk(allFlyBlanks(fly, :, 1, batBlock), ampSpecs);
            [troughValues, troughIndex] = mink(allFlyBlanks(fly, :, 1, batBlock), ampSpecs);

            blankAmpPositions(fly, :, batBlock, 1) = allFlyBlanks(fly, peakIndex, 3, batBlock);
            blankAmpPositions(fly, :, batBlock, 2) = allFlyBlanks(fly, troughIndex, 3, batBlock);

            currentPeak = nanmean(peakValues);
            currentTrough = nanmean(troughValues);

            if useMaxOnly
                blankAmps(fly, batBlock) = currentPeak;
            else
                blankAmps(fly, batBlock) = currentPeak - currentTrough;
            end
        end
    end


    %% amplitude analyses
    %note:these figures will likely only be explanatory - not for report figure, need to add saving or cross ROI script for this
    %LR
    LRAmpSEM = squeeze(nanstd(LRAmps, [], 1)/sqrt(size(LRAmps, 1)));
    if allGraphs
        for block = 1:size(LRAmps, 2)
            figure
            hold on
            bar([1,2,3], [squeeze(nanmean(LRAmps(:,block, :), 1))])
            % errorbar([1,2], [nanmean(LRAmps(:,:,2), 1), nanmean(LRAmps(:,:,3), 1)], squeeze(LRAmpSEM(:,:,[2,3])), 'LineStyle', 'none')
            errorbar([1,2,3], [squeeze(nanmean(LRAmps(:,block, :), 1))], LRAmpSEM(block, :), 'LineStyle', 'none', 'Color', 'k')
            xticks([1,2,3])
            xticklabels({'Combined', 'Left', 'Right'})
            title('Cross-Fly LR amplitudes - block ', num2str(block))
        end
    end

    %lum
    lumAmpSEM = nanstd(lumAmps, [], 1)/sqrt(size(lumAmps, 1));
    if allGraphs
        for block = 1:size(lumAmps, 3)
            figure
            hold on
            bar([1:size(lumAmps, 2)], [nanmean(lumAmps(:,:,block), 1)])
            errorbar([1:size(lumAmps, 2)], [nanmean(lumAmps(:,:,block), 1)], squeeze(lumAmpSEM(:,:,block)), 'LineStyle', 'none', 'Color', 'k')
            title('Cross-Fly lum amplitudes - block ', num2str(block))
        end
    end

    %freq
    freqAmpSEM = nanstd(freqAmps, [], 1)/sqrt(size(freqAmps, 1));
    if allGraphs
        for block = 1:size(freqAmps, 3)
            figure
            hold on
            bar([1:size(freqAmps, 2)], [nanmean(freqAmps(:,:,block), 1)])
            errorbar([1:size(freqAmps, 2)], [nanmean(freqAmps(:,:,block), 1)], squeeze(freqAmpSEM(:,:,block)), 'LineStyle', 'none', 'Color', 'k')
            xticks([1:size(freqAmps, 2)])
            xticklabels({'1Hz', '3Hz', '6Hz', '12Hz', '25Hz'})
            title('Cross-Fly freq amplitudes - block ', num2str(block))
        end
    end

    %% saving crossFly data for crossROI analyses
    overStruct(region).regionName = regionName;
    %haven't bothered with opto because not included in thesis
    %LR
    overStruct(region).LR = allFlyLR;
    overStruct(region).LRSEM = LRSEM;
    overStruct(region).lefts = allFlyLeft;
    overStruct(region).leftSEM = leftSEM;
    overStruct(region).rights = allFlyRight;
    overStruct(region).rightSEM = rightSEM;
    overStruct(region).LRAmps = LRAmps;
    overStruct(region).LRAUC = LRAUC;
    %lum
    overStruct(region).lum = allFlyLum;
    overStruct(region).lumSEM = lumSEM;
    overStruct(region).lumAmps = lumAmps;
    overStruct(region).lumAUC = lumAUC;
    %freq
    overStruct(region).freq = allFlyFreq;
    overStruct(region).freqSEM = freqSEM;
    overStruct(region).freqAmps = freqAmps;
    overStruct(region).freqAUC = freqAUC;

end %roi loop end

%% Cross-ROI analyses
if size(roiList, 2) > 1 
    disp('Cross fly analysis complete! commencing cross ROI analysis')
else
    disp('Cross fly analysis complete!')
end
%% amps
%collecting Amps
regionisedAmps = nan(size(LRAmps, 1), size(lumAmps, 2), size(LRAmps, 2), 3, size(roiList, 2)); %[ fly, conditionFeatures, block, LR(1)/lum(2)/freq(3), ROI ] <- assumes last region+block through loop can be representative of others (pretty sound assumption I believe)

for roi = 1:size(overStruct, 2)
    regionLR = overStruct(roi).LRAmps;
    regionLum = overStruct(roi).lumAmps;
    regionFreq = overStruct(roi).freqAmps;

    % regionLR = overStruct(roi).LRAUC;
    % regionLum = overStruct(roi).lumAUC;
    % regionFreq = overStruct(roi).freqAUC;

    % regionisedAmps(:,1:size(regionLR, 3),:,1,roi) = LRAmps;
    for block = 1:size(regionisedAmps, 3)
        regionisedAmps(:,1:size(regionLR, 3),block,1,roi) = regionLR(:,block,:);
        regionisedAmps(:,:,block,2,roi) = regionLum(:,:,block);
        regionisedAmps(:,1:size(regionFreq, 2),block,3,roi) = regionFreq(:,:,block);
        % regionisedAmps(:,:,block,1,roi) = overStruct(roi).LRAmps(:,block,:);
        % regionisedAmps(:,:,block,2,roi) = overStruct(roi).lumAmps(:,:,block);
        % regionisedAmps(:,:,block,3,roi) = overStruct(roi).freqAmps(:,:,block);
    end
end

regionisedAmpSEM = nanstd(regionisedAmps, [], 1)/sqrt(size(regionisedAmps, 1));

%plotting - will currently make plots comparing amps for each ROI in order of roiList
%LR
%comparing ROIs
for block = 1:size(regionisedAmps, 3)
    figure
    hold on
    % bar( [1:size(regionisedAmps, 5)]', squeeze(nanmean(regionisedAmps(:,:,block,1,:), 1))' )
    % bar( repmat(1:find(~isnan(regionisedAmps(1,:,block,1,1)), 1, 'last'), size(regionisedAmps, 5), 1), squeeze(nanmean(regionisedAmps(:,:,block,1,:), 1))' )
    % bar( repmat(1:3, size(regionisedAmps, 5), 1), squeeze(nanmean(regionisedAmps(:,1:3,block,1,:), 1))' ) % features of left rigth condition hardcoded (safe assumption)
    bar( [1:size(regionisedAmps, 5)]', squeeze(nanmean(regionisedAmps(:,1:3,block,1,:), 1))' ) % features of left rigth condition hardcoded (safe assumption)
    %errorbar
    xticks(1:size(regionisedAmps, 5))
    xticklabels(roiList)
    %spiffy sets
    title(['Cross Fly LR amps - block ' num2str(block)])
end

%roiBased compas
RoiOffset = [-0.2727   -0.0909    0.0909    0.2727]; %hardcoded for sep of 4 bars
RoiOffset = RoiOffset(1:numel(roiList));
for block = 1:size(regionisedAmps, 3)
    % figure
    % hold on
    for stimType = [2,3]
        figure
        hold on
        %stats
        LRData = squeeze(regionisedAmps(:, stimType, block, 1, :));

        nonNormCheck = 0;
        for col = 1:size(LRData, 2)
            normTest = lillietest(LRData(:,col));
            if ~normTest
                nonNormCheck = 1;
            end
        end

        if nonNormCheck
            % disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' LR'])
            % [~,~,LRStats] = kruskalwallis(LRData, [] , 'off');
            [~,~,LRStats] = friedman(LRData, 1 , 'off');
            % close, close %closing obligatory figures that anova makes
            LRSig = multcompare(LRStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            % LRDescription = ['kruskal-wallis test', newline];
            LRDescription = ['Friedman test', newline];
        else
            % [~,~,LRStats] = anova1(LRData, [] , 'off');
            [~,~,LRStats] = ranova(LRData, [] , 'off');
            % close, close %closing obligatory figures that anova makes
            LRSig = multcompare(LRStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            % LRDescription = ['one-way anova', newline];
            LRDescription = ['repeated-measures anova', newline];
        end

        % LRDescription = [];
        for group = 1:size(LRSig, 1)
            LRDescription = [LRDescription, [num2str(LRSig(group,1)), 'v', num2str(LRSig(group,2)), ': ', num2str(LRSig(group,6)), newline]];
        end
            
        %plotting
        bar([1:size(regionisedAmps, 5)], squeeze(nanmean(regionisedAmps(:, stimType, block, 1, :), 1)) )
        barX = RoiOffset + (stimType - 1);
        errorbar([1:size(regionisedAmps, 5)], squeeze(nanmean(regionisedAmps(:, stimType, block, 1, :), 1)), squeeze(regionisedAmpSEM(:,stimType, block, 1,:)), 'LineStyle', 'none', 'Color', 'k')
        if doDots
            scatter([1:size(regionisedAmps, 5)], squeeze(regionisedAmps(:, stimType, block, 1, :)), 10, 'k', 'filled')
        end

        %displaying stats
        figMax = get(gca, 'ylim');
        figMax(1) = [];
        text(barX(end), figMax/2, LRDescription)

        ylabel('Amplitude (peak\DeltaF/F))')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'ylim', [-0.03, 0.04])
        set(gcf, 'Renderer', 'painters')

        if stimType == 2
            title('left Stimulus')
        elseif stimType == 3
            title('right stimulus')
        end
    end
    % xticks([1,2])
    % xticklabels({'left stimuli', 'right stimuli'})
    % 
    % ylabel('Amplitude (peak\DeltaF/F))')
    % set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)

end

regionLRAmps = nan(size(roiList, 2) - 1, 3); %for production of bar with bg plot
for reg = 1:size(regionLRAmps, 1)
    regionLRAmps(reg, :) = [squeeze(nanmean(regionisedAmps(:,[2,3], block, 1 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,[1], block, 1, end)))];
end

%cross-ROI bar - thesis taylored (bg must be last element)
LRNorms = nan(3, size(regionisedAmps,3), size(regionisedAmps,4), size(regionisedAmps,5));
LROffset = [-0.222, 0, 0.222]; %hardcoded for 3 bars/roi
if size(roiList, 2) >= 2 && ( any(ismember(roiList, 'background')) || any(ismember(roiList, 'bg')) )
    % bgIndex = ismember(roiList, 'bg'); %maybe something to account for an unexpected position of bg
    for block = 1:size(regionisedAmps, 3)
        figure
        hold on
        % LRBar = bar([1:size(regionLRAmps, 1)], regionLRAmps);
        for reg = 1:size(regionisedAmps, 5) - 1
            % barX = LRBar(type).XData + LRBar(type).XOffset;
            % errorbar(barX, regionLRAmps(:,type), squeeze(regionisedAmpSEM(:, type, block, 1, end-1)))

            %stats
            LRData = [squeeze(regionisedAmps(:,[2,3],block,1,reg)), squeeze(regionisedAmps(:,1,block,1,end))];
            % LRNorms = [];
            for feat = 1:size(LRData, 2)
                LRNorms(feat, block, 1, reg) = lillietest(LRData(:,feat)); %normality test
            end

            nonNormCheck = 0;
            for col = 1:size(LRData, 2)
                normTest = lillietest(LRData(:,col));
                if ~normTest
                    nonNormCheck = 1;
                end
            end

            if nonNormCheck
                disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' LR'])
                % [~,~,LRStats] = kruskalwallis(LRData, [] , 'off');
                [~,~,LRStats] = friedman(LRData, 1 , 'off');
                % close, close %closing obligatory figures that anova makes
                LRSig = multcompare(LRStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % LRDescription = ['kruskal-wallis test', newline];
                LRDescription = ['Friedman test', newline];
            else
                % [~,~,LRStats] = anova1(LRData, [] , 'off');
                [~,~,LRStats] = ranova(LRData, [] , 'off');
                % close, close %closing obligatory figures that anova makes
                LRSig = multcompare(LRStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % LRDescription = ['one-way anova', newline];
                LRDescription = ['repeated-measures anova', newline];
            end

            % LRDescription = [];
            for group = 1:size(LRSig, 1)
                LRDescription = [LRDescription, [num2str(LRSig(group,1)), 'v', num2str(LRSig(group,2)), ': ', num2str(LRSig(group,6)), newline]];
            end

            % % roiLR = 
            bar([reg], [squeeze(nanmean(regionisedAmps(:,[2,3], block, 1 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,[1], block, 1, end)))]);
            % barX = [currentRegion.XOffset]
            % barX(1), barX(2), barX(3)
            % barX = barX + [currentRegion.XOffset]
            barX = LROffset + reg; %offset hardcoded cause matlab was being irritating
            errorbar(barX, [squeeze(nanmean(regionisedAmps(:,[2,3], block, 1 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,[1], block, 1, end)))], [regionisedAmpSEM(:,[2,3], block, 1, reg), regionisedAmpSEM(:,[1], block, 1, end)], 'LineStyle', 'none', 'Color', 'k')
            if doDots
                scatter(barX, [squeeze(regionisedAmps(:,[2,3], block, 1 , reg)) , squeeze(regionisedAmps(:,[1], block, 1, end))], 10, 'k', 'filled')
            end

            %displaying stats
            figMax = get(gca, 'ylim');
            figMax(1) = [];
            text(barX(end), figMax/2, LRDescription)
        end
        xticks(1:size(regionisedAmps, 5) - 1)
        xticklabels({roiList{1:end-1}})

        ylabel('Amplitude (peak\DeltaF/F)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    end
end

%comparing ROI to background
if size(roiList, 2) == 2 && any(ismember(roiList, 'background'))
    bgIndex = ismember(roiList, 'background'); 
    mainRegionName = overStruct(~bgIndex).regionName;
    for block = 1:size(regionisedAmps, 3)
        figure
        hold on
        bar([1,2], nanmean(regionisedAmps(:,[2,3], block, 1, ~bgIndex), 1)) %amp of left vs right for main ROI
        bar([3], nanmean(regionisedAmps(:,[1], block, 1, bgIndex))) %amp of combined LR for blank
        xticks(1:3)
        xticklabels({['left (', mainRegionName, ')'], ['right (', mainRegionName, ')'], 'background'})
    end
end

% for roi = 1:size()

%lum
%roi based bars
for block = 1:size(regionisedAmps, 3)
    figure
    hold on
    bar( [1:size(regionisedAmps, 5)]', squeeze(nanmean(regionisedAmps(:,1:7, block, 2, :), 1))' ) % number of luminosities hardcoded
    %errorbar
    xticks(1:size(regionisedAmps, 5))
    xticklabels(roiList)
    %spiffy sets
    title(['Cross Fly Lum amps - block ' num2str(block)])
end

%grouping by ROI
for block = 1:size(regionisedAmps, 3)
    figure
    hold on
    for stimType = 1:7
        lumData = squeeze(regionisedAmps(:, stimType, block, 2, :));

        nonNormCheck = 0;
        for col = 1:size(lumData, 2)
            normTest = lillietest(lumData(:,col));
            if ~normTest
                nonNormCheck = 1;
            end
        end

        if nonNormCheck
            disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' lum'])
            % [~,~,lumStats] = kruskalwallis(lumData, [] , 'off');
            [~,~,lumStats] = friedman(lumData, 1 , 'off');
            lumSig = multcompare(lumStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            % lumDescription = ['kruskal-wallis test', newline];
            lumDescription = ['Friedman test', newline];
        else
            % [~,~,lumStats] = anova1(lumData, [] , 'off');
            [~,~,lumStats] = ranova(lumData, [] , 'off');
            lumSig = multcompare(lumStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
            % lumDescription = ['one-way anova', newline];
            lumDescription = ['repeated-measures anova', newline];
        end

        % lumDescription = [];
        for group = 1:size(lumSig, 1)
            lumDescription = [lumDescription, [num2str(lumSig(group,1)), 'v', num2str(lumSig(group,2)), ': ', num2str(lumSig(group,6)), newline]];
        end

        %plotting
        bar(stimType, squeeze(nanmean(regionisedAmps(:, stimType, block, 2, :), 1)) )
        barX = RoiOffset + stimType;
        errorbar(barX, squeeze(nanmean(regionisedAmps(:, stimType, block, 2, :), 1)), squeeze(regionisedAmpSEM(:,stimType, block, 2, :)), 'LineStyle', 'none', 'Color', 'k')
        if doDots
            scatter(barX, squeeze(regionisedAmps(:, stimType, block, 2, :)), 10, 'k', 'filled')
        end

        %displaying stats
        figMax = get(gca, 'ylim');
        figMax(1) = [];
        text(barX(end), figMax/2, lumDescription)
    end
    xticks([1:7])
    % xticklabels(lumTitles)

    ylabel('Amplitude (peak\DeltaF/F))')
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
end

%bg as control - thesis based
% lumOffset = [-0.3500   -0.2500   -0.1500   -0.0500    0.0500    0.1500    0.2500    0.3500]; %hardcoded for 8 bar/roi 
lumOffset = [ -0.3429   -0.2286   -0.1143         0    0.1143    0.2286    0.3429]; %hardcoded for 7 bars
if size(roiList, 2) >= 2 && ( any(ismember(roiList, 'background')) || any(ismember(roiList, 'bg')) )
    for block = 1:size(regionisedAmps, 3)
        figure
        hold on
        for reg = 1:size(regionisedAmps, 5) - 1 %-1 because don't want to plot bg as its own roi
            %stats
            % lumData = [squeeze(regionisedAmps(:,1:7,block,2,reg)), squeeze(regionisedAmps(:,7,block,2,end))];
            lumData = squeeze(regionisedAmps(:,1:7,block,2,reg));

            nonNormCheck = 0;
            for col = 1:size(lumData, 2)
                normTest = lillietest(lumData(:,col));
                if ~normTest
                    nonNormCheck = 1;
                end
            end

            if nonNormCheck
                disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' lum'])
                % [~,~,lumStats] = kruskalwallis(lumData, [] , 'off');
                [~,~,lumStats] = friedman(lumData, 1 , 'off');
                lumSig = multcompare(lumStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % lumDescription = ['kruskal-wallis test', newline];
                lumDescription = ['Friedman test', newline];
            else
                % [~,~,lumStats] = anova1(lumData, [] , 'off');
                [~,~,lumStats] = ranova(lumData, [] , 'off');
                lumSig = multcompare(lumStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % lumDescription = ['one-way anova', newline];
                lumDescription = ['repeated-measures anova', newline];
            end

            % lumDescription = [];
            for group = 1:size(lumSig, 1)
                lumDescription = [lumDescription, [num2str(lumSig(group,1)), 'v', num2str(lumSig(group,2)), ': ', num2str(lumSig(group,6)), newline]];
            end

            % bar(reg, [squeeze(nanmean(regionisedAmps(:,[1:7], block, 2 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,7, block, 2, end)))]); %using bg to highest lum for now
            bar(reg, squeeze(nanmean(regionisedAmps(:,[1:7], block, 2 , reg), 1)));
            barX = lumOffset + reg;
            % errorbar(barX, [squeeze(nanmean(regionisedAmps(:,[1:7], block, 2 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,7, block, 2, end)))],  [regionisedAmpSEM(:,[1:7], block, 2, reg), regionisedAmpSEM(:,7, block, 2, end)], 'LineStyle', 'none', 'Color', 'k') 
            errorbar(barX, squeeze(nanmean(regionisedAmps(:,[1:7], block, 2 , reg), 1)),  squeeze(regionisedAmpSEM(:,[1:7], block, 2, reg)), 'LineStyle', 'none', 'Color', 'k') 
            if doDots
                % scatter(barX, [squeeze(regionisedAmps(:,[1:7], block, 2 , reg)) , squeeze(regionisedAmps(:,7, block, 2, end))], 10, 'k', 'filled')
                scatter(barX, squeeze(regionisedAmps(:,[1:7], block, 2 , reg)), 10, 'k', 'filled')
            end
            % colormap(jet)

            %displaying stats
            figMax = get(gca, 'ylim');
            figMax(1) = [];
            text(barX(end), figMax/2, lumDescription)
        end
        xticks([1:(size(regionisedAmps, 5) -1) ])
        xticklabels({roiList{1:end-1}})

        ylabel('Amplitude (peak\DeltaF/F))')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    end
end

%freq
%roi based bars
for block = 1:size(regionisedAmps, 3)
    figure
    hold on
    bar( [1:size(regionisedAmps, 5)]', squeeze(nanmean(regionisedAmps(:,1:5, block, 3, :), 1))' ) % number of frequencies hardcoded
    %errorbar
    xticks(1:size(regionisedAmps, 5))
    xticklabels(roiList)
    %spiffy sets
    title(['Cross Fly Lum amps - block ' num2str(block)])
end

%ROI comps
for block = 1:size(regionisedAmps, 3)
    figure
    hold on
    for stimType = 2:4
        %stats
        freqData = squeeze(regionisedAmps(:,stimType, block, 3, :));

        nonNormCheck = 0;
        for col = 1:size(freqData, 2)
            normTest = lillietest(freqData(:,col));
            if ~normTest
                nonNormCheck = 1;
            end
        end

        if nonNormCheck
            disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' freq'])
            % [~,~,freqStats] = kruskalwallis(freqData, [] , 'off');
            [~,~,freqStats] = friedman(freqData, 1 , 'off');
            freqSig = multcompare(freqStats, 'Display', 'off');
            % freqDescription = ['kruskal-wallis test', newline];
            freqDescription = ['Friedman test', newline];
        else
            % [~,~,freqStats] = anova1(freqData, [] , 'off');
            [~,~,freqStats] = ranova(freqData, [] , 'off');
            freqSig = multcompare(freqStats, 'Display', 'off');
            % freqDescription = ['one-way anova', newline];
            freqDescription = ['repeated-measures anova', newline];
        end

        % freqDescription = [];
        for group = 1:size(freqSig, 1)
            freqDescription = [freqDescription, [num2str(freqSig(group,1)), 'v', num2str(freqSig(group,2)), ': ', num2str(freqSig(group,6)), newline]];
        end

        %plotting
        bar(stimType - 1, squeeze(nanmean(regionisedAmps(:,stimType, block, 3, :), 1)))
        barX = RoiOffset + (stimType - 1);
        errorbar(barX, squeeze(nanmean(regionisedAmps(:,stimType, block, 3, :), 1)), squeeze(regionisedAmpSEM(:,stimType, block, 3, :)), 'LineStyle', 'none', 'Color', 'k')
        if doDots
            scatter(barX, squeeze(regionisedAmps(:,stimType, block, 3, :)), 10, 'k', 'filled')
        end

        %displaying stats
        figMax = get(gca, 'ylim');
        figMax(1) = [];
        text(barX(end), figMax/2, freqDescription)
    end
    xticks([1:3])
    xticklabels({'3Hz', '6Hz', '12Hz'})

    ylabel('Amplitude (peak\DeltaF/F))')
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
end

%bg as control - thesis based
freqOffset = [-0.2727   -0.0909    0.0909    0.2727]; %hardcoded from tests - for 4 bars/roi
if size(roiList, 2) >= 2 && ( any(ismember(roiList, 'background')) || any(ismember(roiList, 'bg')) )
    for block = 1:size(regionisedAmps, 3)
        figure
        hold on
        for reg = 1:size(regionisedAmps, 5) - 1
            %stats
            freqData = [squeeze(regionisedAmps(:,2:4,block,3,reg)), squeeze(nanmean(regionisedAmps(:,2:4,block,3,end), 2))];

            nonNormCheck = 0;
            for col = 1:size(freqData, 2)
                normTest = lillietest(freqData(:,col));
                if ~normTest
                    nonNormCheck = 1;
                end
            end

            if nonNormCheck
                disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' freq'])
                % [~,~,freqStats] = kruskalwallis(freqData, [] , 'off');
                [~,~,freqStats] = friedman(freqData, 1 , 'off');
                freqSig = multcompare(freqStats, 'Display', 'off');
                % freqDescription = ['kruskal-wallis test', newline];
                freqDescription = ['Friedman test', newline];
            else
                % [~,~,freqStats] = anova1(freqData, [] , 'off');
                [~,~,freqStats] = ranova(freqData, [] , 'off');
                freqSig = multcompare(freqStats, 'Display', 'off');
                % freqDescription = ['one-way anova', newline];
                freqDescription = ['repeated-measures anova', newline];
            end

            % freqDescription = [];
            for group = 1:size(freqSig, 1)
                freqDescription = [freqDescription, [num2str(freqSig(group,1)), 'v', num2str(freqSig(group,2)), ': ', num2str(freqSig(group,6)), newline]];
            end

            bar(reg, [squeeze(nanmean(regionisedAmps(:,[2:4], block, 3 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,[2:4], block, 3, end), [1,2]))])
            barX = freqOffset + reg;
            errorbar(barX, [squeeze(nanmean(regionisedAmps(:,[2:4], block, 3 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,[2:4], block, 3, end), [1,2]))], [regionisedAmpSEM(:,[2:4], block, 3, reg), squeeze(nanmean(regionisedAmpSEM(:,2:4, block, 3, end), 2))], 'LineStyle', 'none', 'Color', 'k')
            if doDots
                scatter(barX, [squeeze(regionisedAmps(:,[2:4], block, 3 , reg)) , squeeze(nanmean(regionisedAmps(:,[2:4], block, 3, end), 2))], 10, 'k', 'filled')
            end

            %displaying stats
            figMax = get(gca, 'ylim');
            figMax(1) = [];
            text(barX(end), figMax/2, freqDescription)
        end
        xticks(1:size(regionisedAmps, 5) - 1)
        xticklabels({roiList{1:end-1}})

        ylabel('Amplitude (peak\DeltaF/F)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
    end
end

%eventually put in an easy access cross block figGen (for induction exps)

%% traces
%LR
regionisedLR = nan(3, size(allFlyLR, 2), size(allFlyLR, 3), size(allFlyLR, 4), size(overStruct, 2)); % [ combined(1)/left(2)/right(3), vol, data(1)/photo(2)/time(3), block, ROI ] 
regionisedLRSEM = nan(3, size(regionisedLR, 2), size(regionisedLR, 4), size(overStruct, 2)); % [ combined(1)/left(2)/right(3), vol, block, ROI ] 

%collection
for roi = 1:size(overStruct, 2)
    %combined
    regionisedLR(1, :, :, :, roi) = nanmean(overStruct(roi).LR, 1);
    regionisedLRSEM(1,:,:,roi) = overStruct(roi).LRSEM;
    %left
    regionisedLR(2, :, :, :, roi) = nanmean(overStruct(roi).lefts, 1);
    regionisedLRSEM(2,:,:,roi) = overStruct(roi).leftSEM;
    %right
    regionisedLR(3, :, :, :, roi) = nanmean(overStruct(roi).rights, 1);
    regionisedLRSEM(3,:,:,roi) = overStruct(roi).rightSEM;
end

%plotting
sideNames = {'combined LR', 'left', 'right'};
for block = 1:size(regionisedLR, 4)
    for side = 1:size(regionisedLR, 1)
        figure
        hold on
        plot(squeeze(regionisedLR(side, :, 3, block, :)), squeeze(regionisedLR(side, :, 1, block, :)), 'LineWidth', 3) % will have to specify colours eventually so errobars match - may be easiest to assume ROI order is same as roiList
        errorbar(squeeze(regionisedLR(side, :, 3, block, :)), squeeze(regionisedLR(side, :, 1, block, :)), squeeze(regionisedLRSEM(side, :, block, :)), 'LineStyle', 'none')
        colours = get(gca, 'ColorOrder');
        
        set(gca, 'ylim', [-0.02, 0.015])

        set(gca, 'ColorOrder', colours(1:size(roiList, 2), :))
        yLims = get(gca, 'ylim');
        fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
        legend(roiList)
        title(['block ', num2str(block), ' ', sideNames{side}, ' stimulus'])

        ylabel('mean fluorescence (\DeltaF/F)')
        xlabel('time (s)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        set(gcf, 'Renderer', 'painters')
    end
end


%lum
regionisedLum = nan(size(allFlyLum, 4), size(allFlyLum, 2), size(allFlyLum, 3), size(allFlyLum, 5), size(overStruct, 2)); % [ lum, vol, data(1)/phot(2)/time(3), block, roi ]
regionisedLumSEM = nan(size(regionisedLum, [1,2,4,5])); %[ lum, vol, block, roi ]

%collecting
for roi = 1:size(overStruct, 2)
    for lum = 1:size(regionisedLum, 1)
        regionisedLum(lum, :, :, :, roi) = nanmean(overStruct(roi).lum(:,:,:,lum,:), 1);
        regionisedLumSEM(lum, :, :, roi) = overStruct(roi).lumSEM(:, lum, :);
    end
end

%plotting
lumTitles = {'Luminosity 1', 'Luminosity 2', 'Luminosity 3', 'Luminosity 4', 'Luminosity 5', 'Luminosity 6', 'Luminosity 7'}; %change to lux eventually - maybe just do in illustrator
for block = 1:size(regionisedLum, 4)
    for lum = 1:size(regionisedLum, 1)
        figure
        hold on
        plot(squeeze(regionisedLum(lum, :, 3, block, :)), squeeze(regionisedLum(lum, :, 1, block, :)), 'LineWidth', 3)
        errorbar(squeeze(regionisedLum(lum, :, 3, block, :)), squeeze(regionisedLum(lum, :, 1, block, :)), squeeze(regionisedLumSEM(lum, :, block, :)), 'LineStyle', 'none')
        colours = get(gca, 'ColorOrder');

        set(gca, 'ylim', [ -0.0500    0.0400 ]) %hardcoded form largest lum

        set(gca, 'ColorOrder', colours(1:size(roiList, 2), :))
        yLims = get(gca, 'ylim');
        fill([0 0 1.089 1.089], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
        legend(roiList, 'Location', 'southeast')
        title(['block ', num2str(block), lumTitles{lum}])

        ylabel('mean fluorescence (\DeltaF/F)')
        xlabel('time (s)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        set(gcf, 'Renderer', 'painters')
    end
end

%freq
regionisedFreq = nan(size(allFlyFreq, 4), size(allFlyFreq, 2), size(allFlyFreq, 3), size(allFlyFreq, 5), size(overStruct, 2)); % [ freq, vol, data(1)/photo(2)/time(3), block, roi ]
regionisedFreqSEM = nan(size(regionisedFreq, [1,2,4,5])); % [ freq, vol, block, roi ]

%collection 
for roi = 1:size(overStruct, 2)
    for freq = 1:size(regionisedFreq, 1)
        regionisedFreq(freq, :, :, :, roi) = nanmean(overStruct(roi).freq(:, :, :, freq, :), 1);
        regionisedFreqSEM(freq, :, :, roi) = overStruct(roi).freqSEM(:, freq, :); 
    end
end

%plotting
for block = 1:size(regionisedFreq, 4)
    for freq = 1:size(regionisedFreq, 1)
        figure
        hold on
        plot(squeeze(regionisedFreq(freq, :, 3, block, :)), squeeze(regionisedFreq(freq, :, 1, block, :)), 'LineWidth', 3)
        errorbar(squeeze(regionisedFreq(freq, :, 3, block, :)), squeeze(regionisedFreq(freq, :, 1, block, :)), squeeze(regionisedFreqSEM(freq, :, block, :)), 'LineStyle', 'none')
        colours = get(gca, 'ColorOrder');
        set(gca, 'ColorOrder', colours(1:size(roiList, 2), :), 'ylim', [ -0.0400    0.0300 ])

        yLims = get(gca, 'ylim');
        fill([0, 0, 1.271, 1.271], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
        % fill([0, 0, freqStimEstimate(freq), freqStimEstimate(freq)], [yLims(1), yLims(2), yLims(2), yLims(1)], 'b', 'FaceAlpha', 0.3)
        legend(roiList)
        title(['block ', num2str(block), ' ', num2str(matSave.batFrequencies(freq)), ' Hz'])

        ylabel('mean fluorescence (\DeltaF/F)')
        xlabel('time (s)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12)
        set(gcf, 'Renderer', 'painters')
    end
end

%% behavioural state analyses (crossBlock for all ROIs)
% ^^ might have made more sense to have this earlier (logically), but I currently need and understand regionisedAmps most

if size(roiList, 2) < 2 || ~( any(ismember(roiList, 'background')) || any(ismember(roiList, 'bg')) )
    disp('ROI listing issue stopping multiBlock analysis: this section of the code requires an ROI and background to be listed (background must be last)')
elseif size(regionisedAmps, 3) == 1
    disp('Not enough blocks given for multiBlock amplitude analysis :(')
else
    disp('crossBlock analysis reached! yay')
    %crossBlock stats
    %stats
    %% 
    tic
    LRBlockStats = cell(1,3, size(regionisedAmps, 3)-1); %crossBlock descriptions for [left, right, bg] of each ROI (dimension 3)
    lumBlockStats = cell(1,7, size(regionisedAmps, 3)-1); %crossBlock descritpions for [lum1-7, bg] of each ROI (dimension 3)
    freqBlockStats = cell(1,3, size(regionisedAmps, 5)-1); %crossBlock descriptions for [3Hz, 6Hz, 12Hz, bg] of each ROI (dimension 3)
    conditionNames = {'LR', 'luminosity', 'frequency'}; % for display purposes
    for reg = 1:size(regionisedAmps, 5)
        for condition = 1:size(regionisedAmps, 4)
            featCount = 1; %used later for saving to blockStat variables
            if condition == 1
                featSize = 1:3;
                bgFeat = 1;
            elseif condition == 2
                featSize = 1:7;
                bgFeat = 7;
            elseif condition == 3
                featSize = 2:4;
                bgFeat = 2:4;
            end

            for feature = featSize %1:size(regionisedAmps, 3)
                % featCount = 1; %used later for saving to blockStat variables
                currentData = [];
                for block = 1:size(regionisedAmps, 3)
                    if reg == size(regionisedAmps, 5) % if background
                        currentData = [currentData, [squeeze(nanmean(regionisedAmps(:,bgFeat, block, condition, reg), 2)) ]];
                    else
                        currentData = [currentData, [squeeze(regionisedAmps(:,feature, block, condition, reg)) ]];
                    end
                end
                % currentData = [currentData, squeeze(nanmean(regionisedAmps(:,bgFeat, block, condition, end), 2))]; %adding background as a feature

                %testing normality
                nonNormCheck = 0;
                for col = 1:size(LRData, 2)
                    normTest = lillietest(currentData(:,col));
                    if ~normTest
                        nonNormCheck = 1;
                    end
                end

                %performing statistcs based on normality
                if nonNormCheck
                    disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' ', num2str(conditionNames{condition})])
                    % [~,~,currentStats] = kruskalwallis(currentData, [] , 'off');
                    [~,~,currentStats] = friedman(currentData, 1 , 'off');
                    currentSig = multcompare(currentStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                    % currentDescription = ['kruskal-wallis test', newline];
                    currentDescription = ['friedman test', newline];
                else
                    disp(['all groups normally distributed, one-way ANOVA used - ', num2str(roiList{reg}), ' ', num2str(conditionNames{condition})])
                    % [~,~,currentStats] = anova1(currentData, [] , 'off');
                    [~,~,currentStats] = ranova(currentData, [] , 'off');
                    currentSig = multcompare(currentStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                    % currentDescription = ['one-way anova', newline];
                    currentDescription = ['repeated-measures anova', newline];
                end

                for group = 1:size(currentSig, 1)
                    currentDescription = [currentDescription, [num2str(currentSig(group,1)), 'v', num2str(currentSig(group,2)), ': ', num2str(currentSig(group,6)), newline]];
                end

                %storing statDescriptions
                if condition == 1
                    LRBlockStats{1, featCount, reg} = currentDescription;
                elseif condition == 2
                    lumBlockStats{1, featCount, reg} = currentDescription;
                elseif condition == 3
                    freqBlockStats{1, featCount, reg} = currentDescription;
                end
                featCount = featCount + 1;
            end
        end
    end
    toc

    %LR
    % yLimsBehav = [-0.0200 , 0.0400]; %for batCtrls
    yLimsBehav = [-0.015, 0.03]; %atr
    for reg = 1:size(regionisedAmps, 5) - 1
        figure
        hold on
        for block = 1:size(regionisedAmps, 3)
            %stats
            LRData = [squeeze(regionisedAmps(:,[2,3],block,1,reg)), squeeze(regionisedAmps(:,1,block,1,end))];

            %{
            nonNormCheck = 0;
            for col = 1:size(LRData, 2)
                normTest = lillietest(LRData(:,col));
                if ~normTest
                    nonNormCheck = 1;
                end
            end

            
            if nonNormCheck
                disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' LR'])
                % [~,~,LRStats] = kruskalwallis(LRData, [] , 'off');
                [~,~,LRStats] = friedman(LRData, 1 , 'off');
                % close, close %closing obligatory figures that anova makes
                LRSig = multcompare(LRStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % LRDescription = ['kruskal-wallis test', newline];
                LRDescription = ['Friedman test', newline];
            else
                [~,~,LRStats] = anova1(LRData, [] , 'off');
                LRSig = multcompare(LRStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % LRDescription = ['one-way anova', newline];
                LRDescription = ['repeated-measures anova', newline];
            end

            % LRDescription = [];
            for group = 1:size(LRSig, 1)
                LRDescription = [LRDescription, [num2str(LRSig(group,1)), 'v', num2str(LRSig(group,2)), ': ', num2str(LRSig(group,6)), newline]];
            end
            %}


            bar(block, [squeeze(nanmean(regionisedAmps(:,[2,3], block, 1, reg), 1)), squeeze(nanmean(regionisedAmps(:,1,block,1,end), 1))])
            barX = LROffset + block;
            errorbar(barX, [squeeze(nanmean(regionisedAmps(:,[2,3], block, 1, reg), 1)), squeeze(nanmean(regionisedAmps(:,1,block,1,end), 1))], [regionisedAmpSEM(:,[2,3], block, 1, reg), regionisedAmpSEM(:,[1], block,1,end)], 'LineStyle', 'none', 'Color', 'k')
            if doDots
                scatter(barX, [squeeze(regionisedAmps(:,[2,3], block, 1, reg)), squeeze(regionisedAmps(:,1,block,1,end))], 10, 'k', 'filled')
            end
            %displaying stats
            %{
            figMax = get(gca, 'ylim');
            figMax(1) = [];
            text(barX(end), figMax/2, LRDescription)
            %}
        end
        xticks(1:size(regionisedAmps, 3))
        xticklabels({'Baseline', 'Red-light', 'Recovery'})
        title([roiList{reg}, ' cross block LR amps vs bg'])

        ylabel('Amplitude (peak\DeltaF/F)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12) % 'ylim', yLimsBehav)
    end

    %lum
    for reg = 1:size(regionisedAmps, 5) -1
        figure
        hold on
        for block = 1:size(regionisedAmps, 3)
            %stats
            lumData = [squeeze(regionisedAmps(:,1:7,block,2,reg)), squeeze(regionisedAmps(:,7,block,2,end))];

            %{
            nonNormCheck = 0;
            for col = 1:size(lumData, 2)
                normTest = lillietest(lumData(:,col));
                if ~normTest
                    nonNormCheck = 1;
                end
            end

            if nonNormCheck
                disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' lum'])
                % [~,~,lumStats] = kruskalwallis(lumData, [] , 'off');
                [~,~,lumStats] = friedman(lumData, 1 , 'off');
                lumSig = multcompare(lumStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % lumDescription = ['kruskal-wallis test', newline];
                lumDescription = ['Friedman test', newline];
            else
                % [~,~,lumStats] = anova1(lumData, [] , 'off');
                [~,~,lumStats] = ranova(lumData, [] , 'off');
                lumSig = multcompare(lumStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
                % lumDescription = ['one-way anova', newline];
                lumDescription = ['repeated-measures anova', newline];
            end

            % lumDescription = [];
            for group = 1:size(lumSig, 1)
                lumDescription = [lumDescription, [num2str(lumSig(group,1)), 'v', num2str(lumSig(group,2)), ': ', num2str(lumSig(group,6)), newline]];
            end
            %}

            bar(block, [squeeze(nanmean(regionisedAmps(:,[1:7], block, 2 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,7, block, 2, end)))]); %using bg to highest lum for now
            barX = lumOffset + block;
            errorbar(barX, [squeeze(nanmean(regionisedAmps(:,[1:7], block, 2 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,7, block, 2, end)))],  [regionisedAmpSEM(:,[1:7], block, 2, reg), regionisedAmpSEM(:,7, block, 2, end)], 'LineStyle', 'none', 'Color', 'k')
            if doDots
                scatter(barX, [squeeze(regionisedAmps(:,[1:7], block, 2 , reg)) , squeeze(regionisedAmps(:,7, block, 2, end))], 10, 'k', 'filled')
            end
            %displaying stats
            figMax = get(gca, 'ylim');
            figMax(1) = [];
            text(barX(end), figMax/2, lumDescription)
        end
        xticks(1:size(regionisedAmps, 3))
        xticklabels({'Baseline', 'Red-light', 'Recovery'})
        title([roiList{reg}, ' cross block Lum amps vs bg'])

        ylabel('Amplitude (peak\DeltaF/F)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12) % 'ylim', yLimsBehav)
    end

    %freq
    for reg = 1:size(regionisedAmps, 5) -1
        figure
        hold on
        for block = 1:size(regionisedAmps, 3)
            %stats
            freqData = [squeeze(regionisedAmps(:,2:4,block,3,reg)), squeeze(nanmean(regionisedAmps(:,2:4,block,3,end), 2))];

            %{
            nonNormCheck = 0;
            for col = 1:size(freqData, 2)
                normTest = lillietest(freqData(:,col));
                if ~normTest
                    nonNormCheck = 1;
                end
            end

            if nonNormCheck
                disp(['a group was not normally distributed, kruskal-wallis test being used for stats - ', num2str(roiList{reg}), ' freq'])
                % [~,~,freqStats] = kruskalwallis(freqData, [] , 'off');
                [~,~,freqStats] = friedman(freqData, 1 , 'off');
                freqSig = multcompare(freqStats, 'Display', 'off');
                % freqDescription = ['kruskal-wallis test', newline];
                freqDescription = ['Friedman test', newline];
            else
                % [~,~,freqStats] = anova1(freqData, [] , 'off');
                [~,~,freqStats] = ranova(freqData, [] , 'off');
                freqSig = multcompare(freqStats, 'Display', 'off');
                % freqDescription = ['one-way anova', newline];
                freqDescription = ['repeated-measures anova', newline];
            end

            % freqDescription = [];
            for group = 1:size(freqSig, 1)
                freqDescription = [freqDescription, [num2str(freqSig(group,1)), 'v', num2str(freqSig(group,2)), ': ', num2str(freqSig(group,6)), newline]];
            end
            %}
            
            bar(block, [squeeze(nanmean(regionisedAmps(:,[2:4], block, 3 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,[2:4], block, 3, end), [1,2]))])
            barX = freqOffset + block;
            errorbar(barX, [squeeze(nanmean(regionisedAmps(:,[2:4], block, 3 , reg), 1)) , squeeze(nanmean(regionisedAmps(:,[2:4], block, 3, end), [1,2]))], [regionisedAmpSEM(:,[2:4], block, 3, reg), squeeze(nanmean(regionisedAmpSEM(:,2:4, block, 3, end), 2))], 'LineStyle', 'none', 'Color', 'k')
            if doDots
                scatter(barX, [squeeze(regionisedAmps(:,[2:4], block, 3 , reg)) , squeeze(nanmean(regionisedAmps(:,[2:4], block, 3, end), 2))], 10, 'k', 'filled')
            end
            %displaying stats
            %{
            figMax = get(gca, 'ylim');
            figMax(1) = [];
            text(barX(end), figMax/2, freqDescription)
            %}
            clear freqDescription freqStats
        end
        xticks(1:size(regionisedAmps, 3))
        xticklabels({'Baseline', 'Red-light', 'Recovery'})
        title([roiList{reg}, ' cross block freq amps vs bg'])

        ylabel('Amplitude (peak\DeltaF/F)')
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12) % 'ylim', yLimsBehav)
    end

end

%% testing

% %code to generate array of lumVol sizes [luminosity, repitition]
% lumSize = [];
% for l = 1:7
%     currentSizes = [];
%     for n = 1:140
%         currentSizes = [currentSizes, sum(~isnan(lumVols(n, :, 2, l)'))];
%     end
%     lumSize = [lumSize;currentSizes];
% end

if singleTesting

    % %creating time-scale with full photoDiode # of points
    % %WARNING: takes a very long time to run! will be commented out indefinitely unless wanted
    % fullPhotFitted = [];
    % for t = 1:size(volTimes, 2) - 1
    %     fullPhotFitted = [fullPhotFitted, linspace(volTimes(1, t), volTimes(1, t+1), volTimes(3, t+1) - volTimes(3, t) )];
    % end

    % % plot(fullPhotFitted-volTimes(1,1), smoothPhot(volTimes(3,1):volTimes(3, end)-1))

    %timings currently incorrect for below (I think they work off iSequence)
    % Code for single block testing
    %blanks comparison plots and calcs
    blankPhoto = smoothPhot([volTimes(3, find(volTimes(1, :) == batStruct.blocks(1).blankRawTimes(1))) : volTimes(3, find(volTimes(1, :) == batStruct.blocks(1).blankRawTimes(end)))]);
    optoPhoto5s = [];
    for x = 1:size(batStruct.blocks(1).optoRawTimes, 2)
        photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).optoRawTimes(x) ) );
        photo = FLIES.BLOCKS(b).photData(photoIndex);
        optoPhoto5s = [optoPhoto5s, photo];
    end
    fullOptoPhoto = [];
    for x = 1:size(batStruct.blocks(1).fullOptoRawTimes, 2)
        photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).fullOptoRawTimes(x) ) );
        photo = FLIES.BLOCKS(b).photData(photoIndex);
        fullOptoPhoto = [fullOptoPhoto, photo];
    end
    LRBlankPhoto = [];
    for x = 1:size(batStruct.blocks(1).LRBlankRawTimes, 2)
        photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).LRBlankRawTimes(x) ) );
        photo = FLIES.BLOCKS(b).photData(photoIndex);
        LRBlankPhoto = [LRBlankPhoto, photo];
    end
    freqBlankPhoto = [];
    for x = 1:size(batStruct.blocks(1).freqBlankRawTimes, 2)
        photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).freqBlankRawTimes(x) ) );
        photo = FLIES.BLOCKS(b).photData(photoIndex);
        freqBlankPhoto = [freqBlankPhoto, photo];
    end
    lumBlankPhoto = timedPhoto(1, [find(volTimes(1, :) == batStruct.blocks(1).lumBlankRawTimes(1)) : find(volTimes(1, :) == batStruct.blocks(1).lumBlankRawTimes(end))]);
    % blankPhoto = [];
    % for x = 1:size(batStruct.blocks(1).blankRawTimes, 2)
    %     photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).blankRawTimes(x) ) );
    %     photo = FLIES.BLOCKS(b).photData(photoIndex);
    %     blankPhoto = [blankPhoto, photo];
    % end
    % optoPhoto5s = [];
    % for x = 1:size(batStruct.blocks(1).optoRawTimes, 2)
    %     photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).optoRawTimes(x) ) );
    %     photo = FLIES.BLOCKS(b).photData(photoIndex);
    %     optoPhoto5s = [optoPhoto5s, photo];
    % end
    % fullOptoPhoto = [];
    % for x = 1:size(batStruct.blocks(1).fullOptoRawTimes, 2)
    %     photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).fullOptoRawTimes(x) ) );
    %     photo = FLIES.BLOCKS(b).photData(photoIndex);
    %     fullOptoPhoto = [fullOptoPhoto, photo];
    % end
    % LRBlankPhoto = [];
    % for x = 1:size(batStruct.blocks(1).LRBlankRawTimes, 2)
    %     photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).LRBlankRawTimes(x) ) );
    %     photo = FLIES.BLOCKS(b).photData(photoIndex);
    %     LRBlankPhoto = [LRBlankPhoto, photo];
    % end
    % freqBlankPhoto = [];
    % for x = 1:size(batStruct.blocks(1).freqBlankRawTimes, 2)
    %     photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).freqBlankRawTimes(x) ) );
    %     photo = FLIES.BLOCKS(b).photData(photoIndex);
    %     freqBlankPhoto = [freqBlankPhoto, photo];
    % end
    % lumBlankPhoto = [];
    % for x = 1:size(batStruct.blocks(1).lumBlankRawTimes, 2)
    %     photoIndex = volTimes(3, find( volTimes(1, :) == batStruct.blocks(1).lumBlankRawTimes(x) ) );
    %     photo = FLIES.BLOCKS(b).photData(photoIndex);
    %     lumBlankPhoto = [lumBlankPhoto, photo];
    % end
    %add photodiode to opto trace
    plot(batStruct.blocks(1).optoRawTimes, optoPhoto5s*200 + 700, 'Color', 'k')
    plot(batStruct.blocks(1).blankRawTimes, blankPhoto*200 + 700, 'Color', 'k')

    % clunky - won't work unless only single block being analysed
    figure
    subplot(2, 2, [1])
    plot(batStruct.blocks(1).blankRawTimes, blankPhoto)
    title('opto')
    subplot(2,2,[2])
    plot(batStruct.blocks(1).LRBlankRawTimes, LRPhoto)
    title('LR')
    subplot(2,2,[3])
    plot(batStruct.blocks(1).lumBlankRawTimes, lumBlankPhoto)
    title('lum')
    subplot(2,2,[4])
    plot(batStruct.blocks(1).freqBlankRawTimes, freqBlankPhoto)
    title('freq')
    sgtitle('Blank photodiodes')

    %same as above (clunky) but plots full smooth photodiode with blanks
    figure
    hold on
    plot(volTimes(1,:), timedPhoto)
    fill([batStruct.blocks.blankRawTimes(1), batStruct.blocks.blankRawTimes(1), batStruct.blocks.blankRawTimes(end), batStruct.blocks.blankRawTimes(end)], [0.15, 0.3, 0.3, 0.15] ,'g', 'FaceAlpha', 0.25 )
    fill([batStruct.blocks.LRBlankRawTimes(1), batStruct.blocks.LRBlankRawTimes(1), batStruct.blocks.LRBlankRawTimes(end), batStruct.blocks.LRBlankRawTimes(end)], [0.15, 0.3, 0.3, 0.15] ,'g', 'FaceAlpha', 0.25 )
    fill([batStruct.blocks.lumBlankRawTimes(1), batStruct.blocks.lumBlankRawTimes(1), batStruct.blocks.lumBlankRawTimes(end), batStruct.blocks.lumBlankRawTimes(end)], [0.15, 0.3, 0.3, 0.15] ,'g', 'FaceAlpha', 0.25 )
    fill([batStruct.blocks.freqBlankRawTimes(1), batStruct.blocks.freqBlankRawTimes(1), batStruct.blocks.freqBlankRawTimes(end), batStruct.blocks.freqBlankRawTimes(end)], [0.15, 0.3, 0.3, 0.15] ,'g', 'FaceAlpha', 0.25 )
    title('raw photodiode with blank locations in green')
    xlabel('seconds')

    %plot with full photodiode dataset
    figure
    hold on
    plot(FLIES.BLOCKS.photData)
    fill([volTimes(3, volTimes(1, :) == LR(1)), volTimes(3, volTimes(1, :) == LR(1)), volTimes(3, volTimes(1, :) == LR(end)), volTimes(3, volTimes(1, :) == LR(end))], [0.15, 0.4, 0.4, 0.15], 'g', 'FaceAlpha', 0.25)
    fill([volTimes(3, volTimes(1, :) == lum(1)), volTimes(3, volTimes(1, :) == lum(1)), volTimes(3, volTimes(1, :) == lum(end)), volTimes(3, volTimes(1, :) == lum(end))], [0.15, 0.4, 0.4, 0.15], 'g', 'FaceAlpha', 0.25)
    fill([volTimes(3, volTimes(1, :) == freq(1)), volTimes(3, volTimes(1, :) == freq(1)), volTimes(3, volTimes(1, :) == freq(end)), volTimes(3, volTimes(1, :) == freq(end))], [0.15, 0.4, 0.4, 0.15], 'g', 'FaceAlpha', 0.25)
    title('raw photodiode with blank locations in green')
    xlabel('photodiode index')

    %could be put into loop but would take long time to run
    %code to plot iTimings from BT
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

    %blank photo - from allIVols
    %could be added to real code eventually

    %4x4 photodiode for each blank
    figure
    for blank = 1:size(blankVols, 1) - 1
        subplot(2,2,blank)
        plot(blankVols(blank, :, 2))
        title(num2str(blank))
    end

    figure
    for blank = 1:size(blankVols, 1) - 1
        subplot(2,4,blank)
        plot(blankVols(blank, :, 2))
        title(num2str(blank))

        subplot(2, 4, blank + 4)
        plot(blankVols(blank, :, 1)')
    end


    %code for adding photodiode to same plot as traces - will be added universally, but added here for reference
    plot(normalize(LRVols(:,:,2)', 'range', [0, max(range(LRVols(:,:,1)))/2]) + max(LRVols(:,:,1)) + max(range(LRVols(:,:,1)))/4 ) % raw photos
    plot(nanmean((normalize(LRVols(:,:,2)', 'range', [0, max(range(LRVols(:,:,1)))/2]) + max(LRVols(:,:,1)) + max(range(LRVols(:,:,1)))/4)'), 'LineWidth', 3, 'Color', 'b' ) %avg photo

    %plot scatter of amplitude positions in time (peaks=2, troughs=1)
    figure
    hold on
    for a = 1:3
        for f = 1:5
            scatter(LRAmpPositions(f,:,a,1), 2*ones(size(LRAmpPositions, 2)),'.', 'Color', 'b')
            scatter(LRAmpPositions(f,:,a,2), ones(size(LRAmpPositions, 2)),'.', 'Color', 'r')
        end
    end
    ylim([0,3])

    %bar plots of amplitudes (LR - crossfly)
    figure
    bar([1,2,3], [[nanmean(LRAmps(:,1))], [nanmean(LRAmps(:,2))], [nanmean(LRAmps(:,3))]])
    %bar plots of amplitudes (lum-individualFly)
    for fly = 1:size(allFlyLum, 1)
        for block = 1:size(allFlyLum, 5)
            figure
            bar([1:size(allFlyLum, 4)], [lumAmps(fly,:,block)])

            title(['fly ', num2str(chosenFlies(fly)), ' block ', num2str(block)])
        end
    end
    %barplot (lum-crossFly)
    lumAmpSEM = squeeze(nanstd(lumAmps, [], 1)/sqrt(size(allFlyLum, 1)));
    for block = 1:size(allFlyLum, 5)
        figure
        hold on
        bar(1:size(allFlyLum, 4), nanmean(lumAmps(:,:,block), 1))
        errorbar(1:size(allFlyLum, 4), nanmean(lumAmps(:,:,block), 1), lumAmpSEM(:,block), '.')
    end
    %barplot (freq-crossFly)
    freqAmpSEM = squeeze(nanstd(freqAmps, [], 1)/sqrt(size(allFlyFreq, 1)));
    for block = 1:size(allFlyFreq, 5)
        figure
        hold on
        bar(1:size(allFlyFreq, 4), nanmean(freqAmps(:,:,block), 1))
        errorbar(1:size(allFlyFreq, 4), nanmean(freqAmps(:,:,block), 1), freqAmpSEM(:,block), '.')
    end

    %play videos
    figure
    testFly = batStruct.blocks(1);
    volNum = find(isnan(nanmean(testFly.LRFrames, [1,2])), 1, 'first'); %inlcudes single nan for ease of viewing
    while true
        for vol = 1:volNum
            imagesc(testFly.LRFrames(:,:,vol))
            title(['fly ', num2str(testFly.fly), ' vol ', num2str(vol)])
            drawnow
            pause(0.1)
        end
    end

    %saving crossBlock stats - put here because must be done manually
    fileName = 'batCTRL';
    statDirectory = 'C:\Users\s4749814\2p\2p\RESULTS\DataSets';
    save([statDirectory, filesep, fileName, '.mat'], 'LRBlockStats', 'lumBlockStats', 'freqBlockStats', 'regionisedAmps')


end
durationTime = seconds(datetime('now') - startTime)

%finis