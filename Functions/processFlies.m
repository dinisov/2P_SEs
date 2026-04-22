function processFlies(flyRecord, chosenFlies, chosenBlocks, chosenZ, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, ...
    separateByState, doRolling, recordPath)
%UNTITLED Summary of this function goes here
disp([char(10),'-------------------------------------'])
%% collate, reduce, filter and concatenate pre-aligned data

%FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory,...
%    'alternateUseCase', 0, 'reqZ', chosenZ); %Moved Matt additions to options
FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory,...
    'alternateUseCase', 0, 'reqZ', chosenZ,...
    'dynamicSG',1,'dynamicSGTimeWidth',10,...
    'doRolling',doRolling); %Moved Matt additions to options
        %Note: Dynamic SG may not function correctly for legacy 512x512 data

%% update records
recordUpdater(FLIES,flyRecord,recordPath,1,'analysisState',0)

%% Interrupt flow for rolling datasets
%try
    %[FLIES] = syncMaster( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', 0, 'rollingAnalysis', -1,...
    %    'disregardRollingDesign', 0, 'overwriteShortcut', 1, 'unsiphonedSEs', 0, 'disregardBattery', 1);
    [FLIES] = syncMaster( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', 0, 'rollingAnalysis', -1,...
        'disregardRollingDesign', 0, 'overwriteShortcut', 0, 'unsiphonedSEs', 0, 'disregardBattery', 1,...
        'outputDirectory',outputDirectory,...
        'shiftImTime',0);
%catch
%    ['## syncMaster error; Continuing with next fly ##']
%    return %Skips rest of blocks for this fly
%end

%And check for accidental battery inclusion
if isfield( FLIES.BLOCKS, 'stimulus') && ~isempty( strfind( [FLIES.BLOCKS.stimulus], 'battery' ) )
    ['## Alert: Battery blocks erroneously(?) included in analysis; Aborting ##']
    %crash = yes
    return
end

%% analyse SEs
% separates images according to preceding sequence of stimuli and
% calculates mean images as a function of the sequence

R = analyse2P(FLIES, chosenFlies, outputDirectory, groupedBlocks);

%% make movies of transients as differences to mean
if analysisToggle(1)
    transientMovies(R, chosenFlies, outputDirectory);
end

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS (per volume/time and collapsed across time)
if analysisToggle(2)
    R = componentFits2P(R, groupedBlocks);
end
%% make movies of fits over time
if analysisToggle(3)
    fitMovies(R, outputDirectory, gridSize, chosenFlies);
end
%% calculate mass t-tests
if analysisToggle(4)
    R = ttests2P(R, groupedBlocks);
end
%% RRRR-RRRA and AAAA-AAAR (collapsed and over time videos)
if analysisToggle(5)
    patternPlots(R, chosenFlies, outputDirectory);
end
%% L vs R analysis (t-tests, L-R, L and R, L and R movies)
if analysisToggle(6)
    analyseLvsR(R, chosenFlies, outputDirectory);
end

%% PCA analysis
if analysisToggle(7)
    %flyPCA2D(R, chosenFlies, outputDirectory);
    flyPCA2D(R, chosenFlies, outputDirectory, 'showModel', 0);
end

%% global transient analysis
if analysisToggle(8)
    globalTransients(R, chosenFlies, outputDirectory);
end

%% plotting
if analysisToggle(7)
    disp('Plotting stuff');
    tic;
    % plot results per block
    for fly = 1:length(FLIES)
        subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))]);
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        plotFly(R(fly), groupedBlocks, subDirectory,'off');
    end
    toc;
end

%% update records again
recordUpdater(FLIES,flyRecord,recordPath,1,'analysisState',1)

end

