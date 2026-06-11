function processFlies(flyRecord, chosenFlies, chosenBlocks, chosenZ, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, ...
options)

%Wrapper for the primary 2P SEs analyses

arguments
    flyRecord table
    chosenFlies double
    chosenBlocks cell
    chosenZ cell
    gridSize double
    dataDirectory char
    sequenceDirectory char
    outputDirectory char
    analysisToggle double
    groupedBlocks double
    options.separateByState double = 0
    options.doRolling double = 0
    options.recordPath char
    options.syncManUnsiphonedSEs double = 0
    options.syncManOverwriteShortcut double = 0
    options.syncManDoVid double = 0
end

separateByState = options.separateByState;
doRolling = options.doRolling;
recordPath = options.recordPath;

overwriteShortcut = options.syncManOverwriteShortcut;
unsiphonedSEs = options.syncManUnsiphonedSEs;
doVid = options.syncManDoVid;

disp([char(10),'-------------------------------------'])
%% collate, reduce, filter and concatenate pre-aligned data

%FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory,...
%    'alternateUseCase', 0, 'reqZ', chosenZ); %Moved Matt additions to options
FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory,...
    'alternateUseCase', 0, 'reqZ', chosenZ,...
    'dynamicSG',1,'dynamicSGTimeWidth',10,...
    'doRolling',doRolling,...
    'separateByState',separateByState); %Moved Matt additions to options
        %Note: Dynamic SG may not function correctly for legacy 512x512 data

%% update records
if ~isempty(recordPath) && ~unsiphonedSEs
    recordUpdater(FLIES,flyRecord,recordPath,1,'analysisState',0.5)
else
    disp(['(Skipping record updating)'])
end

%% Interrupt flow for rolling datasets
try
    %[FLIES] = syncMaster( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', 0, 'rollingAnalysis', -1,...
    %    'disregardRollingDesign', 0, 'overwriteShortcut', 1, 'unsiphonedSEs', 0, 'disregardBattery', 1);
    [FLIES] = syncMaster( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', doVid, 'rollingAnalysis', -1,...
        'disregardRollingDesign', 0, 'overwriteShortcut', overwriteShortcut, 'unsiphonedSEs', unsiphonedSEs, 'disregardBattery', 1,...
        'outputDirectory',outputDirectory,...
            'photStorageMode',2);
        %'sequenceObliteration',{[0,0,0,0,0],[1,1,1,1,1]});
        %'shiftImTime',0);
        %'forceNoIterator',0);
catch
    ['## syncMaster error; Continuing with next fly ##']
    return %Skips rest of blocks for this fly
end

%Some QAs/checks
if unsiphonedSEs
    disp(['Unsiphoned SEs requested; Ceasing (Dinis) processing'])
    return
end
if isempty(FLIES.BLOCKS) %Might crash under normal circumstances?
    ['-# No apparent data supplied and/or made it through syncMaster/etc #-']
    return
end
%And check for accidental battery inclusion
if isfield( FLIES.BLOCKS, 'stimulus') && ~isempty( strfind( [FLIES.BLOCKS.stimulus], 'battery' ) )
    ['## Alert: Battery blocks erroneously(?) included in analysis; Aborting ##']
    %crash = yes
    return
end

%% analyse SEs
% separates images according to preceding sequence of stimuli and
% calculates mean images as a function of the sequence

R = analyse2P(FLIES, chosenFlies, outputDirectory, groupedBlocks, 1); %New argument: saveFull

%% make movies of transients as differences to mean
if analysisToggle(1)
    transientMovies(R, chosenFlies, outputDirectory, 'normalisation', 'pixel');
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
%recordUpdater(FLIES,flyRecord,recordPath,1,'analysisState',1)
if ~isempty(recordPath) && ~unsiphonedSEs
    recordUpdater(FLIES,flyRecord,recordPath,1,'analysisState',1)
else
    disp(['(Skipping record updating)'])
end

%disp(['subFin'])

end

