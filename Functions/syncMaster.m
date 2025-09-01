function FLIES = syncMaster(FLIES, flyRecord, options)

%Mk ???
%Mk 6 - Support for battery, and blanks during bendy block design
%Mk 7 - Better support for DAQ data, highspeed recording support
%Mk 8 - Photodiode (initial) support
%Mk 8.5 - Better photodiode support, minor graphs for SEs data
%Mk 9 - Support for singular Z

%function BLOCKS = syncMaster(BLOCKS, flyRecord, options)
%Script/Function for synchronising newtype (2025+) 2p data synchronised with BT/ThorSync

%To do: Check support for follow-on analysis of true bendy rolling data

%           - Integrate andre branch timing mechanics (Use FrameData together with inferTimes?)

%           - Support for blanks

%           - Stimulus intrusion check (Attempt)

%To improve: Absolute best frame specificity for im/stim calcs
    %Specifically:
        %   - Bendy block design has some by-necessity inaccuracies with regards to stimulus aliasing (see stimCollInds/temp4)

    %Also, ability to analyse LED data as if rolling?

%{{
arguments
    %BLOCKS struct
    FLIES struct
    flyRecord table
    options.dataDirectory string = "I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\"
    options.doPlot double = 1
    options.doVid double = 0
    options.rollingAnalysis double = -1
    options.dataSource double = -1 %Whether to use data from function call (-1) or to look external for registered data (1)
    options.allowRandomSequenceEmpty double = 0 %Whether to allow randomSequence to end up empty
    options.disregardRollingDesign double = 0 %Whether to discard bendy rolling design blocks (e.g. When comparing LEDs to bendy block design)
    options.btInterpolationMethod string = 'intelligent' %Whether to use mode, median, last, or an intelligent means to identify stimuli preceding stimulation for block design
        %Note: The 'intelligent' part of said interpolation only fully applies to bendy block design data; Rolling data just uses the stimulus-imaging interpolation method
    options.nBack double = 5 %Functions same as elsewhere; Will override flyRecord nStimuli if it comes to it (e.g. bendy block design)
    options.saveShortcut double = 1 %Whether to save a faster structure than original h5 files to the ThorSync folder
    options.useShortcut double = 1 %Whether to actually use said shortcut
    options.daqFramespikeVoltage double = 1 %Approximate voltage for the FrameSpike voltage in DAQ/TS
    options.blankHandleMode double = 2 %How to deal with blanks when includeBlanks used (1 - Siphon blanks separately, 2 - Treat as normal and siphon at end [Safer])
    options.disregardBattery double = 0 %Whether to discard battery blocks at end, on account of SEs analysis incompatibility
    options.postHocCorrectInferTimes double = 1 %If DAQ data is present, uses DAQ/PTB timings to adjust inferTimes (Since inferTimes is interpolated, not true clock)
    %options.disregardNonBattery double = 0 %Whether to discard NON-battery blocks at end, for simpler battery analysis; Not really functional/useful due to fact that entire analysis has to process before this procs
    options.cleanPhotData double = 1 %Whether to apply esoteric photodiode cleaning operations
    options.doBendyTransientGraph double = 1 %Whether to do an improvised transient (+ phot) graph for bendy block design data
    options.overwriteShortcut double = 0 %Whether to forcibly overwrite shortcut files (Useful after syncMaster changes)
    options.simulationRun double = 0 %Whether to skip the writing of data to FLIES.BLOCKS structure (Prevents bleedover when debugging)
    options.unsiphonedSEs double = 0 %Whether to *not* do any data siphoning/etc for SEs, and just calculate volTimes/etc and leave data trimmed
    options.doBatteryIFICheck double = 0 %Whether to force IFI checks to be done on framespike data for battery (Standard for block design, omitted usually for battery cos freq. condition)
    options.savePhasePlots double = 1 %Whether to force draw/save of phase plots where applicable
end
functionAlity = 1;
%}
%{
if exist('functionAlity')
    ['## Error: Both function and non-function params specified ##']
    crash = yes
end
%BLOCKS = FLIES(fly).BLOCKS;
FLIES = FLIES;
flyRecord = flyRecord;
options.dataDirectory = dataDirectory;
options.doPlot = 1;
options.doVid = 0;
options.rollingAnalysis = -1;
options.dataSource = -1;
options.allowRandomSequenceEmpty = 0;
options.disregardRollingDesign = 0;
options.btInterpolationMethod = 'intelligent';
options.nBack = 5;
options.saveShortcut = 0;
options.useShortcut = 0;
options.daqFramespikeVoltage = 1;
options.blankHandleMode = 2;
options.disregardBattery = 0;
options.postHocCorrectInferTimes = 1;
%options.disregardNonBattery = 0;
options.cleanPhotData = 1;
options.doBendyTransientGraph = 1;
options.overwriteShortcut = 1;
options.simulationRun = 1;
options.unsiphonedSEs = 0;
options.doBatteryIFICheck = 0; 
options.savePhasePlots = 1;
options
%}


%% Prepare

doPlot = options.doPlot;
doVid = options.doVid;
rollingAnalysis = options.rollingAnalysis;
dataSource = options.dataSource;
allowRandomSequenceEmpty = options.allowRandomSequenceEmpty;
disregardRollingDesign = options.disregardRollingDesign;
btInterpolationMethod = options.btInterpolationMethod;
    %Note: Only last is the one to guarantee <block length> number of preceding stimuli will be grabbed
nBack = options.nBack;
saveShortcut = options.saveShortcut;
useShortcut = options.useShortcut;
daqFramespikeVoltage = options.daqFramespikeVoltage; %Called into being even if not applicable
blankHandleMode = options.blankHandleMode;
disregardBattery = options.disregardBattery;
%disregardNonBattery = options.disregardNonBattery;
cleanPhotData = options.cleanPhotData;
doBendyTransientGraph = options.doBendyTransientGraph;
overwriteShortcut = options.overwriteShortcut;
simulationRun = options.simulationRun;
unsiphonedSEs = options.unsiphonedSEs;
doBatteryIFICheck = options.doBatteryIFICheck;
savePhasePlots = options.savePhasePlots;

%Pre-loop preparation
flagParamSaveList = who;
flagParamSaveList = [flagParamSaveList;'flagParamSaveList'; 'thisBlock' ; 'BLOCKS'; 'fly'; 'FLIES']; %Note semicolons, not commas

for fly = 1:length(FLIES) %Need to check this actually does multiple flies

    %Check if actually applicable to synchronise
    if any([FLIES(fly).BLOCKS.isRolling]) || rollingAnalysis == 1%doRolling %May be misusing rollingAnalysis flag
        disp(['(Synchronisation applicable)'])
    else
        disp(['(Synchronisation not applicable)'])
        continue
    end

    BLOCKS = FLIES(fly).BLOCKS;
    disp(['-- Synchronising Fly #',num2str(BLOCKS(1).flyNum),' data --'])
    
    for thisBlock = [BLOCKS]
        clearvars('-except', flagParamSaveList{:}) %Clear variables from previous loop
    
        %thisFlyRowInd = find( flyRecord.Fly == thisBlock.flyNum & flyRecord.Block == thisBlock.blockNum ); %Wrong; Only valid with N=1 fly, due to flyRecord size
        thisFlyRowInd = find( [BLOCKS.flyNum] == thisBlock.flyNum & [BLOCKS.blockNum] == thisBlock.blockNum ); %Fixed; Will also be used later to overwrite data
            %Note: Only relates to position within BLOCKS (i.e. This fly)
        %thisFlyRecord = flyRecord( thisFlyRowInd ,:) %Only correct if running on one fly at a time
        thisFlyRecord = flyRecord( find( flyRecord.Fly == thisBlock.flyNum & flyRecord.Block == thisBlock.blockNum ) ,: )
        %QA
        if isempty(thisFlyRecord) || size(thisFlyRecord,1) > 1
            ['## Alert: Desynchronisation between block and flyRecord information ##']
            crash = yes
        end
        currentDate = char(datetime(thisFlyRecord.Date(1),'Format','dMMMyy'));
        %flyID = ['fly' num2str(thisFlyBlock.FlyOnDay) '_exp' num2str(thisFlyBlock.Block) '_' currentDate];
            
        flyID = thisBlock.flyID;
    
        %Check to see if rolling analysis needs to actually be done
        if isfield(thisBlock,'isRolling') && ~isempty(thisBlock.isRolling) && thisBlock.isRolling == 1
            disp(['-- Proceeding with rolling analysis for ',flyID,' --'])
        else
            disp(['-# Rolling analysis not applicable for ',flyID,' #-'])
            continue
        end

        %Identify ahead of time whether arbitrary phase shift will need to be applied
        if isfield( thisBlock, 'ArbSyncVolShift' ) && ~isempty(thisBlock.ArbSyncVolShift)
            arbPhaseShift = thisBlock.ArbSyncVolShift;
            disp(['Arbitrary phase shift requested (',num2str(arbPhaseShift),' volumes)'])
        end
    
        %Identify folder to look for ThorSync (h5) data in
        %dataFolder = "I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\30Jan25"
        %expName = "fly1_exp2_30Jan25"
        %dataFolder = [strcat("I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\", currentDate)];
        %dataFolder = [strcat(options.dataDirectory, currentDate)];
        dataFolder = char( [strcat(options.dataDirectory, currentDate)] );
        expName = flyID;
    
        %Try find h5 data
        %hFiles = dir( strcat(dataFolder, filesep, 'TS', filesep, expName, filesep, '*.h5') )
        hFiles = dir( strcat(dataFolder, filesep, '**/*.h5') );
        hTarget = [];
        k = 0;
        hFolder = [];
        for i = 1:size(hFiles,1) 
            if contains( hFiles(i).folder, expName ) == 1
                hTarget = strcat( hFiles(i).folder, filesep, hFiles(i).name );
                hFolder = hFiles(i).folder;
                disp(['Appropriate h5 file found in ',char(10),hFiles(i).folder])
                hFolder = hFiles(i).folder; %Necessary for also acquiring settings
                k = k + 1;
            end
        end
        %QA for over/underfind
        if k > 1 || k == 0
            ['## Error: Too many h5 (or none) found ##']
            crash = yes
        end
    
    
        %% Load

        %------------------------------
    
        isShortcutting = 0;
        if useShortcut == 1 && overwriteShortcut ~= 1
            %shortFile = dir( [hFolder,filesep,'h5Shortcut.mat'] );
            shortFileList = dir( strcat(dataFolder,filesep,'SHORT',filesep,'*h5Shortcut.mat') );
            shortFile = [];
            for i = 1:size(shortFileList,1)
                if contains( shortFileList(i).name, expName) && contains( shortFileList(i).name, ['_exp',num2str(thisBlock.blockNum),'_'])
                        %Second portion of boolean required to disambiguate "exp1" vs "exp12" etc
                        %Note: This section hardcodes assumptions about naming structure, obviously
                    shortFile = [shortFile; shortFileList(i).name];
                end
            end
            
            if ~isempty( shortFile )
                %load([shortFile.folder,filesep,shortFile.name])
                load([dataFolder,filesep,'SHORT',filesep,shortFile])
                canShortcut = 1;
                isShortcutting = 1
                disp(['-- Shortcut file successfully loaded --'])
            end

            %Quick QA for whether shortcut file still applicable
            if exist('arbPhaseShift') && ~isfield( shortStruct,'arbPhaseShift' )
                ['-# Alert: Phase shift (now) requested in record, but was not used when shortcut file generated #-']
                crash = yes %Theoretically overkill, but probably good practice
            end
        elseif overwriteShortcut == 1
            saveShortcut = 1;
            disp(['-# Overwriting shortcut files by request #-'])
        end
    
        if useShortcut == 0 || isShortcutting == 0
        %Read ThorSync
        %fileName = "I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\30Jan25\TS\fly1_30Jan25_exp2\Episode001.h5"
        fileName = hTarget;
    
        %h5disp(fileName)
        %heaven
    
        syncStruct = struct;
    
        disp(['-- Loading H5 data --'])
        tic
    
        %syncStruct.AI.piezoData = h5read(fileName, '/AI/PiezoMonitor'); %Use in future?
        syncStruct.CI.frameData = h5read(fileName, '/CI/FrameCounter');
        %syncStruct.DI.bleachComData = h5read(fileName, '/DI/BleachComplete');
        syncStruct.DI.bleachOutData = h5read(fileName, '/DI/BleachOut');
        %syncStruct.DI.frameInData = h5read(fileName, '/DI/FrameIn');
        syncStruct.DI.frameOutData = h5read(fileName, '/DI/FrameOut');
        %syncStruct.DI.pmtData = h5read(fileName, '/DI/PMTShutter');
        %syncStruct.Freq.freqFitData = h5read(fileName, '/Freq/FitHz');
        %syncStruct.Freq.freqHzData = h5read(fileName, '/Freq/Hz');
        %syncStruct.Global.globData = h5read(fileName, '/Global/GCtr');
        try
            syncStruct.AI.FrameSpike = h5read(fileName, '/AI/FrameSpike');
            syncStruct.AI.Iterator = h5read(fileName, '/AI/Iterator');
            hasDaqData = 1;
            disp(['-- TS DAQ data retrieved --'])
        catch
            ['-# Could not retrieve DAQ-related TS data #-']
            hasDaqData = 0;
        end
        try
            syncStruct.AI.Photodiode = h5read(fileName, '/AI/Photodiode');
            hasPhotData = 1;
            disp(['-- TS Photodiode data retrieved --'])
        catch
            ['-# Could not retrieve photodiode data #-']
            hasPhotData = 0;
        end
    
        toc
        %plastic
    
        %Also attempt to acquire TS settings
        %h5SetFile = dir([ hFiles(i).folder, filesep, 'ThorRealTimeDataSettings.xml' ]) %Relies on i not being used between this and H5 finding
        h5SetFile = dir([ hFolder, filesep, 'ThorRealTimeDataSettings.xml' ]); %Relies on i not being used between this and H5 finding
        %QA
        if isempty(h5SetFile)
            ['## Error in finding H5 settings file ##']
            crash = yes
        end
        h5Sets = readstruct( [h5SetFile.folder,filesep,h5SetFile.name] );
        sampRate = [];
        for daqInd = 1:size( h5Sets.DaqDevices.AcquireBoard,2 )
            if h5Sets.DaqDevices.AcquireBoard(daqInd).activeAttribute == 1
                disp(['Active ThorSync DAQ device: ',h5Sets.DaqDevices.AcquireBoard(daqInd).typeAttribute{1} ])
                for sampInd = 1:size( h5Sets.DaqDevices.AcquireBoard(daqInd).SampleRate,2 )
                    if h5Sets.DaqDevices.AcquireBoard(daqInd).SampleRate(sampInd).enableAttribute == 1
                        sampRate = h5Sets.DaqDevices.AcquireBoard(daqInd).SampleRate(sampInd).rateAttribute;
                        disp([ 'Sampling rate: ', num2str( sampRate ),...
                            ' ("',h5Sets.DaqDevices.AcquireBoard(daqInd).SampleRate(sampInd).nameAttribute{1},'")' ])
                    end
                end
            end
        end
        %QA
        if isempty(sampRate) 
            ['-# Sampling rate could not be derived; Using default 30kHz #-']
            sampRate = 30000;
        end
        %Infer TS timescale
        inferTimes = linspace( 0, size( syncStruct.CI.frameData , 2 ) / sampRate, size( syncStruct.CI.frameData , 2 ) );
            %Linearly space time from 1st to last element of frameData
            %The division of total frame count by sampling rate is weak to (TS) framedrops/etc, but no better option present
        end
    
        %Read experiment information file from ThorImage
        %Note: Requires >=MATLAB 2020b
        eFile = dir( strcat(dataFolder, filesep, expName,filesep,'Experiment.xml') );
        %QA
        if isempty( eFile )
            ['## Experiment.xml file not found! ##']
            crash = yes
        end
        expStruct = readstruct( [eFile.folder,filesep,eFile.name] ); %Structure contains all information ThorImage knows about the recording
    
        imStartTime = expStruct.Date.uTimeAttribute;
            %Note: This is confirmed to be when Start button was pressed, not when imaging commenced
                %i.e. Disparity of 1 - 30s
    
        if useShortcut == 0 || isShortcutting == 0        
        %ThorImage (apparently)
        frameStates = unique( syncStruct.DI.frameOutData );
        %QA
        if size(frameStates,2) > 3
            ['-# Caution: Potentially aberrant number of unique frameOut states found #-']
        end
        frameEndInd = find( syncStruct.DI.frameOutData == nanmax(frameStates), 1, 'last' ); %Find last instance of high state
    
        %Arduino
        bleachStates = unique( syncStruct.DI.bleachOutData );
        %QA
        if size(bleachStates,2) > 3
            ['-# Caution: Potentially aberrant number of unique bleachOut states found #-']
        end
        bleachEndInd = find( syncStruct.DI.bleachOutData == nanmax(bleachStates), 1, 'last' );
        end

        %------------------------------
    
        %Read MATLAB params and BT data
        %Find MATLAB params
        %The following block relies heavily on two assumptions:
        %   Image folder naming architecture is " flyX_expY_<Date><Month><Year> " e.g. fly1_exp2_30Jan25
        %   MAT param naming system is " MAT_<datecode>_fX_bY " e.g. MAT_300125_f1_b2
        expNameMAT = strrep( expName, 'fly', 'f' );
        expNameMAT = strrep( expNameMAT, 'exp', 'b' );
        expNameMAT = strsplit(expNameMAT,'_');
        %MATDate = datestr( datetime("30Jan25",'InputFormat','ddMMMyy'), 'ddmmyy' );
        MATDate = datestr( datetime(expNameMAT{3},'InputFormat','ddMMMyy'), 'ddmmyy' );
        %Read params
        targetMATName = strcat( 'MAT_', MATDate, '_', expNameMAT{1}, '_', expNameMAT{2} );
        matParamFile =dir( strcat(dataFolder, filesep, '**/', targetMATName, '.mat') );
        altTargetMATName = strcat( MATDate, '_', expNameMAT{1}, '_', expNameMAT{2},'_MAT' );
        altMatParamFile =dir( strcat(dataFolder, filesep, '**/', altTargetMATName, '.mat') );
        ancMatName = strcat( MATDate, '_', expNameMAT{1}, '_', expNameMAT{2},'_btAncillary' );
        ancMatFile = dir( strcat(dataFolder, filesep, '**/', ancMatName, '.mat') );
        %QA
        if isempty( matParamFile ) && isempty(altMatParamFile)
            ['## PTB MATLAB parameters not found! ##']
            %crash = yes
            hasPTB = 0;
        else
            hasPTB = 1;
            disp(['PTB MATLAB parameters found'])
        end
        %QA for ancillary/load
        if ~isempty( ancMatFile )
            btAncillary = load( [ancMatFile.folder,filesep,ancMatFile.name] );
            btAncillary = btAncillary.btAncillary;
            disp(['btAncillary loaded'])
        end
        progVerNum = NaN; %Default, will be calculated later if possible
        if hasPTB
            %matParamStruct = load( [matParamFile.folder,filesep,matParamFile.name] );
            if ~isempty( matParamFile )
                disp(['(Legacy name)'])
                matParamStruct = load( [matParamFile.folder,filesep,matParamFile.name] );
            else
                disp(['(Newtype name)'])
                matParamStruct = load( [altMatParamFile.folder,filesep,altMatParamFile.name] );
            end
    
            %Quick check to guard against accidental unitary analysis
            if isequal( matParamStruct.matSave.stimuli, 'unitary' )
                ['-# Caution: Unitary stimulus identity present #-']
                crash = yes %Technically forgivable, but probably want to not analyse for meaningless SEs
            end
    
            %Read BT
            targetBTName = strcat( MATDate, '_', expNameMAT{1}, '_', expNameMAT{2},'_btData' );
            btFile =dir( strcat(dataFolder, filesep, '**/', targetBTName, '.csv') );
            %QA
            if isempty( btFile )
                ['## BT data not found! ##']
                crash = yes
            end
            btData = csvread( [btFile.folder,filesep,btFile.name] );
            %Values correspond to (As of v8.1):
            % 1 - BTtime, 2 - pMarkTime, 3 - dWLE, 4 - btHash, 5 - i, 6 - currentTime, 7 - repCount, 8 - randomSequence(i), 9 - onOff
            %  BTtime (BT time [Since start?]), pMarkTIme (Time since script start?), dWLE (Movement?), btHash (Hash), i (Sequence element #), currentTime, repCount, actual sequence element, flip status (where applicable)
    
            %Clean up
            %if nansum( isnan( btData(end,:) ) ) == numel( btData(end,:) )
            while nansum( isnan( btData(end,:) ) ) == numel( btData(end,:) )
                btData(end,:) = [];
                disp(['Terminal btData element removed due to all-NaN nature'])
            end
            %QA
            if nansum( isnan(btData), 'all' ) > 0
                disp(['-# Caution: NaN elements found in btData #-'])
                %crash = yes %Again, not critical but problematic if not at end
                    %Nowadays likely if not battery etc
            end
            randomSequence = matParamStruct.matSave.randomSequence( unique(btData(:,5)) );
                %Note: Since posthoc derived, may deviate from theoretical
    
            %Acquire information with regards to potential bendy panel block design
            bendyBlockDesign = 0; %Default no (Rolling)
            if isfield( matParamStruct.matSave, 'blockDesign' ) && matParamStruct.matSave.blockDesign == 1
                bendyBlockDesign = 1; %Block
                disp(['Bendy panel block design detected'])
            else
                bendyBlockDesign = 0; %Rolling
                disp(['Bendy panel rolling design (presumably) detected'])
            end

            %Check if battery
            batteryDesign = 0; %Default no
            if isequal( matParamStruct.matSave.stimuli , "battery")
                batteryDesign = 1;
                disp(['~ Battery design detected ~'])
            end

            %Get script version
            if isfield( matParamStruct.matSave, 'progIdent' )
                progIdent = matParamStruct.matSave.progIdent;
            else
                ['-# Could not detect script version #-']
                progIdent = [];
            end

            %Check if blanks
            if isfield( matParamStruct.matSave, 'includeBlanks' )
                includesBlanks = matParamStruct.matSave.includeBlanks;
            else
                includesBlanks = 0; %Default no
            end
    
            ptbStartTime = btData(1,2)/1000; %Divide by 1000 because PTB stores down to millisecond unix time
            %ptbEndTime = btData( find( ~isnan( btData(:,2) ) , 1, 'last'  ) ,2)/1000; %Find necessary because NaN last element
            ptbEndTime = btData( [end-find( diff( flip( btData(:,4) )) ~= 0, 1, 'first')] ,2)/1000; %More accurate, use last non-repetitious time
                %Estimated accuracy: Effectively 100% (~millisecond), no drift
                    %Note: Datetime/etc does not show milliseconds, but this value is (seemingly) still stored
                        %i.e. There is a ~0.4s difference between the last non-repetitious and the 'end' of the data
                    %Secondary note: Relies heavily on the assumption that arduino low -> bt repetitiveness instantly
    
             if ~isShortcutting        
                 if ~batteryDesign
                     estimatedPTBEndTime = imStartTime + inferTimes( bleachEndInd ); %Use arduino low to estimate PTB end time in TS reference
                 else
                     estimatedPTBEndTime = ptbEndTime; %Use PTB, because arduino likely not applicable
                 end
             else
                estimatedPTBEndTime = shortStruct.estimatedPTBEndTime;
             end
                %This is a time guessed from the self-reported imaging start time
                %It is intended to be used to QA ptbEndTime in case of not ending because repetitive, not a true comparison
              %QA
              if abs( ptbEndTime - estimatedPTBEndTime ) > 0.05 * ( estimatedPTBEndTime - imStartTime )  %"Did PTB end more than 5% +- Arduino time?"
                  ['## Alert: PTB self-reported end time differs significantly from Arduino estimated end time ##']
                  crash = yes %Not technically critical, but probably worrying
              end
    
              %Find out if rolling experiment (unless battery)
              if rollingAnalysis == -1 && isfield( matParamStruct.matSave, 'blockDesign' ) && ~batteryDesign
                  if matParamStruct.matSave.blockDesign == 1
                      rollingAnalysis = 0; %Override
                      disp(['Block design detected in PTB params; Using block analysis'])
                  else
                      rollingAnalysis = 1;    
                      disp(['Rolling design detected in PTB params; Using rolling analysis'])
                  end
              end

              %Derive script version
              progIdent = matParamStruct.matSave.progIdent;
              progSplit = split( progIdent, '_');
              progVerNum = NaN;
              for i = 1:size( progSplit,1 )
                  if isequal( progSplit{i}(1) , 'v' )
                      %citadel
                      temp =  progSplit{i};
                      temp = strrep( temp, 'v', '' );
                      temp = strrep( temp, 'dot', '.' );
                      progVerNum = str2num( temp );
                  end
              end
              if ~isempty(progVerNum) && ~isnan(progVerNum)
                  disp(['Derived stimulus script version: ',num2str(progVerNum)])
              else
                  disp(['Failure to derive script version #'])
              end
    
        end

        %------------------------------

        %Little follow-on
        if isShortcutting
            hasDaqData = shortStruct.hasDaqData;
            hasPhotData = shortStruct.hasPhotData; %If crashes, delete old shortcute files and regen
        end
        %DAQ
        if hasDaqData && ~isShortcutting
            %freak

            %---------------------------------
            %Notes on choreography:
            %   - Iterator/FrameSpike: The framespike should always occur at the *end* of an iterator 'block' (i.e. Iterator sits at 135 [Most recent framespike was 134], then framespike #135 occurs and iterator changes to 136, and so on)
            %   - Iterator terminal values: Often the last iterator change will be detected as an iteration, and this is mostly excluded
            %   - Proximal framespikes: With 8.45, a framespike is sent at the start of the script, but since this causes confusion, it is typically removed. For battery, even more are removed (relating to opto)
            %---------------------------------

            %Acquire post-hoc counts of framespike and iterator (Script >v8.55)
            daqFrameSpikeCount = NaN;
            daqIteratorCount = NaN;
            if exist('btAncillary') && isfield( btAncillary, 'daqFrameSpikeCount' )
                daqFrameSpikeCount = btAncillary.daqFrameSpikeCount;
                daqIteratorCount = btAncillary.daqIteratorCount;
            end

            %Precalculate useful number
            if progVerNum < 8.45 %Legacy script, no framespike/iterator initialisation at while ~quit_now commencement
                targetINum = btData(end,5)-1;
            else %Newtype script, with initialisation framespike/etc
                targetINum = btData(end,5); 
            end
            if batteryDesign
                targetINum = targetINum - matParamStruct.matSave.optoElements; %Because opto no framespikes
            end

            nonLinearityPresent = 0;

            %Find framespike peaks
                %NOTE: Current system optimised towards local peak of framespike, whereas 'true' timing of stimulus onset probably lies at onset of voltage increase
            if hasPTB && ~isequal( matParamStruct.matSave.stimuli, 'battery' )
                [framePKS, frameLOCS] = findpeaks( syncStruct.AI.FrameSpike, 'MinPeakHeight', 0.5*daqFramespikeVoltage,...
                   'MinPeakDistance', ( sampRate / matParamStruct.matSave.frequency )*0.45 ); %Look for peaks separated by at 45% of a flip; Adjusted down to account for some rare instances where framespike shifted forwards in time
                   %'MinPeakDistance', ( sampRate / matParamStruct.matSave.frequency )*0.66 ); %Look for peaks separated by at 66% of a flip
            else
                [framePKS, frameLOCS] = findpeaks( syncStruct.AI.FrameSpike, 'MinPeakHeight', 0.5*daqFramespikeVoltage, ...
                    'MinPeakWidth', 200); %Use an empirical spike width based on pilot data (Probably good up to 60Hz framespikes)
                        %Note: Theoretically could make this dynamic for battery by multiplying fastest freq element by sampRate, but cbf
            end

            %Omit first framespike if Mk 8.45+ (Because 8.45 added a preliminary framespike immediately prior to while loop commencement)
            omittedFrameLOCPK = [];
            if progVerNum >= 8.45 %&& ~batteryDesign %Adjusted to now apply to battery
                omittedFrameLOCPK = [frameLOCS(1),framePKS(1)];
                frameLOCS(1) = [];
                framePKS(1) = [];
                disp(['-- First framespike omitted wrt v8.45+ --'])
            end
            %QA
            if ~isnan( daqFrameSpikeCount ) && daqFrameSpikeCount ~= size(frameLOCS,2)+1
                    %Only 80% confident the initialisation framespike is counted as one...
                ['## Alert: Disparity between post-hoc reported framespike count and detected framespike count ##']
                %if ~batteryDesign
                %crash = yes %Removed this so that later iterator QA can back up whether to crash
                %end
            end

            %Plot
            if doPlot && isfield( syncStruct.AI, 'FrameSpike' )
                figure
                subplot(2,1,1)
                theseInds = [1:frameLOCS(10)];
                plot( inferTimes(theseInds) , syncStruct.AI.FrameSpike(theseInds) )
                hold on
                plot( inferTimes(theseInds) , syncStruct.AI.Iterator(theseInds), 'r' )
                scatter( inferTimes(frameLOCS(1:10)), framePKS(1:10) )
                for i = 1:10
                    text( inferTimes(frameLOCS(i)), framePKS(i)*1.1, [num2str(i)], 'Color', 'r' )
                end
                %{
                if progVerNum >= 8.45
                    scatter( inferTimes(frameLOCS(1)), framePKS(1)*1.05, 'Color', [1,0,0])
                    textVals = [0:9];
                else
                    textVals = [1:10];
                end
                for i = 1:size( textVals,2 )
                    text( inferTimes(frameLOCS(i)), framePKS(i)*1.1, [num2str(textVals(i))], 'Color', 'r' )
                end
                %}
                title([strrep(expName,'_',' '),' - First 10 framespike elements + Iterator'])
                xlabel('Times (s)')
                ylabel('Voltage (V)')
                subplot(2,1,2)
                theseInds = [frameLOCS(end-10):frameLOCS(end)+sampRate];
                plot( inferTimes(theseInds) , syncStruct.AI.FrameSpike(theseInds) )
                hold on
                plot( inferTimes(theseInds) , syncStruct.AI.Iterator(theseInds), 'r' )
                scatter( inferTimes(frameLOCS(end-10:end)), framePKS(end-10:end) )
                for i = length(frameLOCS)-10:length(frameLOCS)
                    text( inferTimes(frameLOCS(i)), framePKS(i)*1.1, [num2str(i)], 'Color', 'r' )
                end
                title([strrep(expName,'_',' '),' - Last 10 framespike elements + Iterator'])
                xlabel('Times (s)')
                ylabel('Voltage (V)')
            end

            %Find iterator values
                %Note: As currently written (Mk 8.55), the iterator value will increase 'aberrantly' as a last act before the while loop quits
            %taste
            daqRange = matParamStruct.matSave.daqRange;
            daqItDiff = nanmedian( diff( daqRange ) );
            if round( daqItDiff, 1) == 0.1
                %This section works by effectively smoothing/cleaning the iterator values then attempting to match them to daqRange and then bootleg unwrapping that
                disp(['DAQ iterator values distributed by 0.1'])
                daqItDiff = 0.1;
                %itVals = round( syncStruct.AI.Iterator, 1 );
                itVals = smooth( syncStruct.AI.Iterator, 100 ); %Value of 100 empirically determined based on apparent DAQ voltage kinematics (~500 sample rise time for standard step)
                itVals( itVals <= (daqRange(1) - 0.5*daqItDiff)  ) = NaN;
                %itVals( itVals <= (daqRange(1) )  ) = NaN; %Aggressive cleaning
                %itVals = floor( syncStruct.AI.Iterator*10 )/10;
                %{
                itVals = round( itVals, 3, 'significant' );
                itVals = smooth( itVals, 100 );
                itVals = round( itVals, 3, 'significant' );
                itVals( itVals <= (matParamStruct.matSave.daqRange(1) - daqItDiff)  ) = NaN;
                %}
                itValsDiff = abs( diff(itVals) );
                if hasPTB && ~isequal( matParamStruct.matSave.stimuli, 'battery' )
                    [itPKS, itLOCS] = findpeaks( itValsDiff, 'MinPeakHeight', 0.25e-3,... %Empirical
                       'MinPeakDistance', ( sampRate / matParamStruct.matSave.frequency )*0.66 ); %Look for peaks separated by at 66% of a flip
                else
                   % [itPKS, itLOCS] = findpeaks( itValsDiff, 'MinPeakHeight', 0.25e-3, ...
                   %     'MinPeakWidth', 200); %Use an empirical spike width based on pilot data
                   [itPKS, itLOCS] = findpeaks( itValsDiff, 'MinPeakHeight', 0.25e-3, ...
                        'MinPeakWidth', 25, 'MinPeakDistance',  (1/nanmax( matParamStruct.matSave.freqSequence ))*sampRate*0.5); %Use an empirical spike width based on pilot data together with a derived inter-interval change distance
                end
                    %Note: Will be 1 element longer than true iteration value (See above)
                %QA for emptiness
                if isempty( itLOCS )
                    ['## Alert: No discernible values found in itVals differences (-> itLOCS) ##']
                    crash = yes
                end

                 %Remove last itLOCS value, since extremely high likelihood script ended in this period
                 terminalItDropped = 0; %Will store whether last iterator value dropped (And thus, flipOnsetIndices will likely need trimming too)
                 if frameLOCS(end) < itLOCS(end) %But only do it (more or less safely) if last index less than last index of framespikes
                     itLOCS(end) = []; %Note: If PTB ended mid iterator block, it is possible that the last framespike could precede *2* (apparent) changes in iterator (1 true and 1 EOF)
                     disp(['Terminal itLOCS value dropped'])
                     terminalItDropped = 1; %Note: Must represent how many elements removed so that QA below works
                 end

                 %battery delete first 2(?) elements of itLOCS cos opto
                    %Based on empirical of May 29th starting sequence
                 if batteryDesign && progVerNum >= 8.45
                    itLOCS(1:2) = []; %Empirical, under observation that Left-Right begins with a rapid 2x framespike/iteration combo
                                                    %Iterator stable -> Iterator low -> Iterator 'real' -> Next iterator, in quick succession
                    frameLOCS(1) = []; 
                    disp(['Battery: Proximal 2 (for now) elements of itLOCS and 1 elements of frameLOCS deleted'])
                 elseif batteryDesign && progVerNum < 8.45
                     ['manual case curation needed wrt no initialisation framespike']
                     todo = yes
                 end

                %QA for sanity
                if abs( targetINum - length(itLOCS) ) > 0.05*targetINum
                    ['-# Alert: >5% difference in detected iterator changes and target number of i values #-']
                    crash = yes
                end
                itProc = []; %Will be 1 shorter than itLOCS eventually, but hopefully the same number of elements as max of i
                itProcIdealised = [];
                itProcRangeInd = [];
                for i = 1:length(itLOCS)%-1
                    %itProc(i) = mode( round(itVals( itLOCS(i):itLOCS(i+1) ),3 ,'significant') ); %Not valid
                    if i < length(itLOCS) %Normal
                        itProc(i) = mode( floor(itVals( itLOCS(i):itLOCS(i+1) )*10)/10 ); %Note: Use of floor only valid where 2nd decimal place meaningless and does not exceed 0.1 (Empirically true within DAQ range of 0.5 to 4.5)
                    else %Last element
                        itProc(i) = mode( floor(itVals( itLOCS(i):itLOCS(i)+(itLOCS(i)-itLOCS(i-1)) )*10)/10 ); %Use assumption of normal size and forward-extrapolate (Note: Relies on NaNing of last elements to keep this value sane)
                            %Note: With current itLOCS last element dropping there might be a better way to do this, but cbf currently
                    end
                    [~,itProcRangeInd(i)] = min( abs( repmat( itProc(i), 1, length(daqRange) ) - daqRange ) );
                    itProcIdealised(i) = daqRange( itProcRangeInd(i) );
                end
                %QA to see if discontinuity
                if any( diff( itProcIdealised ) > 0.11 )
                    ['-# Alert: Apparent discontinuity (Jump >0.11) in idealised itProc data #-']
                    crash = yes
                    %This is most likely to be because itLOCS failed to detect a change in itVals, or detected too many etc
                end                
                %itProc(i+1) = mode( round(itVals( itLOCS(end):itLOCS(end)+nanmedian( diff(itLOCS) ) ),3 ,'significant') ); %Note: V. high likelihood of this number being aberrant if itVals go low state immediate after; Round not good
                %itProc(i+1) = mode( floor(itVals( itLOCS(end):itLOCS(end)+nanmedian( diff(itLOCS) ) )*10)/10 ); %Note: V. high likelihood of this number being aberrant if itVals go low state immediate after
                %[~,itProcRangeInd(i+1)] = min( abs( repmat( itProc(i+1), 1, length(daqRange) ) - daqRange ) );
                %itProcIdealised(i+1) = daqRange( itProcRangeInd(i+1) );
                    %By not doing these ~3 lines the last, 'aberrant' element is ommitted
                itProcRangeIndUnwrapped = [];
                itProcRangeIndUnwrapped(1) = itProcRangeInd(1);
                roller = 0;
                rollCount = 0;
                for i = 2:length(itProcRangeInd)
                    if itProcRangeInd(i) < itProcRangeInd(i-1)
                        roller = roller + length(daqRange);         
                        rollCount = rollCount + 1;
                    end
                    itProcRangeIndUnwrapped(i) = itProcRangeInd(i) + roller;
                end

                %more battery stuff: subtract first element (then +1) due to non-zero ("0.5") starting point
                if batteryDesign
                    itProcRangeIndUnwrapped = itProcRangeIndUnwrapped - itProcRangeIndUnwrapped(1) + 1;
                    %Note: May behave unpredictably if first iterator is a loop (4.5 -> 0.5)
                    disp(['Battery: Iterator values reset to 1 at first non-opto position'])
                end

                %QA
                if numel(unique(itProcRangeIndUnwrapped)) ~= numel(itProcRangeIndUnwrapped)
                    ['-# Alert: Discontinuity detected in unwrapped iterator values #-'] %Ideally shouldn't happen, but may
                        %Likely reasons are lag, etc
                    crash = yes
                end
                if ~batteryDesign && exist('btAncillary') && isfield(btAncillary,'daqLoopCount') && rollCount ~= btAncillary.daqLoopCount %May need to allow disparity; Also, not applicable to battery (Due to opto silent rolling?)
                    ['-# Asynchrony between reported number of iterator loops and detected number of loops #-']
                    [num2str(btAncillary.daqLoopCount),' vs ', num2str(rollCount)]
                end
                if ~isnan( daqIteratorCount ) && ( daqIteratorCount ~= itProcRangeIndUnwrapped(end)+terminalItDropped && daqIteratorCount ~= itProcRangeIndUnwrapped(end) ) %Allow for disparity of one (See below terminal cases)
                        %Only 80% confident the initialisation framespike is counted as one...
                    ['## Alert: Disparity between post-hoc reported iterator count and detected iterator count ##']
                    crash = yes
                end

                %Plot
                if length(itLOCS) == length(itProcRangeIndUnwrapped)
                    plotInds = [1:10; itProcRangeIndUnwrapped(end)-10+1:itProcRangeIndUnwrapped(end) ];
                    figure
                    for i = 1:size(plotInds,1)
                        subplot( size(plotInds,1), 1, i )
                        theseCoords = plotInds(i,:);
                        firstLast = [ itLOCS( theseCoords(1) ),...
                            itLOCS( theseCoords(end) )];
                        %firstLastCoords = [ find( inferTimes >= inferTimes( firstLast(1) )-1, 1 , 'first' ),...
                        %    find( inferTimes >= inferTimes( firstLast(2) )+1, 1 , 'first' ) ];
                        %plot( [firstLastCoords(1):firstLastCoords(2)] , itVals( firstLastCoords(1):firstLastCoords(2) ) )
                        plot( [ inferTimes( firstLast(1):firstLast(2) ) ],...
                            [ itVals( firstLast(1):firstLast(2) ) ] )
                        hold on
                        for subI = 1:length(theseCoords)
                            thisSub = theseCoords(subI);
                            text( [inferTimes(itLOCS(thisSub))], [ itVals(itLOCS(thisSub))+1*daqItDiff ], [ num2str(thisSub) ], 'Color', 'r' )
                        end
                        xlabel(['Time (s)'])
                        ylabel(['Voltage (V)'])
                    end
                else
                    ['-# Cannot do iterator plot with dissimilar length itLOCS/unwrapped values #-']
                end

                %QA to compare iterator against framespike
                if abs( frameLOCS(1) - itLOCS(2) ) > nanmedian( diff(frameLOCS) ) * 0.75
                    ['-# Alert: Apparent >n+1 gap between first framespike and supposedly following iterator change #-']
                    crash = yes %Highly like if 'first' framespike detected is actually a little way into experiment or similar 
                    %Goes hand in hand with framespike detection numbers being wrong
                end

                %Make simpler
                lastIt = itProcRangeIndUnwrapped(end);

                %Report
                disp(['Final calculated iterator value: ',num2str(lastIt)])
                
            else
                ['-# Cannot detect DAQ iterator step value! #-']
                crash = yes %Maybe overkill
            end

            %shadow

            %OMEGA PLOT
                %DANGER - 12:30 (As in, it takes 12.5m [Not really] to generate this plot or change the limits)
                %{
            figure
            plot( inferTimes, syncStruct.AI.FrameSpike )
            hold on
            plot( inferTimes, syncStruct.AI.Iterator )
            for i = 1:length(frameLOCS)
                text( inferTimes(frameLOCS(i)), framePKS(i)*1.1, [num2str(i)], 'Color', 'r' )
            end
            if length(itLOCS) == length(itProcRangeIndUnwrapped)
                for i = 1:length(itLOCS)
                    text( inferTimes(itLOCS(i)), itVals(itLOCS(i))+daqItDiff, [num2str(i)], 'Color', 'g' )
                end
            else
                ['-# Cannot (safely) report iterator values with dissimiliarity #-']
            end
                %}

            %--------------------
            %Check for IFI loss and correct first
            temp2 = diff( frameLOCS ); %Inter-framespike index difference; Again, should be (relatively) stable
            ifiLoss = 0;
            if nanmax( abs(temp2) ) > nanmean(temp2) + 4*nanstd(temp2) %Note: Calcs may not cover symmetrical conditions of frameSpike shift
                %Note: No explicit inter-framespike interval rectification currently implemented, only phase loss
                %errorSpikes = find( abs(temp2) > nanmean(temp2) + 4*nanstd(temp2) ); %SD-based system
                errorSpikes = find( abs(temp2 - nanmean(temp2))  >  0.1*nanmean(temp2)  );
                errorSpikesMag = abs(temp2(errorSpikes)) - nanmean(temp2);
                if exist('daqFrameSpikeCount')
                    ['-# Alert: ',num2str(length(errorSpikes)),' significant deviation/s in inter-framespike interval detected (of ~',num2str(daqFrameSpikeCount),' total) #-']
                else
                    ['-# Alert: ',num2str(length(errorSpikes)),' significant deviation/s in inter-framespike interval detected #-']
                end
                if length(errorSpikes) < 10
                    disp([num2str(errorSpikes)])
                else
                    disp([num2str(errorSpikes(1:5)),'...',num2str(errorSpikes(end-5:end))])
                end
                if ~batteryDesign || ( batteryDesign && doBatteryIFICheck )
                    disp([ '(Error spike/s deviated by ',num2str(round(errorSpikesMag)),' TS frames from +-10% mean IFI threshold [',num2str(round(0.1*nanmean(temp2))),'])' ])
                end
                ifiLoss = 1;
            end
            if exist('errorSpikes') && size( errorSpikes, 2 ) > 0 && ( ~batteryDesign || ( batteryDesign && doBatteryIFICheck ) )
                disp(['IFI loss in framespikes detected'])

                %Plot
                figure
                if size( errorSpikes, 2 ) > 10
                    eroInds = errorSpikes( floor([1:length(errorSpikes)/10:length(errorSpikes)]) );
                    ['-# Excessive number of error spikes; Using 10 across range #-']
                    disp(num2str(eroInds))
                else
                    eroInds = 1:size( errorSpikes, 2 );
                end
                a = 1;
                for eros = eroInds%1:size( errorSpikes, 2 ) 
                    %subplot( 1, size( errorSpikes, 2 ) , eros  )
                    subplot( 1, size( eroInds, 2 ) , a  )
                    thisErrorSpike = errorSpikes(eros);
                    spikeRange = [ thisErrorSpike-2, thisErrorSpike+2 ];
                    %spikeRange( spikeRange < 0 ) = [];
                    spikeRange( spikeRange < 0 ) = 1; %Use framespike 1, rather than null
                    spikeRange( spikeRange > size(frameLOCS,2) ) = [];
                    plot( [inferTimes( frameLOCS(spikeRange(1)):frameLOCS(spikeRange(2)) )],...
                        [syncStruct.AI.FrameSpike( frameLOCS(spikeRange(1)) : frameLOCS(spikeRange(2)) )] )
                    hold on
                    scatter( [inferTimes( frameLOCS( spikeRange(1):spikeRange(2) ) )],...
                        [syncStruct.AI.FrameSpike( frameLOCS(spikeRange(1):spikeRange(2))) ] )
                    scatter( [inferTimes( frameLOCS( thisErrorSpike ) )], [syncStruct.AI.FrameSpike( frameLOCS( thisErrorSpike ) )]+0.1, 'Color', [1,0,0] )
                    title(['IFI loss frameSpike (#',num2str(errorSpikes(eros)),') location'])
                    xlabel(['Time (s)'])
                    a = a + 1;
                end

                %Identify
                errorLOCS = []; %Will store deleted LOCS for posterity
                rectifiedErrorSpikes = []; %Will store errorSpikes that are no longer valid
                    %Note: Currently not done with a while loop/etc, so there is no guarantee that framespikes following rectification will actually be valid
                for eros = 1:size( errorSpikes,2 )
                    thisErrorSpike = errorSpikes(eros);
                    disp(['Processing error spike ', num2str(eros),' (#',num2str(thisErrorSpike),')'])
                    %Pre-QA
                    if thisErrorSpike == 1
                        ['## Alert: Unsafe to perform IFI rectification on first framespike ##']
                        %I mean, theoretically it's no different, but in practice one would want to do it with care
                        crash = yes
                    end
                    if size( errorSpikes,2 ) >= 2 && eros <= size( errorSpikes,2 )-1 && thisErrorSpike+1 ==  errorSpikes(eros+1) %"More than 1 error spike AND not last AND this error spike is directly followed (in frameLOCS) by another error spike"
                        %This is designed for the case where the framespike writing lagged out and took ~100+ms to write and unwrite, thus leading to an aberrantly broad framespike
                            %Note: This is not designed to deal with 3+ error spikes in quick succession
                            errorLOCS = thisErrorSpike;
                            frameLOCS( thisErrorSpike ) = []; %"Delete the first of the aberrantly double-detected framespikes"
                            framePKS( thisErrorSpike ) = [];
                            rectifiedErrorSpikes = [rectifiedErrorSpikes,errorSpikes(eros+1)];
                            disp(['-# Framespike #', num2str(thisErrorSpike),' deleted to rectify phase #-'])
                                %Note: Currently as written the actual framespike deleted here doesn't matter **theoretically**
                                    %This is because aside from Framespike 1 and last, the physical position of individual framespikes is not used for anything
                                        %Secondary note: This is not true wrt phase, and so it is possible that picking one of an aberrant double could be better/worse for preserving phase
                    elseif ismember( thisErrorSpike, rectifiedErrorSpikes )
                        disp(['-# Error spike already rectified #-'])
                        continue
                    elseif abs( errorSpikesMag(eros) ) > 1.25*nanmean(temp2)
                        ['## likely missed framespike; case not written yet ##']
                        crash = yes
                    elseif abs( (frameLOCS( thisErrorSpike+1 ) - frameLOCS( thisErrorSpike )) - (itLOCS( thisErrorSpike+2 ) - itLOCS( thisErrorSpike+1 )) ) < 0.01*nanmean(temp2)
                        %Note: This case may have to be broadened for 'normal' phase loss (i.e. use itLocs +1 and 0 respectively, or vice versa)
                        disp(['-# IFI stutter and iterator stutter appear identical; Not rectifying #-'])
                            %Phase loss checks will be important to ensure this does not cause issues
                    else %CHECK FOR PHASE?
                        ['## case not written yet ##']
                        %This will be something like a framespike being just really delayed or ahead in time; TBD what to do
                        crash = yes
                    end
                end

            end

            %--------------------

            %Check for phase loss/other irregularities in data
            temp = []; %Index difference between framespikes and iterator changes; Should be stable
            for i = 1:nanmin( [size(frameLOCS,2), size( itLOCS,1)] ) %Note: itLOCS not 100% guaranteed 1:1 with itProcRangeIndUnwrapped
                temp(i) = frameLOCS(i) - itLOCS(i);
            end
            meanPhaseDiff = nanmean(temp);
            stdPhaseDiff = nanstd(temp);
            %temp = diff(temp); %Diff, better for mean/SDing

            phaseLoss = 0;
            if any(temp < 0)%nanmax( abs(temp) ) > nanmean(temp) + 4*nanstd(temp)
                ['-# Alert: Phase loss detected between framespike and iterator #-']
                phaseLoss = 1;
            end

            %Rectify
            if phaseLoss %&& ifiLoss %FrameSpike positions both undertake/overtake iterator AND unequally spaced
                not modified yet wrt ifi prerectification
                %errorSpikes = find( abs(temp2) > nanmean(temp2) + 4*nanstd(temp2) );
                phaseLossSpike = errorSpikes( find( frameLOCS(errorSpikes) - itLOCS(errorSpikes)' < 0 ) );
                    %NOTE: CURRENTLY ONLY FINDS FIRST INSTANCE, SINCE NOT RECURSIVE
                disp(['Apparent phase loss frameSpike/s: #',num2str(phaseLossSpike),...
                    ' (~',num2str(floor( inferTimes( frameLOCS( phaseLossSpike) ) )),'s)'])
                phaseLossSpikeInd = frameLOCS( phaseLossSpike );

                %Plot
                figure
                for phasmo = 1:size( phaseLossSpikeInd,2 )
                    subplot( 1, size( phaseLossSpikeInd,2 ), phasmo )
                    %plot( [inferTimes( frameLOCS(phaseLossSpike-2) : frameLOCS(phaseLossSpike+2) )], [syncStruct.AI.FrameSpike( frameLOCS(phaseLossSpike-2) : frameLOCS(phaseLossSpike+2) )] )
                    plot( [inferTimes( frameLOCS(phaseLossSpike(phasmo)-2) : frameLOCS(phaseLossSpike(phasmo)+2) )],...
                        [syncStruct.AI.FrameSpike( frameLOCS(phaseLossSpike(phasmo)-2) : frameLOCS(phaseLossSpike(phasmo)+2) )] )
                    hold on
                    %scatter( [inferTimes( frameLOCS(phaseLossSpike-2:phaseLossSpike+2) )], [syncStruct.AI.FrameSpike( frameLOCS(phaseLossSpike-2:phaseLossSpike+2) )] )
                    %scatter( [inferTimes( frameLOCS(phaseLossSpike) )], [syncStruct.AI.FrameSpike( frameLOCS(phaseLossSpike) )]+0.1, 'Color', [1,0,0] )
                    scatter( [inferTimes( frameLOCS(phaseLossSpike(phasmo)-2:phaseLossSpike(phasmo)+2) )], [syncStruct.AI.FrameSpike( frameLOCS(phaseLossSpike(phasmo)-2:phaseLossSpike(phasmo)+2) )] )
                    scatter( [inferTimes( frameLOCS(phaseLossSpike(phasmo)) )], [syncStruct.AI.FrameSpike( frameLOCS(phaseLossSpike(phasmo)) )]+0.1, 'Color', [1,0,0] )
                    %title(['Phase loss frameSpike (#',num2str(phaseLossSpike),') location'])
                    title(['Phase loss frameSpike (#',num2str(phaseLossSpike(phasmo)),') location'])
                    xlabel(['Time (s)'])
                end

                %Check if phase loss spike within believable proximity to intended 'true' position
                if itLOCS(phaseLossSpike+1) - phaseLossSpikeInd < meanPhaseDiff+2*stdPhaseDiff && ...
                        itLOCS(phaseLossSpike+1) - frameLOCS( phaseLossSpike-1 ) > meanPhaseDiff+2*stdPhaseDiff %"phase loss spike within normal distance of next iterator change AND n-1 spike not within believable range"
                    line([inferTimes(frameLOCS( phaseLossSpike - 1 )),inferTimes(frameLOCS( phaseLossSpike - 1 ))], [0.9,1.1],'LineStyle', ':', 'Color','k', 'LineWidth', 1.5) %Amend plot to show which element being removed
                    frameLOCS( phaseLossSpike - 1 ) = []; %Remove n-1 spike
                        %EXTRAORDINARILY EMPIRICAL
                        %This is based on assumption that two framespikes were detected in aberrantly close proximity OR framespike aberrantly long due to daq write stutter (e.g. 14May Fly 2 Exp 3 framespike #8820)
                        %Note: Actual element removal also heavily based on phase of diff with regards to length-1 nature etc
                    disp(['Framespike #', num2str(phaseLossSpike - 1), ' removed from data to restore phase'])
                    nonLinearityPresent = 1;
                else
                    ['-# Cannot (currently) rectify #-']
                end
            end

            %Recheck for phase loss
            temp = []; 
            for i = 1:nanmin( [size(frameLOCS,2), size( itLOCS,1)] )
                temp(i) = frameLOCS(i) - itLOCS(i);
            end
            if any(temp < 0)
                ['## Alert: Phase loss persists ##']
                crash = yes
            end

            %To add: More rectification methods depending on different data fuckery types




            %--------------------

            nSpikes = size(frameLOCS,2);
            disp(['Found ',num2str(nSpikes),' apparent peaks in DAQ framespike data'])
            disp(['(Approximate peak/presumed flip rate: ',...
                num2str(1 / ( nanmean( diff( frameLOCS ) ) / sampRate )),'Hz w/ ~',...
                num2str(nanstd( diff( frameLOCS ) ) / sampRate),'s SD)'])
            if hasPTB
                disp(['(Expected fliprate: ',num2str(matParamStruct.matSave.frequency),'Hz)'])
            end

            %Report on instability in framespike locations
            temp = diff( frameLOCS );
            disp(['Mean inter-spike interval: ', num2str(nanmean(temp)), ', +- ', num2str(nanstd( temp )), '( ',num2str( (nanstd( temp ) / nanmean(temp))*100 ),'% variability)'])

            %Interim QA
            if nSpikes > targetINum
                ['## Alert: Detected number of framespikes exceeds target ##'] %This should never happen
                crash = yes
            end

            %QA
            if nSpikes ~= lastIt
                disp(['-# Caution: nSpikes (',num2str(nSpikes),') and last iterator value (',num2str(lastIt),') differ #-']) %May actually be normal condition
                %crash = yes %Probably most likely if weirdness with iterator lasting for 1 frame/etc
            end

            %disp(['Number of detected framespikes: ',num2str(nSpikes),'; Reported last valid i value: ',num2str(btData(end,5)-1)])
            %Match framespikes to btData
            if ~batteryDesign
                disp(['Number of detected framespikes: ',num2str(nSpikes),'; Target based on btData (Script ver ',num2str(progVerNum),'): ',num2str(targetINum),' or ',num2str(targetINum-1)])
                    %-1 possibility comes from situations where iterator halfway through an element when script ended, thus leaving iterator at +1 but framespike not having happened yet
            %else
            %    disp(['Number of detected framespikes: ',num2str(nSpikes),'; Target based on btData (Script ver ',num2str(progVerNum),'): ',num2str(targetINum),' or ',num2str(targetINum-1)])
            end
            %if (batteryDesign == 0 o&& nSpikes == btData(end,5)-1 ) || ( batteryDesign == 1 )%&& (nSpikes == btData(end,5)) ) %-1 based on exactly N=1 testing
            %if (batteryDesign == 0 && nSpikes == targetINum ) || ( batteryDesign == 1 )%&& (nSpikes == btData(end,5)) ) %-1 based on exactly N=1 testing
            if (batteryDesign == 0 && (nSpikes == targetINum || nSpikes == targetINum-1) ) || ( batteryDesign == 1 )
                    %Note: As of v8.2, it may be the norm for btData at end to match nSpikes without - 1
                if nSpikes == targetINum%btData(end,5)-1
                    %disp(['Perfect match between number of detected framespikes (',num2str(nSpikes),') and reported last valid i value (',num2str(btData(end,5)-1),')'])
                    disp(['Perfect match between number of detected framespikes and reported last valid i value']) %Note: This may actually be bad for flipOnsetIndices below
                %else
                %    disp(['Number of detected framespikes: ',num2str(nSpikes),'; Reported last valid i value: ',num2str(btData(end,5)-1)])
                elseif nSpikes == targetINum - 1
                    disp(['n-1 match between number of detected framespikes and reported last valid i value'])
                end
            
                %Standard calcs
                if ~batteryDesign
                    %temp = bwlabel( btData(:,9) == 2 ); %Use onOff == 2 to find flips; Assumption of having btData, but is this analysis possible without it?
                    temp = bwlabel( [0; diff( btData(:,5) )] ); %Use changes in i to find 'flips'
                        %Note: This value may technically be always 1 less than actual stimuli presented, since current DAQ writing is only done at a change, not at initialisation
                    [~,flipOnsetIndices] = ismember( [1:nanmax(temp)], temp ); %First position of state switch in onOff
                        %Necessary to find first position in case of different duty cycles?
                        %Reminder: FrameSpike seems to occur at position of switch to next (e.g. FrameSpike #1 occurs at i transition from 1 to 2 in btData)
                            %(This is also where the iterator changes, although note that the iterator value represents the 'true' i while the framespike count will be 1 behind)
                                %((As in, when i changes to 33, the iterator will almost immediately represent 33 but framespike 33 will only occur at the end of the i event))

                    %Report
                    disp([num2str(length(flipOnsetIndices)),' changes in i detected in btData'])

                    %QA
                    %if ( batteryDesign ~= 1 && length( flipOnsetIndices ) ~= nSpikes ) || ( batteryDesign == 1 && length( flipOnsetIndices ) ~= nSpikes-1 ) %HIGHLY EMPIRICAL
                    if length( flipOnsetIndices ) ~= btData(end,5)-1 %Probably very unlikely; Potentially redundant QA
                        ['## ERROR: CRITICAL DESYNC BETWEEN flipOnsetIndices AND btData ##']
                        crash = yes
                    end
                    %More QAs
                        %These QAs currently disabled until correctness identified
                    %{
                    if progVerNum >= 8.45
                        if length( flipOnsetIndices )+1 ~= targetINum
                            ['## Alert: Disparity between target final i value (',num2str(targetINum),') and detected flip positions (',num2str(length( flipOnsetIndices )+1),') ##']
                        end
                        if length( flipOnsetIndices )+1 ~= itProcRangeIndUnwrapped(end)
                            ['## Alert: Disparity between final reported i value (',num2str(lastIt),') and detected flip positions (',num2str(length( flipOnsetIndices )+1),') ##']
                        end
                    end
                    if length( flipOnsetIndices ) ~= nSpikes
                        ['## Alert: Disparity between detected framespikes and detected flip positions ##']
                        crash = yes
                    end
                    %}
                    %chester

                    %Establish which framespike to match with last element of flipOnsetIndices
                    %Pre QA
                    if nanmin( [lastIt, length(flipOnsetIndices),  length(frameLOCS), length(itLOCS)] ) - nanmax( [lastIt, length(flipOnsetIndices),  length(frameLOCS), length(itLOCS)] ) > 2
                        ['-# Alert: Non-insignificant disparity between smallest and largest iterator/frame/etc counts #-']
                        crash = yes %Might be 'normal' cases where this occurs
                    end
                    lastUseful = [];
                    %---------
                    %Empirical booleans
                        %"The Empire Strikes Back"
                    if lastIt == length(flipOnsetIndices)+1 && frameLOCS(end) > itLOCS(end)
                            %" Ends with a full iterator block (inc. framespike) immediately followed by an iterator increase BUT NOT a change in btData "
                                %Note that this is actually really only the case for a PTB end happening coincidentally right on an i change (Since btData saving occurs 1 loop cycle after iterator increase)
                        lastUseful = length(frameLOCS)-1;
                        disp(['Last flip-matching framespike calculated as #',num2str(lastUseful)])
                        disp(['(Terminal pattern Foxtrot 2 Beta 2 [frameSpike n - 1])']) %So named because it's how 14May F2 B2 ends
                    elseif lastIt == length(flipOnsetIndices) && frameLOCS(end) > itLOCS(end)
                            % "Ends midway through an iterator block, no framespike at very end, with a change in btData at last framespike"
                            % Presumable common case?
                        lastUseful = length(frameLOCS);
                        disp(['Last flip-matching framespike calculated as #',num2str(lastUseful)])
                        disp(['(Terminal pattern Foxtrot 2 Beta 3 [frameSpike n])']) %So named because it's how 14May F2 B3 ends
                    elseif lastIt == length(flipOnsetIndices)+1 && frameLOCS(end) < itLOCS(end)
                            % Similar to F2B3, except iterator larger(?)
                        lastUseful = length(frameLOCS);
                        disp(['Last flip-matching framespike calculated as #',num2str(lastUseful)])
                        disp(['(Terminal pattern Foxtrot 1 Beta 1 [frameSpike n])']) %So named because it's how 14May F1 B1 ends
                    else
                        ['-# Uncommon ending case detected; Specification does not exist #-']
                        todo = yes
                    end
                    %---------
                    %mushroom
                    disp(char(10))
                    %QA for pre-drift between inferTimes and btData
                    preDrift = (inferTimes( frameLOCS(lastUseful) ) - inferTimes( frameLOCS(1) )) - (btData( flipOnsetIndices(lastUseful), 6) - btData( flipOnsetIndices(1), 6));
                    disp(['Initial apparent drift magnitude: ',num2str(preDrift),'s']) %"Tokyo Drift"
                    if abs(preDrift) > 2*nanmean( diff( btData( flipOnsetIndices, 6) ) ) %Unlikely to be possible if post-hoc inferTimes correction deployed
                        ['## Caution: Apparent starting drift (',num2str(preDrift),') between inferTimes and btData exceeds 2*mean inter-framespike interval (',num2str(2*nanmean( diff( btData( flipOnsetIndices, 6) ) )),') ##']
                        %Empirically this may occur 'naturally' for a dataset exceeding 1h
                    end

                    %SUPER EXPERIMENTAL METHOD TO POST-HOC CORRECT inferTimes
                        %Note: Absolutely critically dependent on lastUseful being correct
                    if options.postHocCorrectInferTimes
                        corrFa = (inferTimes( frameLOCS(lastUseful) ) - inferTimes( frameLOCS(1) )) - (btData( flipOnsetIndices(lastUseful), 6) - btData( flipOnsetIndices(1), 6));
                        corrFa = corrFa + 0.013; %Utterly empirical additional tweak based on N=1 (14May F2B3)

                        %Reporter plot
                            %Basically plots drift between inferTimes and btData
                                %Note that theoretically btData shouldn't 'drift', whereas inferTimes, as an interpolated set of values may
                        blirg = nan( lastUseful, 1 );
                        blorg = nan( lastUseful, 1 );
                        iferTi = inferTimes;
                        iferTi = iferTi * ( iferTi(end) / (iferTi(end)+corrFa) );
                        %iferTi = iferTi / (1+cur.p1+cur.p2);
                        for i = 1:lastUseful
                            blirg(i) = (iferTi( frameLOCS(i) ) - iferTi( frameLOCS(1) )) - (btData( flipOnsetIndices(i), 6) - btData( flipOnsetIndices(1), 6));
                            blorg(i) = (inferTimes( frameLOCS(i) ) - inferTimes( frameLOCS(1) )) - (btData( flipOnsetIndices(i), 6) - btData( flipOnsetIndices(1), 6));
                        end
                        figure
                        plot(blorg)
                        hold on
                        plot(blirg)
                        title(['Uncorrected and corrected inferTimes against btData time'])
                        %If the orange line isn't flat-ish then there are problems

                        inferTimes = inferTimes * ( inferTimes(end) / (inferTimes(end)+corrFa) );
                        estimatedPTBEndTime = imStartTime + inferTimes( bleachEndInd ); %Recalculate
                        
                        disp(['inferTimes (, estimatedPTBEndTime) drift corrected by ',num2str(corrFa),'s based on btData'])
                        
                    end

                    %QA for drift between inferTimes and btData
                    %maomao
                    finalDrift = (inferTimes( frameLOCS(lastUseful) ) - inferTimes( frameLOCS(1) )) - (btData( flipOnsetIndices(lastUseful), 6) - btData( flipOnsetIndices(1), 6));
                    disp(['Final apparent drift magnitude: ',num2str(finalDrift),'s'])
                    if abs(finalDrift) > nanmean( diff( btData( flipOnsetIndices, 6) ) ) %Unlikely to be possible if post-hoc inferTimes correction deployed
                        ['## Alert: Apparent drift between inferTimes and btData exceeds mean inter-framespike interval ##']
                        crash = yes
                    end

                    estPTBInferTimeStart = inferTimes( frameLOCS(1) ) - btData( flipOnsetIndices(1), 6); %Theoretical time (inferTimes reference) PTB started at
                        %Use framespike to find first flip position in inferTimes, then subtract known duration since PTB start from that
                        %Note: Only valid as long as pre-loop framespike being dropped
                    %%estPTBInferTimeEnd = inferTimes( frameLOCS(end) ) + ( btData(end,6) - btData( flipOnsetIndices(end), 6) 
                    estPTBInferTimeEnd = inferTimes( frameLOCS(lastUseful) ) + ( btData(end,6) - btData( flipOnsetIndices(end), 6) );  
                        %Similarly, find inferTime[s] where last framespike happened, then add time distance (btData self-report) between that and last btData element

                    %Add more code
                    %if progVerNum > 8.55
                    %   ['to do : add timing improvement with new timeSinceLastFlip column in btData']
                    %  %to do = yes
                    %end
                    
                    %Report if necessary
                    if nonLinearityPresent
                        ['-# Caution: Non-linearity now present in data without appropriate countermeasures #-']
                        %Refers to fact that removal of frameLOCS element without also slicing data is likely to lead to 1/inter-framespike interval timing issues compared to 'true' timing
                        %This can probably be fixed by adjustments to frame interpolation below
                    end

                else

                    preTemp =  [0; diff( btData(:,5) )];
                    temp =  zeros( size(btData,1), 1 );
                    temp( find(btData(:,5) == 1, 1, 'first') ) = 1;
                    temp( preTemp == 1 ) = btData( find(preTemp == 1), 5 ); %Don't use bwlabel, because bwlabel doesn't quite work if there are no repetitions of elements (i.e. Speedy optomotor)
                        %Related: This may encounter issues if stimuli switching faster than BT data being saved
                    [~,flipOnsetIndices] = ismember( [1:nanmax(temp)], temp );
                    
                    %g11

                    %Remove opto-related i elements
                    if matParamStruct.matSave.optoElements > 0
                        flipOnsetIndices( 1:matParamStruct.matSave.optoElements ) = [];
                        disp(['Battery: Proximal ', num2str(matParamStruct.matSave.optoElements),' btData flip detects removed due to opto relation'])
                    end

                    %Check for approximate timing between first flipOnset and first non-opto iteration (Should work even if no opto)
                    if abs( btData( flipOnsetIndices(1) , 6 ) - inferTimes( itLOCS(1) ) ) > 5 %Allow 5 startup PTB time effectively
                        ['-# Alert: Potential error in first (DAQ) iteration <--> first BT flip detection timing #-'] %Probably most likely if first itLOCS is at start of PTB whilst flipOnset is post-opto, or vice versa
                        if exist( 'omittedFrameLOCPK' )
                            ['First btData time: ',num2str(btData(1,6))]
                            ['Initialisation framespike inferTime: ',num2str(inferTimes( omittedFrameLOCPK(1) ))]
                        end
                        if ~batteryDesign
                            crash = yes
                        else
                            ['-# Allowing to proceed because battery design #-']
                        end
                    end

                    %New
                    firstUseful = NaN;
                    lastUseful = NaN;

                    firstUseful = 1; %Note: Assumes that initialisation framespike/etc have been cleared
                    disp(['First flip-matching framespike assumed as #',num2str(firstUseful)])

                    %---------
                    %Empirical booleans
                        %"The Empire Strikes Back"
                    if lastIt == length(flipOnsetIndices)+1 && frameLOCS(end) > itLOCS(end)
                            %" Ends with a full iterator block (inc. framespike) immediately followed by an iterator increase BUT NOT a change in btData "
                                %Note that this is actually really only the case for a PTB end happening coincidentally right on an i change (Since btData saving occurs 1 loop cycle after iterator increase)
                        lastUseful = length(frameLOCS)-1;
                        disp(['Last flip-matching framespike calculated as #',num2str(lastUseful)])
                        disp(['(Terminal pattern Foxtrot 2 Beta 2 [frameSpike n - 1])']) %So named because it's how 14May F2 B2 ends
                    elseif lastIt == length(flipOnsetIndices) && frameLOCS(end) > itLOCS(end)
                            % "Ends midway through an iterator block, no framespike at very end, with a change in btData at last framespike"
                            % Presumable common case?
                        lastUseful = length(frameLOCS);
                        disp(['Last flip-matching framespike calculated as #',num2str(lastUseful)])
                        disp(['(Terminal pattern Foxtrot 2 Beta 3 [frameSpike n])']) %So named because it's how 14May F2 B3 ends
                    elseif lastIt == length(flipOnsetIndices)+1 && frameLOCS(end) < itLOCS(end)
                            % Similar to F2B3, except iterator larger(?)
                        lastUseful = length(frameLOCS);
                        disp(['Last flip-matching framespike calculated as #',num2str(lastUseful)])
                        disp(['(Terminal pattern Foxtrot 1 Beta 1 [frameSpike n])']) %So named because it's how 14May F1 B1 ends
                    else
                        ['-# Uncommon ending case detected; Specification does not exist #-']
                        todo = yes
                    end
                    %---------
                    estPTBInferTimeStart = inferTimes( frameLOCS(firstUseful) ) - btData( flipOnsetIndices(1), 6); %Theoretical time (inferTimes reference) PTB started at
                    estPTBInferTimeEnd = inferTimes( frameLOCS(lastUseful) ) + ( btData(end,6) - btData( flipOnsetIndices(end), 6) );  

                    %Old
                    %{
                    %if exist( 'progIdent' ) && ( ~isequal( progIdent, 'FLY_SEQUENTIAL_DEPENDENCIES_v8dot3_XM' ) && ~isequal( progIdent, 'FLY_SEQUENTIAL_DEPENDENCIES_v8dot4_XM' ) ) 
                    if progVerNum >= 8.45
                        estPTBInferTimeStart = inferTimes( frameLOCS(1) ) - btData( 1 , 6); %Assume that first framespike occurred immediately prior to first BT element
                            %Note that this is only valid as long as first element of frameLOCS is referring to that spike
                        disp(['Using legacy assumptions for battery first framespike position'])
                    else
                        estPTBInferTimeStart = inferTimes( frameLOCS(1) ) - btData( flipOnsetIndices( matParamStruct.matSave.optoElements+1 ), 6); %For older (battery) data, assume that first framespike occurs immediately after optomotor
                    end
                    estPTBInferTimeEnd = inferTimes( frameLOCS(end) ) + ( btData(end,6) - btData( flipOnsetIndices(end), 6) );  %This is generally a safe assumption, regardless of framespike position
                    %}
                    %QA
                    if estPTBInferTimeStart > 10
                        ['## Alert: PTB estimated to have started >10s (',num2str(estPTBInferTimeStart),'s) after TS initiation; Error? ##'] %Likely to either be a very lagged PTB start or an error in framespike/frameLOCS attribution
                        if ~batteryDesign
                            crash = yes
                        else
                            disp(['(Again; Allowing to proceed because battery)'])
                        end
                    end

                end

                %estPTBInferTimeStart = inferTimes( frameLOCS(1) ) - btData( flipOnsetIndices(1), 6); %Theoretical time (inferTimes reference) PTB started at; Moved above for separation
                    %Use framespike to find first flip position in inferTimes, then subtract known duration since PTB start from that
                [inferPTBStartDisp, inferPTBStartInd] = min( abs(inferTimes - estPTBInferTimeStart) );  %Find closest match in inferTimes for estimated PTB start time based on framespike
                    %Returns respectively the time error between the inferTimes point and PTB start, as well as the index of inferTimes (/TS) where PTB started (But not necessarily first stim)
                 
                %estPTBInferTimeEnd = inferTimes( frameLOCS(end) ) + ( btData(end,6) - btData( flipOnsetIndices(end), 6) ); 
                    %Similarly, find inferTime[s] where last framespike happened, then add time distance (btData self-report) between that and last btData element
                [inferPTBEndDisp, inferPTBEndInd] = min( abs(inferTimes - estPTBInferTimeEnd) );

                %QA
                if abs(inferPTBStartDisp) > 5 || abs(inferPTBEndDisp) > 5 %Check if closest match >5s (Likely indicates that PTB ran before/after TS)
                    ['-# Alert: PTB apparent start/end time significantly (>5s) precedes/postcedes TS start/end time ##']
                    crash = yes
                end

            else %Note: Not tested with v8.2 new position of Data_Array saving
                ['I M P E R F E C T I O N']
                crash = yes

                %If this actually happens, add code to allow for (hopefully) minor disparities
            end
            
        end

        %------------------------------
    
    
        %% Dynamic calcs
        guessIndex = [{'PTB'},{'Exp*.xml'},{'Framespike + PTB'}]; %List of ways to guess the start time
        if ~isShortcutting
            guessMode = 0;
            if hasPTB
                %if ~isShortcutting
                if ~hasDaqData
                    bestGuessCommenceTime = ptbEndTime - inferTimes( bleachEndInd ); %Use PTB end as a known point in posix, match with inferred relative TS time
                    guessMode = 1;
                else
                    bestGuessCommenceTime = ( btData( flipOnsetIndices(1), 2) / 1000 ) - inferTimes( frameLOCS(1) ); %Subtract inferTimes at 1st framespike from known PTB posix at same point
                    guessMode = 3;
                end
                %else
                %bestGuessCommenceTime = shortStruct.bestGuessCommenceTime;
                %end
                %guessMode = 1;
                %disp(['Best guess commencement time (from PTB): ',datestr(datetime( bestGuessCommenceTime , 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10'))])
            else
                bestGuessCommenceTime = imStartTime; %Use self-reported start time blindly
                guessMode = 2;
                %disp(['Best guess commencement time (from Exp*.xml): ',datestr(datetime( bestGuessCommenceTime , 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10'))])
            end
        else
            bestGuessCommenceTime = shortStruct.bestGuessCommenceTime;
            guessMode = shortStruct.guessMode;
        end
        disp([char(10),'Best guess commencement time (from ',guessIndex{guessMode},'): ', ...
            datestr(datetime( bestGuessCommenceTime , 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10'))])
                %Note: Relates to estimated experiment start time, not just PTB or imaging
    
        if ~isShortcutting %Only calculate if actually loaded h5 data
            %inferTimesPosix = inferTimes + bestGuessCommenceTime; %Posix form of inferTimes, for simplicity; Moved below
        
            %Do these here so shortcutting can happen
            tsLastFrameTime = inferTimes( frameEndInd );
            tsLastArduinoTime = inferTimes( bleachEndInd );
            tsLastTime = inferTimes( end );
            if hasPTB %&& hasData 
                if ~hasDaqData
                    inferTimesPosix = inferTimes + bestGuessCommenceTime; %Posix form of inferTimes, for simplicity
                    lastImStimFrameInd = find( inferTimesPosix - ptbEndTime >= 0 , 1, 'first' ) - 1; %The last frame to receive a full stimulation; PTB posix inference
                    firstImStimFrameInd = find( inferTimesPosix - ptbStartTime >= 0 , 1, 'first' ); %The first imaging frame to receive stimulation
                else
                    lastImStimFrameInd = inferPTBEndInd;
                    firstImStimFrameInd = inferPTBStartInd;
                end
            end
        else
            tsLastFrameTime = shortStruct.tsLastFrameTime;
            tsLastArduinoTime = shortStruct.tsLastArduinoTime;
            tsLastTime = shortStruct.tsLastTime;
            if hasPTB %&& hasData 
                lastImStimFrameInd = shortStruct.lastImStimFrameInd;
                firstImStimFrameInd = shortStruct.firstImStimFrameInd;
            end
            sampRate = shortStruct.sampRate;
        end
    
        %here
    
        %% Report
    
        if hasPTB
            disp(['Known PTB start time: ', datestr(datetime( ptbStartTime, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ))])
            disp(['Known PTB end time: ', datestr(datetime( ptbEndTime, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ))])
            disp(['PTB self-calculated duration: ',...
                num2str( (ptbEndTime - ptbStartTime)/60 ),' mins (',...
                num2str( (ptbEndTime - ptbStartTime) ),'/',num2str( (ptbEndTime - ptbStartTime)+btData(1,1)),' uncorrected/corrected s)'])
                %Corrected s comes from adding PTB spinup time (First btData row) to final duration
            disp(['Imaging - PTB start time disparity: ',num2str( imStartTime - ptbStartTime ),'s'])
            %Negative means imaging started first, which is typical situation
                %See above comment re: Imaging start time
    
            %{
                %Currently disabled on account of questions over btData column 1 identity
            ptbEstimatedArdStartTime = ( btData(1,2)/1000 ) - btData(1,1);
                %Estimated accuracy: Unknown (Maybe +- repThresh duration? Also, depends on placement of startTime in PTB script)
            disp([ 'PTB estimated arduino commencement time: ',...
                datestr(datetime( ptbEstimatedArdStartTime, 'ConvertFrom', 'posixtime', 'Timezone', 'UTC+10')) ])
            %}
        end
    
        disp( ['ThorSync-calculated imaging duration: ',...
        num2str( (tsLastFrameTime/60) ),' mins (',...
        num2str( (tsLastFrameTime) ),' s)'] )  
        %disp(['Estimated imaging end time (Exp self report -> Last frame high): ', ...
        %    datestr( datetime( imStartTime + (inferTimes( frameEndInd )), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])
        disp(['Estimated imaging end time (',guessIndex{guessMode},' -> Last frame high): ', ...
            datestr( datetime( bestGuessCommenceTime + tsLastFrameTime, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])
    
        disp([ 'Estimated arduino duration: ', num2str((tsLastArduinoTime)/60),'m' ])
        %disp(['Estimated arduino end time (Exp self report -> Last bleach high): ', ...
        %    datestr( datetime( imStartTime + (inferTimes( bleachEndInd )), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])
        disp(['Estimated arduino end time (',guessIndex{guessMode},' -> Last bleach high): ', ...
            datestr( datetime( bestGuessCommenceTime + tsLastArduinoTime, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])
    
        disp(['Estimated TS end time: ',...
            datestr( datetime( bestGuessCommenceTime + tsLastTime, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])
    
        %% Plot
    
        if doPlot && ~isShortcutting
            %Plot
            %datList = [{'AI'},{'CI'},{'DI'},{'Freq'},{'Global'}]; %Could probably be dynamicised
            datList = fieldnames( syncStruct )'; 
    
            for dat = datList
                thisFieldNames = fieldnames(syncStruct.(dat{1}));
    
                figure
                for row = 1:size(thisFieldNames,1)
                    thisFiel = thisFieldNames{row};
    
                    subplot( size(thisFieldNames,1) , 1, row )
                    plot( syncStruct.(dat{1}).(thisFiel) )
                    title(thisFiel)
                end
    
            end
        end
    
        %% Load actual data
        hasData = 0;
        if dataSource == -1
            thisImData = thisBlock.greenChannel;
            disp(['Function call data loaded'])
            disp(['Imaging size: ', num2str(size(thisImData))])
            hasData = 1;
        else
            imDataFile = [strcat( dataFolder,filesep,expName,filesep,'avg_z_green_aligned.mat' )];
            %QA
            if isempty(dir(imDataFile))
                ['## Imaging data not found! ##']
                %crash = yes
                hasData = 0;
            else
                hasData = 1;
                load( imDataFile );
                disp(['External data loaded'])
                thisImData = avg_z_green_aligned;
                clear avg_z_green_aligned
                %disp(['Imaging size: ', num2str(size(avg_z_green_aligned))])
                disp(['Imaging size: ', num2str(size(thisImData))])
            end
            %load( imDataFile );
            %disp(['Imaging size: ', num2str(size(avg_z_green_aligned))])
        end
    
        %% Derive more metrics
    
        %PTB sequence
            %Moved above
        %{
        if hasPTB
            if nansum( isnan( btData(end,:) ) ) == numel( btData(end,:) )
                btData(end,:) = [];
            end
            %QA
            if nansum( isnan(btData), 'all' ) > 0
                ['-# Warning: NaN elements found in btData #-']
                crash = yes %Again, not critical but problematic if not at end
            end
            randomSequence = matParamStruct.matSave.randomSequence( unique(btData(:,5)) );
                %Note: Not tethered to PTB/Imaging/etc reference frame
        end
        %}
    
        if hasData
            %Imaging frames
            if ~isShortcutting
                bwFrameOut = bwlabel( syncStruct.DI.frameOutData );
                [~,frameOnsetIndices] = ismember( [1:nanmax(bwFrameOut)], bwFrameOut );
                tsImFrameCount = numel(frameOnsetIndices);
                %Calculate best guess time for each frame
                imFrameTimes = inferTimes( frameOnsetIndices );
            else
                frameOnsetIndices = shortStruct.frameOnsetIndices;
                tsImFrameCount = shortStruct.tsImFrameCount;
                imFrameTimes = shortStruct.imFrameTimes; %If this crashes, delete SHORT files and regenerate
            end
            disp(['ThorSync calculated number of imaging frames: ', num2str(tsImFrameCount)])
            %disp(['Expected number of frames: ',num2str(size( avg_z_green_aligned , 3 ) * (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames))])
            disp(['Expected number of frames: ',num2str(size( thisImData , 3 ) * (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames))])
            %QA
            if  tsImFrameCount ~= size( thisImData , 3 ) * (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames) %Note: All of these vars need to be manually imported/defined
                ['-# Alert: TS framecount and imaging framecount differ #-']
                crash = yes %Eventually will probably have systems to handle this case
            end
        end
    
        %% Some more plots
    
        if hasData && doPlot% && ~isShortcutting
            %Imaging frame timing instability (courtesy of TS)
                %If there are many unique values (i.e. >2) on this graph, timing was unstable 
            figure
            hist( diff(frameOnsetIndices)/sampRate, 128 )
            title(['Inter-frame interval hist'])
            xlabel(['Time (s)'])
            ylabel(['Count'])
        end

        %% Photodiode stuff
        %Quick allocation for battery design currently
        if ~isShortcutting && hasPhotData %&& batteryDesign 
        %Incorporate photodiode data
                %NOTE: INEFFICIENT; PROBABLY REPLACE WITH INTERPOLATED/SMALLER VERSION SOON
            %if hasPhotData
            photData = syncStruct.AI.Photodiode;
            %end
        elseif isShortcutting && hasPhotData %&& batteryDesign 
            photData = shortStruct.photData;
            if isfield( shortStruct, 'photProc' )
                photProc = shortStruct.photProc;
            end
        end

        if ~isShortcutting && hasPhotData && cleanPhotData
            %destroyer
            disp([char(10),'Cleaning photodiode data'])

            %gachiakuta

            %Advanced phot rectification
                %Non-functional
            %{
            temp = -photData;

            temp(1:69418) = nan;
            temp(51441900:end) = nan;


            minPeriodGap = (( numel( matParamStruct.matSave.daqRange ) / matParamStruct.matSave.frequency ) * sampRate ) * 0.9;


            [PKS,LOCS,W] = findpeaks( temp , 'MinPeakDistance', minPeriodGap );

            temprange = [ nanmean( mink( temp, 15 ) ) , nanmean( maxk( temp, 15 ) ) ];


            %Saw
            blirg = -sawtooth( [0:1/minPeriodGap*5.2:pi*numel(PKS)] ); %5.2 arb?

            blerg = normalize( blirg , 'Range', temprange );

            templerg = temp - blerg;

            %Or
            blorg = nan( 1, numel(temp) );
            for bInd = 2:size(LOCS,2)
                coords = LOCS(bInd-1) : LOCS(bInd);
                %blorg( coords ) = -0.01*bInd;
                %recfunc = linspace( temprange(2) , temprange(1), numel(coords) );
                cutcoords = floor( linspace( nanmin(coords), nanmax(coords), 100 ) );
                cutdata = temp( cutcoords );
                p = polyfit( cutcoords, cutdata, 2 );
                f = polyval(p,coords);
                
                blorg(coords) = f;
            end

            plot( )








            %}
            
            photTemp = photData;
            photTemp( photTemp < 0.3 ) = 0; %Hopefully remove all iterator shared signal (Note: Will delete low luminosity true phot events)
            photTemp( photTemp >= 0.3 ) = 1;

            %Smoothing method
            smoothVal = 100;
            photTemp = smooth(photTemp, smoothVal ); %Empirical numbers
                %Note: This is a bit unconventional, smoothing a binary trace
            photTemp( photTemp ~= 0 ) = 1;
            photTemp( photTemp ~= 1 ) = 0;

            photLabel = bwlabel( photTemp ); %Find post-smoothed stimulus events

            %'Fix' artificial lead/lag induced by smoothing
            disp(['Fixing artificial lead/lag in phot data'])
            tic
            flatFails = 0;
            for i = 1:nanmax( photLabel )
                if nansum( photLabel == i ) > smoothVal 
                    startEnd = [ find( photLabel == i , 1, 'first' ), find( photLabel == i , 1, 'last' ) ]; %Find putative start/end of this stimulus event
                    coords = [ startEnd(1):startEnd(1)+floor(0.5*smoothVal), startEnd(2)-floor(0.5*smoothVal):startEnd(2) ]; %Identify artificial lead/lag coords
                    coords( coords < 1 ) = []; coords( coords > length(photTemp) ) = []; %Clean
                    photTemp( coords ) = 0; %Flatten
                else
                    flatFails = flatFails + 1;
                    if ~batteryDesign
                        ['-# error: cannot flatten smoothing lead/lag for apparent stim. event #',num2str(i),' due to small size (',num2str(nansum( photLabel == i )),') #-']
                    end
                end
            end
            toc
            if batteryDesign && flatFails > 0
                ['-# battery design encountered ',num2str(flatFails),' failures to flatten lead/lag (',num2str(nanmax( photLabel )),' total detected events)']
            end

            %darkness

            %Gap ligation method on original photodiode data
                %Deprecated because v. slow
            %{  
            %Find min. valid gap
            if hasPTB
                if ~batteryDesign
                    %maxFreq = matParamStruct.matSave.frequency; %Use stimuli to decide gap
                    maxFreq = matParamStruct.matSave.panelFrequency; %Use panel frequency to decide gap
                else
                    not coded yet; %Ask andre
                end
            else
                %maxFreq = 6; %Assume default 6Hz
                maxFreq = 60; %Assume default 60Hz
            end
            minValidGap = sampRate / maxFreq;
    
            %Label gaps
            invPhotLabel = bwlabel( ~photTemp );
            crimson
            %For loop ligate gaps
                %VERY SLOW (Like, 1h per block)
            tic
            ligatedPhotGaps = zeros( 1, nanmax(invPhotLabel) );
            parfor i = 1:1000%nanmax( invPhotLabel ) 
                if nansum( invPhotLabel == i ) < 1.1*minValidGap
                    %disp([num2str(nansum( invPhotLabel == i ))])
                    %photTemp( invPhotLabel == i ) = 1;
                    ligatedPhotGaps(i) = 1;
                end
            end
            disp([num2str(nansum(ligatedPhotGaps)),' gaps in phot trace ligated (',num2str(toc),'s)'])
            %}

            %Prepare for wider use
            photProc = photTemp;

            %Some potentially useful post-hoc reporting
            temp = bwlabel( photProc );
            disp(['Phot indicates ',num2str(nanmax(temp)),' stimulus events'])            

        end
    
            %% Save data, to increase speed for subsequent runs
        if saveShortcut == 1 && isShortcutting ~= 1
            shortStruct = struct;
    
            temp = strsplit( hFolder, filesep );
            shortStruct.h5Name = temp{end}; %A little more safe than just assuming expName is correct 
            
            shortStruct.sampRate = sampRate;
    
            %shortStruct.inferTimes = inferTimes;
            %shortStruct.inferTimesPosix = inferTimesPosix;
    
            shortStruct.estimatedPTBEndTime = estimatedPTBEndTime;
            shortStruct.bestGuessCommenceTime = bestGuessCommenceTime;
            shortStruct.guessMode = guessMode;
    
            shortStruct.tsImFrameCount = tsImFrameCount;
            shortStruct.frameOnsetIndices = frameOnsetIndices;
            shortStruct.imFrameTimes = imFrameTimes;
    
            shortStruct.tsLastFrameTime =tsLastFrameTime;
            shortStruct.tsLastArduinoTime = tsLastArduinoTime;
            shortStruct.tsLastTime = tsLastTime;
    
            if hasPTB %&& hasData 
                shortStruct.lastImStimFrameInd = lastImStimFrameInd;
                shortStruct.firstImStimFrameInd = firstImStimFrameInd;
            end

            shortStruct.hasDaqData = hasDaqData;
            shortStruct.hasPhotData = hasPhotData;
            if hasPhotData
                shortStruct.photData = photData;
                if exist('photProc')
                    shortStruct.photProc = photProc; %Yet larger shortcut files...
                end
            end

            if exist('arbPhaseShift')
                shortStruct.arbPhaseShift = arbPhaseShift;
            end
    
            %save( [hFolder,filesep,'h5Shortcut.mat'], 'shortStruct' )
            %Check/Make folder to put shortcut files into
            shortFolder = strcat( dataFolder,filesep,'SHORT' );
            if exist(shortFolder) ~= 7
                mkdir( shortFolder )
                disp(['Shortcut folder made at ',shortFolder])
            end
            save( [shortFolder,filesep,shortStruct.h5Name,'_h5Shortcut.mat'], 'shortStruct' )
            %sword
            disp(['Shortcut data saved'])
    
        end
        %fighti    
    
        %% Tether image to known reference frame
    
        %To do: Non-linear frame placement based on TS frame timings
            %Probably not a huge issue, but likely to be useful
        if hasData && hasPTB
            %lastImStimFrameInd = find( frameOnsetIndices - bleachEndInd >= 0 , 1, 'first' ) - 1; %The last frame to receive a full stimulation; Arduino low inference
            %%lastImStimFrameInd = find( inferTimesPosix - ptbEndTime >= 0 , 1, 'first' ) - 1; %The last frame to receive a full stimulation; PTB posix inference; Moved above
                %Note: May be at any point in a volume
            %%firstImStimFrameInd = find( inferTimesPosix - ptbStartTime >= 0 , 1, 'first' ); %The first imaging frame to receive stimulation; Moved above
                %Again, can be any point
                %Note: Both of these in TS reference frame
            imStimStart = find( frameOnsetIndices - firstImStimFrameInd > 0 , 1 , 'first' ); %Might theoretically be better to find min, rather than first after?
            imStimEnd = find( frameOnsetIndices - lastImStimFrameInd > 0 , 1 , 'first' );
            %QA
            if isempty( imStimStart ) || isempty( imStimEnd )
                if isempty( imStimStart )
                    ['-# Alert: Failure to identify imaging start #-']
                elseif isempty( imStimEnd )
                    ['-# Alert: Failure to identify imaging end #-']
                    if exist('inferTimes')
                        ['Apparent PTB duration: ', num2str((inferTimes( lastImStimFrameInd )/60) ),'m']
                    end
                end
                if isShortcutting ~= 1
                ['Imaging started ', num2str( inferTimes( frameOnsetIndices(1) ) ),'s after TS start']
                ['Imaging ended apparently after ', num2str( inferTimes( frameOnsetIndices(end) ) / 60 ),'m']
                ['PTB self-reported ended after ',num2str(btData(end,6)/60),'m']
                ['(TS ended after ', num2str( inferTimes( end ) / 60 ),'m)']
                end
                crash = yes
            end
            disp([ 'Stimulation comprised ', num2str( (imStimEnd - imStimStart) / numel( frameOnsetIndices )*100 ), '% of imaging duration (',...
                num2str(imStimEnd - imStimStart),' of ',num2str(numel( frameOnsetIndices )),' imaging frames)' ])
    
            %Make list of frames -> vols
            imFrameList = [1:tsImFrameCount]; %All frames, irrespective of Z (or imaging)
            %imFrameListZ = reshape(imFrameList, (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames), size( avg_z_green_aligned , 3 )); %Z in rows, Vol # in cols
            imFrameListZ = reshape(imFrameList, (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames), size( thisImData , 3 )); %Z in rows, Vol # in cols
            %g11
    
            %Find first and last usable volume
            %Start
            [tempX,tempY] = find( imFrameListZ == imStimStart );
            if tempX ~= 1
                adjImStimStart = imFrameListZ( 1 , tempY+1 ); %Frame was not first of vol, find start of next full vol
            else
                adjImStimStart = imStimStart; %Frame was first of vol (Lucky lucky)
            end
            %End
            %{
            [tempX,tempY] = find( imFrameListZ == imStimEnd );
            if tempX ~= 1
                adjImStimEnd = imFrameListZ( 1 , tempY+1 ); %Frame was not first of vol, find start of next full vol
            else
                adjImStimEnd = imStimEnd; %Frame was first of vol (Lucky lucky)
            end
            %}
            adjImStimEnd = imStimEnd; %Use last stim frame as end, regardless of vol position
                %This simplifies some sequence calcs
            %Express in volume reference
            [~,adjImStimStartVol] = find( imFrameListZ == adjImStimStart ); %Note that in 1 Z-plane situations these will actually be the same it seems
            [~,adjImStimEndVol] = find( imFrameListZ == adjImStimEnd ); %Note that in 1 Z-plane situations these will actually be the same it seems
            clear adjImStimEnd adjImStimStart %Clear these to remind of fact that they may be inaccurate if end volumes adjusted because NaNs etc
            disp(['Adjusted imaging volumes of stimulation interest: ', num2str(adjImStimStartVol),' : ', num2str(adjImStimEndVol),...
                ' (',num2str(adjImStimEndVol-adjImStimStartVol+1),' total)' ])
                %Note: For multi-Z plane imaging, this section may desync imStimStart/adjImStimStartVol/etc [Vol ref.] and firstImStimFrameInd [TS ref] by up to a vol's worth of time (Probably) 

            %Testatory figure of transient against phot (if existing)
            %{
            if hasPhotData
                dataCoords = [imStimStart:imStimEnd];
                photCoords = [firstImStimFrameInd:lastImStimFrameInd];
                blirg = squeeze(nanmean(FLIES.BLOCKS.greenChannel(:,:,dataCoords),[1,2]));
                newc = interp1([1:size(blirg,1)], [1:size(blirg,1)] , linspace(1,size(blirg,1), size(photCoords,2) ) )'; %Interpolated original volume number
                blirg = interp1([1:size(blirg,1)], blirg , linspace(1,size(blirg,1), size(photCoords,2) ) )'; %Bit inefficient to bring trans data up to phot size but eh
                figure
                plot( blirg )
                hold on
                plot( (photData(photCoords)*10)+130)
                title( ['Interpolated transient + phot'] )
            end
            %}
    
            %Interpolate sequence element across (usable) volume space
            %Old
            %{
            btSeqPosInterp = interp1([1:size(btData,1)], btData(:,5) , linspace(1,size(btData,1), tsImFrameCount ), 'previous' )'; %Interpolate btData 5th column to be imaging size
            btSeqInterp = randomSequence( btSeqPosInterp ); %Grab appropriate randomSequence elements for the timing
            btSeqInterpZ = reshape( btSeqInterp, (thisFlyBlock.Steps + thisFlyBlock.FlybackFrames), size( avg_z_green_aligned , 3 )  ); %Reshape as above
                %Basically represents for each volume the stimulus being presented (Changes mid-column represent the stim changing mid z-stack)
            %}
            %New
                %Both of these arrays generated by interpolation, not true time matching of frames to BT
            btSeqPosInterpZ = nan( size(imFrameListZ,1), size(imFrameListZ,2) ); %Make empty array for randomSequence to go into, of *total* imaging size
            btSeqPosInterpZ( imStimStart: imStimEnd) = interp1([1:size(btData,1)], btData(:,5) , linspace(1,size(btData,1), imStimEnd-imStimStart+1 ), 'previous' )'; %Interpolate and transplant randomSequence into 'frame' list
            btSeqPosInterpZ = btSeqPosInterpZ( : , adjImStimStartVol:adjImStimEndVol );
            %Repeat for raw btData row
            btRowPosInterpZ = nan( size(imFrameListZ,1), size(imFrameListZ,2) ); %Will hold info about which row of btData corresponds with frames, not i value
            btRowPosInterpZ( imStimStart: imStimEnd) = interp1([1:size(btData,1)], [1:size(btData,1)] , linspace(1,size(btData,1), imStimEnd-imStimStart+1 ), 'previous' )';
            btRowPosInterpZ = btRowPosInterpZ( : , adjImStimStartVol:adjImStimEndVol );

            %Check if data coming from forced singular Z plane
            if isfield( thisBlock, 'singularZ' ) && ~isempty( thisBlock.singularZ )
                disp(['Singular Z plane detected; Adjusting seq/vol interpolation for Z=',num2str(thisBlock.singularZ)])
                btSeqPosInterpZ = btSeqPosInterpZ( thisBlock.singularZ, : );
                btRowPosInterpZ = btRowPosInterpZ( thisBlock.singularZ, : );

            end

    
            %Pick a method to identify preceding stimuli by
            disp(['Interpolating bt sequence'])
            switch btInterpolationMethod
                case 'median'
                    %medBTSeqInterpZ = nanmedian( btSeqPosInterpZ , 1 ); %Time-matched (Imaging reference) list of stimulus present during that frame
                        %E.g. pos 360 being 2 means frame 360 was being presented with stimulus #2
                    methodBTSeqInterpZ = nanmedian( btSeqPosInterpZ , 1 );
                    methodBTRowInterpZ = nanmedian( btRowPosInterpZ , 1 );
                    %disp(['Interpolated with median'])
                case 'mode'
                    %modeBTSeqInterpZ = mode( btSeqPosInterpZ , 1 ); %Switch to mode, so that can never be non-integer
                    methodBTSeqInterpZ = mode( btSeqPosInterpZ , 1 );
                    methodBTRowInterpZ = mode( btRowPosInterpZ , 1 );
                    %disp(['Interpolated with mode'])
                case {'last' , 'intelligent'} %"aka, yaoi"
                    methodBTSeqInterpZ = btSeqPosInterpZ( end , : );
                    methodBTRowInterpZ = btRowPosInterpZ( end , : );
                    if isnan( methodBTSeqInterpZ(end) )
                        methodBTSeqInterpZ(end) = []; %May cause desync?
                        methodBTRowInterpZ(end) = [];
                        adjImStimEndVol = adjImStimEndVol - 1;
                        disp(['Adjusted imaging end volume reduced by 1 (',num2str(adjImStimEndVol+1),'->',num2str(adjImStimEndVol),') to account for NaN last element'])
                    end
                    if any( isnan(methodBTSeqInterpZ) ) %Allow no NaNs, because it crashes imStimTerp etc
                    %if nansum( isnan( methodBTSeqInterpZ ) ) > 1 %Allow 1 column (Presumably terminal) of NaN, crash if more
                        ['## Alert: Excessive NaN presence in last-interpolated bt sequence ##']
                        crash = yes
                    end
                    %disp(['Interpolated with last OR intelligent'])
            end
            disp(['Interpolated with ',btInterpolationMethod])
            %QA for critical aliasing issues
            %if numel(unique(medBTSeqInterpZ)) < numel( [medBTSeqInterpZ(1):medBTSeqInterpZ(end)] )
            %if numel(unique(modeBTSeqInterpZ)) < numel( [modeBTSeqInterpZ(1):modeBTSeqInterpZ(end)] )
            if numel(unique(methodBTSeqInterpZ)) < numel( [methodBTSeqInterpZ(1):methodBTSeqInterpZ(end)] )
                    %In theory detects if more elements likely existed than present in imaging frames
                        %The most common cause for this is imaging not being (much) faster than stimulation
                ['-# Warning: Sequence elements apparently lost due to aliasing/imaging framerate #-']
            end

            %kirara
    
             %Find actual stimuli delivered (via randomSequence)   
             if ~batteryDesign %Not battery (Rolling, Block, etc)
                 randomSeqOrig = matParamStruct.matSave.randomSequence; %Note that PTB script saves randomSequence as -1 (Blank), 1, (Stim 1), and 2 (Stim 2), rather than 5, 0, and 1, respectively
                 if ( includesBlanks == 0 && ( isequal( unique( randomSeqOrig ), [-1, 1, 2] ) ) || isequal( unique( randomSeqOrig ), [1, 2] ) ) ||...
                         ( includesBlanks == 1 && ( isequal( unique( randomSeqOrig ), [-2, -1, 1, 2] ) ) || isequal( unique( randomSeqOrig ), [-2, 1, 2] ) )
                     randomSeqCorr = randomSeqOrig;
                     randomSeqCorr( randomSeqCorr == -1 ) = 5;
                     randomSeqCorr( randomSeqCorr == 1 ) = 0;
                     randomSeqCorr( randomSeqCorr == 2 ) = 1; %Note order as critical, lest all 2 -> 1 -> 0

                     %imStimTerp = randomSeqCorr( medBTSeqInterpZ );
                     %imStimTerp = randomSeqCorr( modeBTSeqInterpZ );
                     imStimTerp = randomSeqCorr( methodBTSeqInterpZ );
                        %Note: Due to volume/timing inefficiencies this may skip the first actual presented element/etc
                            %Also, number of imaging frames/volumes per element/block may be inconsistent
                 else
                     ['## Alert: Unrecognised elements present (or absent) in randomSequence! ##']
                     unique( randomSeqOrig )
                     crash = yes
                 end
             else %Battery
                 %gorilla
                 iSeq = [1:matParamStruct.matSave.sequenceLength]; %Note: No guarantee that all elements will have been displayed, even for battery (But this will be QA'd later)
                 iSeqTerp = iSeq( methodBTSeqInterpZ ); %Theoretically represents which elements of i are associated with each (trimmed) volume of the recording
             
                 imStimTerp = matParamStruct.matSave.randomSequence( methodBTSeqInterpZ ); %Calculate this as well, for blank removal purposes

             end  
    
              %Trim data to relevant (i.e. Stimulated) portion
              %dataStimTrim = avg_z_green_aligned( :,:, adjImStimStartVol:adjImStimEndVol );
              dataStimTrim = thisImData( :,:, adjImStimStartVol:adjImStimEndVol );
                %Should we be concerned this doesn't exactly match imStimTerp for size?
              disp(['Data trimmed to stimulated portion only'])
              disp([size(dataStimTrim)])

              %Calculate times of volumes
              %if batteryDesign
              frameToVolList = imFrameListZ(1,adjImStimStartVol:adjImStimEndVol); %Basically which frame is associated with start of each vol
              volTimes = imFrameTimes( frameToVolList ); %What time (inferTimes reference) each volume occurred at [Accurate to inferTimes]
              volTimesDesc{1} = 'InferTimes of vol'; %Will store descriptions of each row
              volTimes(2,:) = methodBTRowInterpZ; %What row of btData each volume is associated with (Depending on method) [Interpolated, not true calculated]
               %If the preceding line crashes, there is desync between methodBTSeqInterpZ and volTimes
              volTimesDesc{2} = 'btData row of vol';
              volTimes(3,:) = frameOnsetIndices( frameToVolList );
              volTimesDesc{3} = 'TS index of vol';
              %end
              %motherbase

              %volTimes 'bootleg' adjustment
              %volTimesAdj = volTimes;
              %volTimesAdj(3,:) = volTimesAdj(3,:) - firstImStimFrameInd; %Manually shift TS reference back by 
              

              %Second testatory figure of transient against phot (if existing)
              %{{
              if hasPhotData && ~isShortcutting

                  %Pre calcs
                  %dataCoords = [imStimStart:imStimEnd]; %Using dataStimTrim this time, which is based on adjImStartVol
                  photCoords = [firstImStimFrameInd:lastImStimFrameInd]; %Note: As above, less guaranteed synchronicity with adjusted timings (Due to no equivalent adjustment to TS start/end timings)
                  blirg = squeeze(nanmean(dataStimTrim,[1,2]));
                  blirg = interp1([1:size(blirg,1)], blirg , linspace(1,size(blirg,1), size(photCoords,2) ) )'; %Bit inefficient to bring trans data up to phot size but eh
                  btlerg = btData( volTimes(2,:) , 5 );

                  %Old
                    %Misleading, as it plots Ca2+/phot data on 1:len indices, then randomSeq on volTimes indices, which are off by 1st element ~index
                    %{
                  figure
                  plot( blirg ) %Mean transient data
                  hold on
                  plot( (photData(photCoords)*10)+130) %Phot data
                  scatter( volTimes(3,:),  (btlerg*0.1)+134 ) %i value (1 in 10 dilution)
                  plot( volTimes(3,:), btData( volTimes(2,:) , 8 )+137 ) %Random sequence status
                  title( ['Interpolated transient + phot + 1/10th i + stim (Index)'] )
                    %}

                  %New
                  phoTimes = inferTimes(photCoords);
                  figure
                  plot( phoTimes, blirg' )
                  hold on
                  plot( phoTimes, (photData(photCoords)*10)+130)
                  plot( volTimes(1,:), btData( volTimes(2,:) , 8 )+135 ) %Use original TS reference (Time)
                  %title( ['Interpolated transient + phot (TS time ref)'] )
                  title( [strrep(flyID,'_',' '),' - interpolated transient + phot (TS time ref)'] )

              end
              %}

              %strobe

              %Testatory plot
              figure
              subplot(2,1,1)
              plot( squeeze( nanmean(dataStimTrim, [1,2] ) ) )
              hold on
              plot( imStimTerp+nanmean(dataStimTrim, 'all' )-7 )
              title('data and imStimTerp pre blank intro trim')

    
              %Siphon out blanks
                    %Start by extracting pre-experiment blank baseline (If existing)
              preBaselineExists = 0;
              if find( imStimTerp ~= 5, 1, 'first' ) ~= 1 && imStimTerp(1) == 5 

                  if batteryDesign
                    ['not coded for battery yet']
                    crash = yes
                  end

                  preBaselineExists = 1;

                  blankBaseline = [];
                  baseInds = [ 1:find( imStimTerp ~= 5, 1, 'first' )-1 ]; %Note: Assumption that there is no pre-baseline of stimuli/etc
                  blankBaseline = dataStimTrim( :,:,baseInds );
                  blankSeq = imStimTerp( baseInds ); %Not to be confused with blankSequence
                  blankVolTimes = volTimes(:,baseInds);
                  %QA
                  if numel( unique( blankSeq ) ) > 1
                      ['## Alert: Critical overfind in blank baseline collection ##']
                      crash = yes
                  end
                  dataStimTrim( :,:,baseInds ) = [];
                  imStimTerp( baseInds ) = [];
    
                  btSeqPosInterpZ( :, baseInds ) = []; %Trim also two relevant matrices
                  %medBTSeqInterpZ( baseInds ) = [];
                  %modeBTSeqInterpZ( baseInds ) = [];
                  methodBTSeqInterpZ( baseInds ) = [];

                  volTimes( :, baseInds ) = []; %Might break things
    
                  disp(['Blank pre-experiment baseline of ',num2str(numel(baseInds)),' elements collected; Data/sequence trimmed'])
              end

              %Finish plot
              subplot(2,1,2)
              plot( squeeze( nanmean(dataStimTrim, [1,2] ) ) )
              hold on
              plot( imStimTerp+nanmean(dataStimTrim, 'all' )-7 )
              title('data and imStimTerp post blank intro trim')

              %Follow up by collecting during-sequence blanks (if existing)
                    %(Only for inline blank handling though)
              if includesBlanks && blankHandleMode == 1
                  disp(['Blank trials will be removed inline (Mode 1)'])
                  %blankInds = [find(imStimTerp == -2)]; 
                  %blankStack = dataStimTrim( :,:, blankInds );
                  %kdrew

                  temp = bwlabel( imStimTerp == -2 ); %Label for all instances where stim was blank
                  invTemp = bwlabel( imStimTerp == 5 ); %Label for all other instances
                    %Note: Might be fragile based on imStimTerp processing/etc

                  disp([num2str(nanmax(temp)),' blank trial periods identified; Siphoning inline'])

                  preBlankStack = []; %Will hold all blank trial data
                  preBlankInds = []; %Ditto, for collection indices
                  preBlankIDs = []; %Will store the associated blank trial index (Used later for subspecification)
                  for row = 1:nanmax( temp )
                      thisBlankCollID = invTemp( find( temp == row, 1, 'last' )+1 ); %Find the ID of the stimulus collection period after this blank block
                        %Note: This ID will typically far exceed row, on account of blank trials being far lower rate than stimulated
                      %QA
                      if thisBlankCollID == 0
                          ['## Error: Fatal desync between labelling for blank collection ##']
                          crash = yes %Honestly not really possible, unless something weird happens with bwlabel
                      end
                      theseBlankCollInds = find( invTemp == thisBlankCollID );
                      preBlankStack = cat(3, preBlankStack, dataStimTrim( :,:, theseBlankCollInds ) ); %Append blank data along 3rd dim
                      preBlankInds = cat(3, preBlankInds, reshape( theseBlankCollInds, 1, 1, size(theseBlankCollInds,2) ) ); %Reshape inds and do same
                      preBlankIDs = cat(3, preBlankIDs, repmat( row, 1, 1, size(theseBlankCollInds,2) ) ); %Similar
                  end

                  blankInds = [find(imStimTerp == -2)]; %Still calculate this, for removal
                  dataStimTrim( :,:, blankInds ) = [];
                  imStimTerp( blankInds ) = [];
                  volTimes( :,blankInds ) = [];
              end
    
              %if thisBlock.blockNum == 3
              %prequel
              %end

              %Battery-specific processing
              if batteryDesign
                    %ichiban
                    clear imStimTerp %Prevent potential bleed
              end 
    
              %Perform bendy rolling specific processing
              if bendyBlockDesign == 0 && batteryDesign == 0
                  disp(['Processing for bendy rolling design'])
                  %Derive 'true' sequence delivered during rolling acquisition
                  randomSeqActual = randomSeqCorr( 1:methodBTSeqInterpZ(end) ); %Collect theoretical randomSequence from start to last element actually displayed
                    %Note: Heavily relies on assumption that imaging pre- and postcedes stimulus delivery
                  %Check if number of sent stimuli is a mod of nBack (unlikely)
                  if mod( length(randomSeqActual), nBack ) ~= 0
                      disp(['Bendy rolling sequence (Length ',num2str(length(randomSeqActual)),') requires trimming to match requested nBack (',num2str(nBack),')'])
                      %Trim
                      newLength =  length(randomSeqActual) - mod( length(randomSeqActual), nBack );
                      randomSeqActual = randomSeqActual( 1:newLength );
    
                      btSeqPosInterpZ( btSeqPosInterpZ > newLength ) = NaN; %Mostly for posterity. Note that *values* larger than newLength are removed, not indices
                      methodBTSeqInterpZ( methodBTSeqInterpZ > newLength ) = [];
    
                      dataStimTrim = dataStimTrim( :,:, 1:length(methodBTSeqInterpZ) ); %Use length of interpolated Z, since we are using stim to define imaging region now
                        %Should this actually be nVol frames after last stim? Or just ditch last few rolling stims?
                      imStimTerp = imStimTerp( 1:length(methodBTSeqInterpZ) ); 
                        %Reminder that imStimTerp technically relates to what was actually experienced by an imaging frame, not the true visible history                
                      volTimes = volTimes( :, 1:length(methodBTSeqInterpZ) );
                  end
    
                  %QA for stimulation amounts outnumbering imaging frames
                  if length(methodBTSeqInterpZ) < length( randomSeqActual )/nBack
                      ['-# Caution: Imaging frames (',num2str(length(methodBTSeqInterpZ) ),') outnumbered by stimulus elements (',num2str( length( randomSeqActual )/nBack),') #-']
                  end
    
                  %stimFrameInds = linspace( 0, size( dataStimTrim,3 ), length(randomSeqActual) );
                  stimFrameInds = linspace( 0, size( dataStimTrim,3 ), length(randomSeqActual) );
    
                  %Calculate nVol, for later use
                  %nVol = floor( size(imageStack,3)/ ( length(randomSequence) / options.nBack ) ); %Stolen from sortSEs2P initial implementation
                  nStimuli = options.nBack;
                  disp(['Using nBack of ', num2str(options.nBack),' to calculate acceptable volume counts'])
                  nVol = floor( size(dataStimTrim,3)/ ( length(randomSeqActual) / nStimuli ) ); %Calculates theoretically optimal number of volumes in imaging time period
                  while ceil( stimFrameInds( nStimuli ) - nVol ) < 0
                      disp(['-# Caution: Initially calculated nVol of ',num2str(nVol),' may be too large; Reducing #-'])   
                      nVol = nVol - 1;
                  end
                  disp(['Calculated nVol: ',num2str(nVol)])
                  %QA
                  if nVol <= 0
                      ['## Error in calculating nVol ##']
                      crash = yes
                  end
    
                  %Make data be exactly nVol * sequence long
                  %  Are there situations where original size might want to be preserved?
                  %stimFrameInds = linspace( 0, size( dataStimTrim,3 ), length(randomSeqActual) ); %Moved above to be more useful
    
                  %Check if mod math can be correctly applied
                  if nStimuli ~= 5
                      ['-# Alert: Mod calculations likely to be incorrect on account of non-standard nBack #-']
                      crash = yes
                  end
                  %samPoints = [ 0 : (length(randomSeqActual)/nStimuli-mod((nStimuli+4),8)) ]; %More or less just a list from 0 to sequenceLength/nStimuli
                  samPoints = [ 0 : (length(randomSeqActual)/nStimuli-mod((nStimuli+4),8)) ]; %More or less just a list from 0 to sequenceLength/nStimuli
    
                  lastVolFrameInds = stimFrameInds( samPoints*nStimuli + options.nBack); %Indices of last volume corresponding to each block (Note actually decimal, therefore needs to be ceiled/etc to function as inds)
                  newFrameInds = nan( nVol, size(lastVolFrameInds,2) );
                  for vol = 1:nVol
                      newFrameInds( vol, : ) = ceil( lastVolFrameInds ) - nVol + vol; %nVol rows, n of blocks cols
                  end
                  newFrameInds = reshape( newFrameInds, 1, nVol*size(newFrameInds,2));
                  %QA to ensure unique frames being grabbed
                  if numel( unique( newFrameInds ) ) ~= numel( newFrameInds )
                      ['## Alert: Same frame/s grabbed multiple times for bendy rolling reindicisation ##']
                      crash = yes
                      %May happen naturally with large volume numbers
                  end
    
                  dataStimTrim = dataStimTrim( :,:, newFrameInds ); %If this crashes, NaNs were probably in newFrameInds for some error reason
                  imStimTerp = imStimTerp( newFrameInds );
                  volTimes = volTimes( :, newFrameInds );
                  check validity of imStim/volT reindicising
                  disp( ['Bendy rolling data reindicised to length ',num2str( size(dataStimTrim,3) ),' (',num2str( size(dataStimTrim,3) / nVol ),' "blocks") for analysis (',...
                      num2str( ( 1 - ( size(dataStimTrim,3) / max(newFrameInds) ) ) * 100 ),'% Loss)'] )
                  %Loss in this context is how much of the original imaging data could not be associated with 'blocks' 
                    %It is likely to be highest when imaging rates are low and stimulus frequencies are high (i.e. 1 in 3 is 33%, but 1 in 8 is only 12.5%)
    
              end  

              %kc4k

              hasSiphoned = 0; %Flag to indicate whether blanks siphoned off
              stillRollable = 1; %Flag to indicate whether data can still be analysed rolling-style after blank removal
              phaseShifted = 0; %Default no
              if batteryDesign ~= 1 && ( any( imStimTerp == 5 ) || bendyBlockDesign == 1 ) && ~unsiphonedSEs
                  if bendyBlockDesign == 0 %"Blanks indicate true blank periods, intended for alternative analysis/etc"
                    %Note: This isn't really main rolling analysis here, just blank removal, unlike below for block design bendy
                        disp(['Rolling design; Siphoning blanks for alternative analysis'])
                        not coded for randomSeqActual yet
                        especially in context of reindicising
        
                        %blankStack = dataStimTrim( imStimTerp == 5 );
                        [inds] = find( imStimTerp == 5 );
                            %Note: For rolling bendy design, this can be thought of as an actual 'blank' period
                            %          For bendy block design, blanks are actually the analagous recording period
                        %blankStack = dataStimTrim( :,:, inds );
                        baselineBlankStack = dataStimTrim( :,:, inds );
                        %blankInds = inds;
                        baselineBlankInds = inds;
                        baselineBlankTimes = volTimes( :, inds );
            
                        %Find if blanks occurred anywhere other than start
                            %Deprecated on account of baseline extraction now
                        %{
                            %Note that a rolling paradigm with blanks only at the start doesn't exist yet
                        temp = zeros( 1, size(imStimTerp,2) );
                        temp( imStimTerp == 5 ) = 1;
                        if numel( unique( bwlabel(temp) ) ) > 2 %"If blanks are non-contiguous"
                            stillRollable = 0; %Blanks throughout
                            disp(['Non-contiguous blanks present; Cannot do subsequent rolling analysis'])
                        else
                            stillRollable = 1; %Blanks only at start
                            disp(['Blanks contiguous; Can do rolling analysis'])
                        end
                        %}
                        stillRollable = 0; %If any blanks remaining, by inference they must be during exp
            
                        dataStimTrim( :,:,inds ) = [];
                        imStimTerp( inds ) = [];
                        volTimes( :, inds ) = [];
                        disp(['Blank data siphoned and removed from data'])
                        hasSiphoned = 1;
    
    
                        %Note: imStimTerp duplication due to aliasing not addressed for this condition yet
    
                  else %"Blanks are intended recording period" (a la LED block design)
    
                        %Section readme:
                        %{
                        btData (PTB reference frame)
                        is interpolated into
                        btSeqPosInterpZ (imaging reference, values represent what element of randomSeq being presented by PTB)
                        which is medianed into
                        medBTSeqInterpZ (imaging reference, values indicate randomSequence element)
                        which is used to calculate 
                        imStimTerp (imaging reference, values indicate raw stimulus [with potential aliasing duplication])
                        
                        which is copied as
                        temp (imaging reference, raw stimulus, derived from imStimTerp)
                        is labelled for inter-stimulus periods and used to find
                        intInds (imaging reference, values indicate start of inter-stimulus periods [imaging reference])
                        which is extrapolated out as
                        collInds (imaging reference, forced to be nominal number of frames [e.g. 6] after stimulus cessation)
                        
                        which is also used to make
                        stimInds (imaging reference, values indicate imaging frame #)
                        %}
    
                        disp(['Analysing for bendy block design'])
                        stillRollable = 0; %Get this out of the way initially
    
                        %Find inter-stimulus (i.e. Blank) periods
                        [inds] = find( imStimTerp == 5 );
    
                        %Acquire block length data
                        blockLength = matParamStruct.matSave.blockLength;
                        %Compare against nBack
                        if blockLength ~= nBack
                            ['-# Caution: Derived blockLength (',num2str(blockLength),') differs from requested or default nBack (',num2str(nBack),')']
                        end
    
                        %Quick estimation of number of blocks
                        if exist('matParamStruct') && exist('btData') && size( btData,2 ) >= 6
                            temp = matParamStruct.matSave.randomSequence(1: btData(end,5) );
                            disp(['Estimated number of blocks delivered: ', num2str(nansum( temp ~= -1 ) / blockLength) ])
                            disp(['(~',num2str(1 / (nansum( temp ~= -1 ) / blockLength / btData(end,6))),' s/block)']) %Does not account for pre-baseline throwing off timing slightly
                                %If inter-stim periods are set to equal duration, applies to those as well obvs
                        end
    
                        %Find block-style random sequence
                        imStimTerpPure = imStimTerp( imStimTerp ~= 5 );
                            %Note: No safeguard for interrupted stim/rest block/etc
                        numFullBlocks = floor( size(imStimTerpPure ,2) / blockLength ); %How many full stimulus blocks were presented in the imaging time
    
                        %Back-calculate inter-stimulus period
                        temp = imStimTerp;
                        temp(temp ~= 5) = 0; %Make labellable version of inter-stimulus periods
                        temp = bwlabel( temp ); %Label
                        disp(['There were ',num2str(nanmax(temp)),' inter-stimulus periods']) %Note: May be truncated at end
                        temp2 = repmat( [1:nanmax(temp)]' , 1 , size(temp,2) ); %Assemble comparative list of labels
                        temp3 = temp == temp2; %Query presence of labels in labelled sequence
                            %Note: No (current) guarantee of full imaging frame (Could be truncated)
                        interStimPeriods = nansum( temp3,2 );
                        disp([ 'Period lengths: ', num2str(unique(interStimPeriods)') ]) %More than 1 element means variable inter-stimulus periods
                        disp(['(',num2str(( nansum(interStimPeriods == unique(interStimPeriods)') / numel(interStimPeriods) )*100),' % of all periods, respectively)'])
    
                        %Find inter-stimulus periods and collect minimum 'normal' number of frames
                        interUn = unique(interStimPeriods)';
                        interUn( [nansum(interStimPeriods == interUn) < 0.05*numel(interStimPeriods)] ) = []; %Remove inter-stim periods accounting for only a small fraction of events (Likely to be last element/etc)
                        nomInter = min( interUn ); %Nominal minimal inter-stim period
                            %Note: This can still only be like, 25% or less of all inter-stimulus periods
                        disp(['Selected inter-stim period to collect: ', num2str(nomInter),' (',num2str(( nansum(interStimPeriods == nomInter) / numel(interStimPeriods) )*100),'% of original periods)'])
    
                        %timeoflove

                        %Find indices of inter-stimulus periods
                        [~,intInds] = max(temp3,[],2); %First element of each respective inter-stim period
                        [~,outInds] = max(fliplr(temp3),[],2); %End of each inter-stim period, flipped reference frame
                        outInds = size(temp3,2) - outInds + 1; %Fix reference frame
    
                        %Assemble inds for collection
                        collInds = intInds + repmat( [1:nomInter], size(intInds,1), 1 ) - 1; %Note: Original phase of inter-stimulus obviously obliterated here
                        disp([ 'Last frame of imaging to be collected: ', num2str(collInds(end)) ])
                        disp([ 'Total imaging size: ', num2str( size(dataStimTrim,3) ) ])

                        %Make sure last element of collection not outside data
                        if collInds(end) > size(dataStimTrim,3)
                            collInds(end,:) = []; %Remove one imaging period
                            disp(['Terminal collection period removed; New end: ', num2str(collInds(end)) ])
                        end                       

                        %Another testatory figure
                        meanStimTrimTrans = nanmean(dataStimTrim, 'all' );
                        figure
                        plot( squeeze( nanmean(dataStimTrim, [1,2] ) ) )
                        hold on
                        plot( imStimTerp+meanStimTrimTrans-8 )
                        for row = 1:size(collInds,1)
                            line( [ collInds(row,1),collInds(row,end) ], [meanStimTrimTrans+5,meanStimTrimTrans+5], 'Color', 'g' )
                            if row > 1
                                line( [ collInds(row-1,end),collInds(row,1) ], [meanStimTrimTrans+6,meanStimTrimTrans+6], 'Color', 'r' )
                            end
                        end
                        xlim([1,collInds(5,end)]) %Free to resize
                        title(['Meaned data + collection period (g) + stim period (r) + sequence (bottom) [First 5 trials window]'])
    
                        %Prepare a labelled form of randomSeqCorr for intelligent stimulus collection, if applicable
                        if isequal( btInterpolationMethod, 'intelligent' ) 
                            if ~includesBlanks %No blanks
                                %temp = randomSeqCorr;
                                randomSeqForLabel = randomSeqCorr;
                            else %Blanks
                                %temp = randomSeqCorr( randomSeqCorr ~= -2 );
                                %kairi
                                if blankHandleMode == 1 %Treat differently
                                    disp(['Blanks will be removed at end (Mode 2)'])
                                    randomSeqForLabel = randomSeqCorr( randomSeqCorr ~= -2 );                                    
                                else %Treat 'same'
                                    randomSeqForLabel = randomSeqCorr;
                                end
                            end
                            temp = randomSeqForLabel;

                            temp( temp == 0 ) = 1;
                            if includesBlanks && blankHandleMode == 2
                                temp( temp == -2 ) = 1;                                
                            end
                            temp( temp == 5 ) = 0;
                            %if includesBlanks
                            %    temp( temp == -2 ) = 0; %Incorrect; To match with imStimTerp, these elements must be *removed*
                            %end
                            randomSeqLabel = bwlabel( temp );
                        end
    
                        %melting
    
                        %Assemble stimulus periods for collection
                            %Note: This relies on the assumption that there were no huge issues in theoretical->actual stimulus presentation
                        stimInds = collInds(:,1) - repmat( fliplr([1:blockLength]), size(collInds,1), 1 ); %Identify stimulated frame periods (Not for collection, but to derive sequence)
                            %Use collInds, not intInds, because cleaned
                            %Note: Due to frame attribution inaccuracies/phase differences between imaging and display, there will be variability in number of 'stimulus' frames preceding an inter-stimulus period
                        %Quick QA
                        if nanmax( nanmax( stimInds ) ) > size(methodBTSeqInterpZ,2)%methodBTSeqInterpZ(end)
                                %stimInds max represents last imaging volume *during* stimulation, size of methodBTSeqInterpZ indicates how many valid frames collected (Values indicate associated i element)
                            ['## Alert: Requested post-stimulus imaging frame (',num2str(nanmax( nanmax( stimInds ) )),') beyond apparent imaging duration (',num2str(methodBTSeqInterpZ(end)),') ##']
                            crash = yes
                                %This is marginally more likely with last due to lastNaN removal
                        end
                        %stimSeq = randomSeqCorr( medBTSeqInterpZ( stimInds(:,:) ) ); %Pulls duplicates
                        %stimSeq = randomSeqCorr( [-blockLength+1:0] +  medBTSeqInterpZ( stimInds(:,end) ) );
                        stimSeq = nan( size(stimInds) );
                        stimSeqCorrInds = nan( size(stimInds) );
                        randomSeqBlock = NaN; %Predefine
                        for row = 1:size( stimInds,1 )
                           % burnthem
                           switch btInterpolationMethod
                               case {'median' , 'mode', 'last'}
                                    %stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  medBTSeqInterpZ( stimInds(row,end) ); %Which elements of randomSeqCorr were pulled
                                    %stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  modeBTSeqInterpZ( stimInds(row,end) );
                                    stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  methodBTSeqInterpZ( stimInds(row,end) );
                                        %pull medBTSeqInterpZ values at stimInds, go 5 back, pull those elements of randomSeqCorr
                                    %stimSeq( row, : ) = randomSeqCorr( stimSeqCorrInds( row, : ) ); %Stimuli
                                    stimSeq( row, : ) = randomSeqForLabel( stimSeqCorrInds( row, : ) ); %Stimuli
                               case 'intelligent'
                                   thisStimEnd =  methodBTSeqInterpZ( stimInds(row,end) );
                                   %stimSeqCorrInds( row, : ) = stimInds(row,:); %Note: In this incarnation, this is more of a feelgood approximation of stim volumes than perfect reporter
                                        %Due to difference in how Intelligent mode collates stim identity
                                        %Note: This is subsampled to nBack and thus less than useful?
                                   %Check for accidentally landing in inter-stimulus period
                                        %Note: A rolling phase adjustment might be useful, but would have to be justified from first principles (Why would phase become unaligned over time?)
                                   if randomSeqLabel( thisStimEnd ) == 0
                                       if randomSeqLabel( thisStimEnd-1 ) ~= 0
                                           thisStimEnd = thisStimEnd - 1; %Adjust        
                                       elseif randomSeqLabel( thisStimEnd-2 ) ~= 0
                                           thisStimEnd = thisStimEnd - 2; %Same, but more
                                           disp(['Caution: Trial/Row #',num2str(row),' adjusted by -2 to keep phase'])
                                       elseif randomSeqLabel( thisStimEnd-3 ) ~= 0
                                           thisStimEnd = thisStimEnd - 3; %SAME, BUT EVEN MORE
                                           disp(['CAUTION: TRIAL/ROW #',num2str(row),' ADJUSTED BY -3 TO KEEP PHASE'])
                                       else
                                            ['## Fatal failure to intelligently identify stimulus block identity ##']
                                            crash = yes
                                            %Being 1 stimulus before the end of the block = okay, being 1 after = acceptable, but any more probably indicates a large error has occurred
                                            %Likely some grand loss of phase
                                       end
                                   end
                                   %QA, to make sure same block isn't used twice (Phase check, basically)
                                   if randomSeqLabel( thisStimEnd ) == randomSeqBlock || randomSeqLabel( thisStimEnd ) == 0
                                       ['## Alert: Stimulus period captured twice OR failure to identify trial # at end of stimulus period ##']
                                       ['Row: ',num2str(row),', prev. period #: ',num2str(randomSeqBlock),', this period #: ',num2str(randomSeqLabel( thisStimEnd ) )]
                                   end
                                   randomSeqBlock = randomSeqLabel( thisStimEnd );
                                   %Potentially incorrect QA
                                   if randomSeqBlock ~= row
                                       ['-# Potential phase loss? #-'] %This may be an incorrect assumption for different bendy designs/etc
                                       crash = yes
                                   end
                                   thisStimInds = [ find( randomSeqLabel == randomSeqBlock,1, 'first' ) : find( randomSeqLabel == randomSeqBlock,1, 'last' ) ];
                                   %QA
                                   if numel(thisStimInds) ~= blockLength
                                       ['## Fatal mismatch between identified stimulus block length and theoretical block length ##']
                                       crash = yes
                                   end
                                   %If no issues, use labelled region to get stimuli
                                   %stimSeq( row, : ) = randomSeqCorr( thisStimInds );
                                   stimSeq( row, : ) = randomSeqForLabel( thisStimInds );
                           end
                        end
                        %QA
                        if any( stimSeq == 5, 'all' )
                            ['## Alert: Blanks collected in stimSequence; Probable phase failure ##'] %phailure
                            crash = yes
                                %If using interpolation other than last, this may be result of median/mode 'skipping' stimulus elements 
                                    %i.e. >1 stimuli occurred in the span of a single volume
                        end

                        %Average stim  plot
                        if savePhasePlots
                            %Check/Make folder to put figures into
                            phaseFolder = strcat( dataFolder,filesep,'PHASE' );
                            if exist(phaseFolder) ~= 7
                                mkdir( phaseFolder )
                                disp(['Phase figure folder made at ',phaseFolder])
                            end
                        end
                        %(Old position)
                        %{
                        temp = squeeze( nanmean( dataStimTrim, [1,2] ) );
                        %temp2 = temp( stimSeqCorrInds );
                        temp2 = [];
                        temp2 = nan( size(collInds,1), floor(nanmax(diff(collInds(:,end)))*1.25) );
                        coords = [1:collInds(1,1)-1]; %Collect pre first collInds stim period
                        temp2( 1, 1:length(coords) ) = temp(coords);
                        for row = 2:size(collInds,1) 
                            coords = [collInds(row-1,end)+1:collInds(row,1)-1]; %Collect pre first collInds stim period
                            temp2( row, 1:length(coords) ) = temp(coords);
                        end
                        figure
                        %plot( nanmean(temp2,1) );
                        %Stim
                        errorbar( nanmean(temp2,1), nanstd(temp2,[],1) )
                        %}
                        %Average stim and collection plot
                        temp = squeeze( nanmean( dataStimTrim, [1,2] ) );
                        %temp2 = temp( stimSeqCorrInds );
                        temp2 = [];
                        temp2 = nan( size(collInds,1), floor(nanmax(diff(collInds(:,end)))*1.25) );
                        coords = [1:collInds(1,1)-1]; %Collect pre first collInds stim period
                        temp2( 1, 1:length(coords) ) = temp(coords);
                        for row = 2:size(collInds,1) 
                            coords = [collInds(row-1,end)+1:collInds(row,1)-1]; %Collect pre first collInds stim period
                            temp2( row, 1:length(coords) ) = temp(coords);
                        end
                        figure
                        %plot( nanmean(temp2,1) );
                        %Stim
                        subplot(2,1,1)
                        errorbar( nanmean(temp2,1), nanstd(temp2,[],1) / sqrt( size(temp2,1) ) )
                        legend({'Stim'})
                        %hold on
                        %for row = 1:size(temp2,1)
                        %    plot( temp2(row,:) )                      
                        %end
                        title([strrep(flyID,'_',' '),' - [Whole] Stim period average transient'])
                        xlabel(['Time (vol)'])
                        %Note that this plot uses inter-collInds for 'stim' collection, not technically any stim inds
                        %Coll
                        temp3 = temp(collInds);
                        if exist('arbPhaseShift')
                            coords = collInds + arbPhaseShift;
                            [temp5,~] = nanmax( coords, [], 2);
                            [temp6,~] = nanmin( coords, [], 2);
                            if any(temp5 > size(temp,1))
                                coords( find(temp5 > size(temp,1)), : ) = [];
                            end
                            if any(temp6 < 1)
                                coords( find(temp6 < 1), : ) = [];
                            end
                            temp4 = temp(coords);
                        end
                        subplot(2,1,2)
                        errorbar( nanmean(temp3,1), nanstd(temp3,[],1) / sqrt( size(temp3,1) ) )
                        if exist('arbPhaseShift')
                            hold on
                            errorbar( nanmean(temp4,1), nanstd(temp4,[],1) / sqrt( size(temp4,1) ) )
                        end
                        %title([strrep(flyID,'_',' '),' - Coll period average transient'])
                        titleStr = ['Fly #',num2str(thisBlock.flyNum),'-',num2str(thisBlock.blockNum),' - ',strrep(flyID,'_',' '),' - Coll period average transient'];
                        if exist('arbPhaseShift')
                            titleStr = [titleStr,char(10),'(And phase-shifted [',num2str(arbPhaseShift),'] version)'];
                            legend([{'Orig.'},{'Phase shifted'}])
                        else
                            legend({'Orig.'})
                        end
                        xlabel(['Time (vol)'])
                        title(titleStr)
                        if savePhasePlots
                            drawnow
                            saveas(gcf, [phaseFolder,filesep,flyID,'_phasePlot'], 'png');
                            disp(['Phase plot saved'])
                        end

                        %Interruption to siphon blanks if applicable and moded
                        if includesBlanks && blankHandleMode == 2
                            blankTrialIDs = find( nansum( stimSeq == -2, 2 ) == blockLength ); %Identify which trials were blank
                                %Note: No explicit check for sequences that were mixed or aberrant blanks/ISI/etc
                            disp([num2str(length(blankTrialIDs)),' blank trials found in final sequence; Siphoning outline'])

                            blankCollInds = collInds( blankTrialIDs, : );
                            deImBlankInds = reshape( collInds( blankTrialIDs, : )', 1, size(collInds( blankTrialIDs, : ),1)*nomInter ); %Based on deImInds as done below, obviously
                            %postBlankData = dataStimTrim( :,:, deImBlankInds );
                            blankStack = dataStimTrim( :,:, deImBlankInds );
                            blankInds = deImBlankInds;
                            blankSequence = reshape(stimSeq( blankTrialIDs, : )', 1, size(stimSeq( blankTrialIDs, : ),1)*blockLength);
                            blankTimes = volTimes( :, deImBlankInds ); %Not checked for same size as dataStimTrim

                            %QA
                            if mod( size( blankStack, 3), nomInter ) ~= 0
                                ['## Alert: Potential phase loss in post-stim blank data ##']
                                crash = yes
                            end

                            %Remove blanks from 'real' data
                            stimSeq( blankTrialIDs , : ) = [];
                            collInds( blankTrialIDs , : ) = [];

                            %Testatory plot of blanks
                            figure
                            errorbar( nanmean(temp(blankCollInds),1), nanstd(temp(blankCollInds),[],1) / sqrt(size(blankCollInds,1)) )
                            title([strrep(flyID,'_',' '),' - Blank coll period average transient'])
                            xlabel(['Time (vol)'])
                            if exist('arbPhaseShift')
                                try %Too lazy to write properly
                                    hold on
                                    errorbar( nanmean(temp(blankCollInds+arbPhaseShift),1), nanstd(temp(blankCollInds),[],1) / sqrt(size(blankCollInds,1)) )
                                    legend([{'Orig. blank'},{'Phase shifted blank'}])
                                catch
                                    disp([('(Failure to append phase shifted blank data to testatory plot)')])
                                end
                            end

                        end
    
                        %Reshape derived stimulus sequence
                        deRandomSeq = reshape(stimSeq', 1, size(stimSeq,1)*blockLength); %Transposition v important here for correct rowwise-ness
                         

                        %Apply phase shift equally to all collection indices
                            %Do this immediately prior to actual collection
                        %phaseShifted = 0; %Default no (Moved above in case of unsiphoned)
                        if exist('arbPhaseShift') && ~isempty(arbPhaseShift)
                            collInds = collInds + arbPhaseShift;
                                %Note: Does not shift stimulus inds, by design
                            disp(['Phase shift applied'])
                            [temp5,~] = nanmax( collInds, [], 2);
                            [temp6,~] = nanmin( collInds, [], 2);
                            if any(temp5 > size(dataStimTrim,3))
                                collInds( find(temp5 > size(dataStimTrim,3)), : ) = []; %Delete first collection period
                                deRandomSeq( 1:blockLength ) = []; %And corresponding stim period
                                disp(['(Trial/s ',num2str(find(temp5 > size(dataStimTrim,3))),' had to be ditched due to overrun)'])
                            end
                            if any(temp6 < 1)
                                collInds( find(temp6 < 1), : ) = [];
                                deRandomSeq( end-blockLength+1:end ) = [];
                                disp(['(Trial/s ',num2str(find(temp6 < 1)),' had to be ditched due to underrun)'])
                            end
                            phaseShifted = 1;
                        end

                        %Now to collect frames
                        deImInds = reshape( collInds', 1, size(collInds,1)*nomInter ); %A run-on list of what volumes relate to stimulation
                            %Note that collInds holds the original, segmented version of volume positions
                        postStimData = dataStimTrim( :,:, deImInds ); %Only valid as long as dataStimTrim timing synchronous with other matrices
                            %Note: Unless imaging was 1:1 framerate with display, this will not necessarily match deRandomSeq in size
                        %QA for correct phase-size
                        if mod( size( postStimData, 3), nomInter ) ~= 0 || mod( size( deRandomSeq, 2), blockLength ) ~= 0
                            ['## Alert: Potential phase loss in post-stim data and/or derived random sequence ##']
                            ['Imaging size: ',num2str(size( postStimData ))]
                            ['Random seq size: ',num2str( size(deRandomSeq,2) )]
                            crash = yes
                        end
                        postStimTimes = volTimes( :, deImInds );
    
                        %Report
                        disp(['Final number of imaging events: ',num2str(size( postStimData,3 )/nomInter)]) %Add potential max # imaging events here
                        disp(['Final number of stimulus events: ',num2str(size( deRandomSeq,2 )/blockLength)]) %Add potential total (btData) here
                        if includesBlanks && blankHandleMode == 2
                            disp(['(And ',num2str(length(blankTrialIDs)),' blank trials)'])
                        end

                        %Interim plot if requested
                        if doBendyTransientGraph
                            temp = [];
                            temp{1} = squeeze( nanmean( dataStimTrim, [1,2] ) ); %Mean of imaging 
                            temp{2} = temp{1}( collInds ); %2D array of means (Event, Time)
                            temp{6} = nanmean(temp{2},1); %1D mean of transients
                            temp{8} = nanstd( temp{2}, [], 1 ) / sqrt( size(temp{2},1) ); %1D SEM of transients
                            if hasPhotData && cleanPhotData
                                temp{3} = volTimes(3,:); %All TS timepoints associated with volumes
                                temp{4} = temp{3}( collInds ); %2D array of TS timepoints relating to specific volumes
                                temp{5} = photProc( temp{4} ); %2D array of phot data during specific volumes
                                temp{7} = nanmean( temp{5}, 1 ); %normalize( nanmean( temp{5}, 1 ), 'Range', [nanmin( temp{6} ) - range( temp{6} ) , nanmin( temp{6} ) ] ); %1D mean of phot
                                temp{9} = nanstd( temp{5}, [], 1 ) / sqrt( size(temp{5},1) ); %1D SEM of phot
                            end
                            temp{10} = [];
                            for vol = 1:size( collInds,2 )
                                temp{10}(:,:,vol) = nanmean( dataStimTrim(:,:, collInds(:,vol) ), 3 );
                            end
                            [~, temp{11}(1)] = nanmin( temp{6} ); temp{11}(2) = floor( size(temp{6},2) /2 ); [~, temp{11}(3)] = nanmax( temp{6} );

                            %falling
                            figure
                            subplot( 2,3, [1:3] )
                            errorbar( temp{6}, temp{8} )
                            if hasPhotData && cleanPhotData
                                hold on
                                errorbar( temp{7}*range(temp{6}) + nanmin( temp{6} ) - 0.5*range(temp{6}), temp{9}*range(temp{6}) )
                                legend([{'Transient'},{'Phot'}])
                            else
                                legend([{'Transient'}])
                            end
                            if ~exist('arbPhaseShift')
                                title(['Block design transient plot - ',strrep(flyID,'_',' '),' (','#',num2str(thisBlock.flyNum),'-',num2str(thisBlock.blockNum),')'])
                            else
                                title(['Block design transient plot - ',strrep(flyID,'_',' '),' (Phase shifted)',' (','#',num2str(thisBlock.flyNum),'-',num2str(thisBlock.blockNum),')'])
                            end
                            xlabel(['Volume'])
                            ylabel(['Pixel intensity (a.u.)'])
                            xlim([1,size(temp{6},2)])

                            subplot( 2,3, [4] )
                            imagesc( temp{10}(:,:, temp{11}(1) ) )
                            title(['Min. frame (Vol #',num2str(temp{11}(1)),')'])
                            subplot( 2,3, [5] )
                            imagesc( temp{10}(:,:, temp{11}(2) ) )
                            title(['Mid. frame (Vol #',num2str(temp{11}(2)),')'])
                            subplot( 2,3, [6] )
                            imagesc( temp{10}(:,:, temp{11}(3) ) )
                            title(['Max. frame (Vol #',num2str(temp{11}(3)),')'])
                            if savePhasePlots
                                drawnow
                                saveas(gcf, [phaseFolder,filesep,flyID,'_transientPlot'], 'png');
                                disp(['Phase transient plot saved'])
                            end

                            clear temp
                        end

                        %Process blank data into something valid (if applicable)

                        %to do: blank 'transient' plot (w/ phot)

                        %crawl
                        if includesBlanks && blankHandleMode == 1 && ~isempty(preBlankStack)
                            %elysium
                            disp(['Subselecting ',num2str(nanmax(preBlankIDs)),' blank trial periods to match data length of ',num2str(nomInter)])
                            blankStack = [];
                            blankInds = []; %Reuse (or rather, regenerate from ghost)
                            for row = [1:nanmax(preBlankIDs)]
                                thisData = nan( [size(preBlankStack,[1,2]), nomInter] ); %Preallocate as NaNs
                                theseCoords = find(preBlankIDs == row);

                                safeness = min(length(theseCoords) , nomInter);
                                    %Note: If any one blank collection period shorter than safeness, NaNs will be inserted into blankStack
                                thisData(:,:, [1:safeness] ) = preBlankStack( :,:, [theseCoords(1):theseCoords(safeness)] );

                                blankStack = cat(3, blankStack,thisData ); %Append blank data along 3rd dim
                                blankInds = cat(3, blankInds, reshape([theseCoords(1):theseCoords(safeness)],1,1, safeness) );
                            end

                            %QA
                            if size( blankStack,3 ) ~= nanmax(preBlankIDs)*nomInter
                                ['## Alert: Correct sizing failure in assembly of blank stack ##']
                                crash = yes
                            end
                            
                        end
    
                  end %bendyBlockDesign end

              elseif ~batteryDesign && unsiphonedSEs 
                  disp(['-# Not siphoning SEs data by request #-'])
                    %NOTE: NOT TESTED FOR ROLLING
              elseif batteryDesign
                  disp(['-- Battery design; Siphoning unnecessary --'])
              else
                  ['## unknown case ##']
                  crash = yes
              end

              %Collect appropriate processed photodiode data, now that volTimes is finished being modified presumably
              if hasPhotData && cleanPhotData
                  if ~batteryDesign && ~unsiphonedSEs
                      if bendyBlockDesign == 0
                          photProcVol = photProc( volTimes(3,:) )'; %Collect only phot data from imaging volumes
                            %Note: Will crash if volTimes does not exist, which may occur
                            %Secondary note: A better system might collect *all*
                            %the phot data between imaging vols, rather than just
                            %the point value
                      else
                          photProcVol = photProc( postStimTimes(3,:) )'; %Note different vol/time selection for bendy block                        
                      end
                  elseif batteryDesign || unsiphonedSEs
                      photProcVol = photProc( volTimes(3,:) )'; %Ostensibly same as bendy rolling                      
                  end
              end

              %--------------------------------------------------------------------------------------------
    
              %callgravity

              %Overwrite data (if function source)
              if dataSource == -1
                  %thisBlock.greenChannel = dataStimTrim; %Mostly unnecessary
                  %thisBlock.randomSequence = imStimTerp; %Mostly unnecessary
                  %BLOCKS(thisFlyRowInd,:) = thisBlock;
                  if ~batteryDesign && ~unsiphonedSEs %Rolling, Block, etc
                      %abort
                      if bendyBlockDesign == 0 %Rolling
                          BLOCKS( thisFlyRowInd ).greenChannel = dataStimTrim;
                          %BLOCKS( thisFlyRowInd ).randomSequence = imStimTerp;
                          BLOCKS( thisFlyRowInd ).randomSequence = randomSeqActual; %New, not interpolated
                            %Note: Using randomSequence not of exact same length as imaging may cause issues with rolling implementation in analyseBlock
                          BLOCKS( thisFlyRowInd ).volTimes = volTimes; %Times
                          BLOCKS( thisFlyRowInd ).nVol = nVol;
                          BLOCKS( thisFlyRowInd ).nStimuli = nStimuli;
                          BLOCKS( thisFlyRowInd ).stimulus = 'bendy_rolling';
                          %clear dataStimTrim imStimTerp randomSeqActual nVol nStimuli %Just in case
                      else %Block
                          BLOCKS( thisFlyRowInd ).greenChannel = postStimData; %Data
                          BLOCKS( thisFlyRowInd ).randomSequence = deRandomSeq; %Sequence
                          BLOCKS( thisFlyRowInd ).volTimes = postStimTimes; %Times
                          BLOCKS( thisFlyRowInd ).imagingInds = deImInds; %Vol #s collected
                          BLOCKS( thisFlyRowInd ).nVol = nomInter;
                          BLOCKS( thisFlyRowInd ).stimulus = 'bendy_block';
                          %BLOCKS( thisFlyRowInd ).blankBlocks = 0; %Need to add support later for blank blocks
                          BLOCKS( thisFlyRowInd ).blankBlocks = includesBlanks; %Need to add support later for blank blocks
                          BLOCKS( thisFlyRowInd ).fauxBlockDesign = 1; %Just to keep track
                          disp(['Faux-block design created'])
                          %clear postStimData deRandomSeq nomInter
                      end
                      if isfield( matParamStruct.matSave, 'blockDesign' )
                         BLOCKS( thisFlyRowInd ).bendyBlockDesign =  matParamStruct.matSave.blockDesign;
                      end
                      disp(['Modified data and randomSequence inserted into BLOCKS'])
                      %slyleaf
                      if stillRollable == 0
                          BLOCKS( thisFlyRowInd ).isRolling = 0;                
                          disp(['Rolling status revoked'])
                      end
                      %Blanks (if applicable)
                      if includesBlanks && ~isempty(blankStack)
                          BLOCKS( thisFlyRowInd ).blankImageStack = blankStack; %Note: Architecture ostensibly should be repeating, 'nVol' sized groups of frames                  
                          if blankHandleMode == 2 %blankSequence calcs too annoying to do for inline siphoning
                            BLOCKS( thisFlyRowInd ).blankSequence = blankSequence;
                          end
                          BLOCKS( thisFlyRowInd ).blankImageInds = blankInds; %Note: blankInds nature may differ depending on blank handling mode
                          BLOCKS( thisFlyRowInd ).blankImageTimes = blankTimes; %Might crash depending on handling mode?
                          disp(['Modified blank data (and inds) inserted into BLOCKS'])                          
                      end
                  elseif batteryDesign %Battery
                      BLOCKS( thisFlyRowInd ).greenChannel = dataStimTrim;
                      BLOCKS( thisFlyRowInd ).iSequence = iSeqTerp;
                      BLOCKS( thisFlyRowInd ).volTimes = volTimes;                    
                      BLOCKS( thisFlyRowInd ).stimulus = 'battery';
                      %if hasPhotData
                      %    BLOCKS( thisFlyRowInd ).photData = photData; %Will blow out size a bit probably
                      %end
                  elseif ~batteryDesign && unsiphonedSEs
                      BLOCKS( thisFlyRowInd ).greenChannel = dataStimTrim;
                      BLOCKS( thisFlyRowInd ).volTimes = volTimes;                    
                      BLOCKS( thisFlyRowInd ).stimulus = 'bendy_block_unsiphoned';
                  end
                  BLOCKS( thisFlyRowInd ).hasPhotData = hasPhotData;
                  if hasPhotData
                      BLOCKS( thisFlyRowInd ).photData = photData; %Will blow out size a bit probably
                      if exist('photProcVol')
                          BLOCKS( thisFlyRowInd ).photProcVol = photProcVol;
                      end
                      disp(['~~ Photodiode data saved to BLOCKS ~~'])
                  end

                  if exist( 'preBaselineExists' ) && preBaselineExists %Pre-experiment blank baseline
                      BLOCKS( thisFlyRowInd ).preBaseline = blankBaseline;     
                      BLOCKS( thisFlyRowInd ).preBaseSeq = blankSeq;
                      BLOCKS( thisFlyRowInd ).preBaseTimes = blankVolTimes;
                  end
                  %BLOCKS( thisFlyRowInd ).volTimes = volTimes; %Note: Does not (presently) account for blank removal/etc
                  BLOCKS( thisFlyRowInd ).syncModified = 1;
                  BLOCKS( thisFlyRowInd ).phaseShifted = phaseShifted;
              end
    
              %QA for empty randomSequence data
                %In theory there might be reasons for this to be the case, but none are good for following analysis
              if ~batteryDesign && isempty( BLOCKS( thisFlyRowInd ).randomSequence ) && ( allowRandomSequenceEmpty == 0 && ~unsiphonedSEs )
                  ['## Alert: Failure to acquire/generate randomSequence data for block ',num2str( BLOCKS(thisFlyRowInd).blockNum ),' ##']
                  crash = yes
              end
    
        end
    
        %% Save a video
        if doVid && hasData
            %Make a bootleg 'full' form of the sequence across imaging (accounting for pre/post periods)
            
            %Old
            %fullStimTerp = nan( 1, size(thisImData,3) );
            %fullStimTerp( adjImStimStartVol:adjImStimEndVol  ) = imStimTerp;
            %Make a copy of the frames
            %imCopy = avg_z_green_aligned;
            %imCopy = thisImData;
    
            %New
            if bendyBlockDesign == 0
                imCopy = dataStimTrim;
                stimSeq = randomSeqActual;
            else
                imCopy = postStimData;
                stimSeq = deRandomSeq;
            end        
    
            %Use the sequence data to burn in some sequence representors
            %for unI = unique(medBTSeqInterpZ) %Old
            %for unI = unique(imStimTerp) %New
            for unI = unique(stimSeq) %Newer
                %thisSeqFrameCoords = find( medBTSeqInterpZ == unI ); %Old, assumption of stim==im
                %thisSeqFrameCoords = find( fullStimTerp == unI ); %New, accounts for pre/post non-stim time
                thisSeqFrameCoords = find( stimSeq == unI ); %New, accounts for pre/post non-stim time
                %imCopy( 1:12, 1 + (unI-1)*12:unI*12, thisSeqFrameCoords ) = 255*ones( 12,12, size( thisSeqFrameCoords,2 ) ); %Make little white boxes
                imCopy( 1:12, (unI + 1)*12:(unI+2)*12-1, thisSeqFrameCoords ) = 255*ones( 12,12, size( thisSeqFrameCoords,2 ) ); %Make little white boxes
            end
            imCopy = repmat( imCopy, 1, 1, 1, 3 ); %Add a fourth 'colour' dimension (Fake)
            imCopy = reshape( imCopy, size( imCopy,1 ), size( imCopy,2 ), size( imCopy,4 ), size( imCopy,3 ) ); %Reshape
            %Make a sacrificial array
            vidFrames2 = uint8( repmat( [0] , [size(imCopy, [1:2] ) , 3]) ); %Make a single frame of zeroes as an unsigned 8-bit
            vidFrames2 = repmat( vidFrames2 , 1 , 1 , 1 , size(imCopy,4) ); %Replicate to be full frame number
            vidFrames2( 1:size(imCopy,1), 1:size(imCopy,2), 1:3, : ) = imCopy; %Paste in original vid data
            %Save
            %vidOutObj = VideoWriter([strcat(dataFolder, filesep, expName, filesep, 'rollVid','.mp4')],'MPEG-4')
            vidOutObj = VideoWriter([strcat(dataFolder, filesep, expName, filesep, 'rollVid','.avi')], 'Motion JPEG AVI');
            %vidOutObj.FrameRate = size(avg_z_green_aligned,3) / inferTimes( frameEndInd ); %Estimate framerate
            %vidOutObj.FrameRate = size(thisImData,3) / inferTimes( frameEndInd ); %Estimate framerate
            vidOutObj.FrameRate = size(thisImData,3) / tsLastFrameTime; %Estimate framerate
            vidOutObj.Quality = 100 %Arbitrary
            tic
            open(vidOutObj)
            writeVideo(vidOutObj,vidFrames2);
            close(vidOutObj)
            clear imCopy vidFrames2
            disp(['Written in ',num2str(toc),'s'])
        end
    
        clear dataStimTrim imStimTerp randomSeqActual nVol nStimuli %Just in case, for rolling
        clear postStimData deRandomSeq nomInter %For block
    
    %block iteration end    
    end
    
    %Remove certain blocks, if requested
        %A little slower to do this here, rather than skipping analysis, but preserves order of blocks till end
    if disregardRollingDesign || disregardBattery
        for b = size( BLOCKS ,2 ):-1:1
            if ( disregardRollingDesign && (isfield( BLOCKS, 'bendyBlockDesign' ) && ~isempty(  BLOCKS( b ).bendyBlockDesign  ) && BLOCKS( b ).bendyBlockDesign == 0 )) || ...
                    ( disregardBattery && (isfield( BLOCKS, 'stimulus' ) && ~isempty(  BLOCKS( b ).stimulus  ) && isequal(BLOCKS( b ).stimulus, 'battery' ) ))
                disp([ '-# Block ', num2str( BLOCKS( b ).blockNum ), ' removed from analysis due to disregardation #-' ])
                BLOCKS( b ) = [];
            end
        end
    end

    if dataSource == -1 && ~simulationRun
        FLIES(fly).BLOCKS = BLOCKS;
        disp(['-- Modified BLOCKS reinserted into FLIES --'])
    elseif simulationRun
        ['-# Modified BLOCKS not reinserted, due to simulation #-']
    end
    
    disp(['-- Newtype data synchronised --'])

end

%function end
end