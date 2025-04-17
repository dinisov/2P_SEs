function FLIES = syncMaster(FLIES, flyRecord, options)
%function BLOCKS = syncMaster(BLOCKS, flyRecord, options)
%Script/Function for synchronising newtype (2025+) 2p data synchronised with BT/ThorSync

%To do: Check support for follow-on analysis of true bendy rolling data

%           - Support for blanks

%           - Support for battery


%To improve: Absolute best frame specificity for im/stim calcs
    %Specifically:
        %   - Bendy block design has some by-necessity inaccuracies with regards to stimulus aliasing (see stimCollInds/temp4)

    %Also, support for a bendy panel block design with blank periods (Currently all blanks are assumed to be inter-stim periods)

    %Also, ability to analyse LED data as if rolling?



%{{
arguments
    %BLOCKS struct
    FLIES struct
    flyRecord table
    options.dataDirectory string = "I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\"
    options.doPlot double = 0
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
end
%}
%{
%BLOCKS = FLIES(fly).BLOCKS;
FLIES = FLIES;
flyRecord = flyRecord;
options.dataDirectory = dataDirectory;
options.doPlot = 0
options.doVid = 0
options.rollingAnalysis = -1
options.dataSource = -1
options.allowRandomSequenceEmpty = 0
options.disregardRollingDesign = 0
options.btInterpolationMethod = 'intelligent' 
options.nBack = 5
options.saveShortcut = 1
options.useShortcut = 0
options.daqFramespikeVoltage = 1
%}

%{
clear %Will go away once functionised
close all %Will go away once functionised

%% Prepare

doPlot = 1;
doVid = 0;
rollingAnalysis = -1; %Whether to force rolling analysis (1) or not (0), or go with detection from PTB files if available (default to no)

%Load block data (Unnecessary if functionised)
blocks = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");
blocks = blocks(~logical(blocks.Exclude),:);

chosenFly = 260;
chosenBlock = 2;
%}
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
        if useShortcut == 1
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
        end
    
        if useShortcut == 0 || isShortcutting == 0
        %Read ThorSync
        %fileName = "I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\30Jan25\TS\fly1_30Jan25_exp2\Episode001.h5"
        fileName = hTarget;
    
        %h5disp(fileName)
    
        syncStruct = struct;
    
        disp(['-- Loading H5 data --'])
        tic
    
        syncStruct.AI.piezoData = h5read(fileName, '/AI/PiezoMonitor');
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
        catch
            ['-# Could not retrieve DAQ-related TS data #-']
            hasDaqData = 0;
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
        %QA
        if isempty( matParamFile )
            ['## PTB MATLAB parameters not found! ##']
            %crash = yes
            hasPTB = 0;
        else
            hasPTB = 1;
        end
        if hasPTB
            matParamStruct = load( [matParamFile.folder,filesep,matParamFile.name] );
    
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
            if nansum( isnan( btData(end,:) ) ) == numel( btData(end,:) )
                btData(end,:) = [];
            end
            %QA
            if nansum( isnan(btData), 'all' ) > 0
                ['-# Warning: NaN elements found in btData #-']
                crash = yes %Again, not critical but problematic if not at end
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
    
            ptbStartTime = btData(1,2)/1000; %Divide by 1000 because PTB stores down to millisecond unix time
            %ptbEndTime = btData( find( ~isnan( btData(:,2) ) , 1, 'last'  ) ,2)/1000; %Find necessary because NaN last element
            ptbEndTime = btData( [end-find( diff( flip( btData(:,4) )) ~= 0, 1, 'first')] ,2)/1000; %More accurate, use last non-repetitious time
                %Estimated accuracy: Effectively 100% (~millisecond), no drift
                    %Note: Datetime/etc does not show milliseconds, but this value is (seemingly) still stored
                        %i.e. There is a ~0.4s difference between the last non-repetitious and the 'end' of the data
                    %Secondary note: Relies heavily on the assumption that arduino low -> bt repetitiveness instantly
    
             if ~isShortcutting        
             estimatedPTBEndTime = imStartTime + inferTimes( bleachEndInd ); 
             else
             estimatedPTBEndTime = shortStruct.estimatedPTBEndTime;
             end
                %This is a time guessed from the self-reported imaging start time
                %It is intended to be used to QA ptbEndTime in case of not ending because repetitive, not a true comparison
              %QA
              if abs( ptbEndTime - estimatedPTBEndTime ) > 0.05 * ( estimatedPTBEndTime - imStartTime ) %"Did PTB end more than 5% +- Arduino time?"
                  ['## Alert: PTB self-reported end time differs significantly from Arduino estimated end time ##']
                  crash = yes %Not technically critical, but probably worrying
              end
    
              %Find out if rolling experiment
              if rollingAnalysis == -1 && isfield( matParamStruct.matSave, 'blockDesign' )
                  if matParamStruct.matSave.blockDesign == 1
                      rollingAnalysis = 0; %Override
                      disp(['Block design detected in PTB params; Using block analysis'])
                  else
                      rollingAnalysis = 1;    
                      disp(['Rolling design detected in PTB params; Using rolling analysis'])
                  end
              end
    
        end

        %------------------------------

        %DAQ
        if hasDaqData && ~isShortcutting
            %freak
            %Find framespike peaks
                %NOTE: Current system optimised towards local peak of framespike, whereas 'true' timing of stimulus onset probably lies at onset of voltage increase
            if hasPTB && ~isequal( matParamStruct.matSave.stimuli, 'battery' )
                [framePKS, frameLOCS] = findpeaks( syncStruct.AI.FrameSpike, 'MinPeakHeight', 0.5*daqFramespikeVoltage,...
                   'MinPeakDistance', ( sampRate / matParamStruct.matSave.frequency )*0.66 ); %Look for peaks separated by at 66% of a flip
            else
                [framePKS, frameLOCS] = findpeaks( syncStruct.AI.FrameSpike, 'MinPeakHeight', 0.5*daqFramespikeVoltage, ...
                    'MinPeakWidth', 200); %Use an empirical spike width based on pilot data
            end
            nSpikes = size(frameLOCS,2);
            disp(['Found ',num2str(nSpikes),' apparent peaks in DAQ framespike data'])
            disp(['(Approximate peak/presumed flip rate: ',...
                num2str(1 / ( nanmean( diff( frameLOCS ) ) / sampRate )),'Hz w/ ~',...
                num2str(nanstd( diff( frameLOCS ) ) / sampRate),'s SD)'])
            if hasPTB
                disp(['(Expected fliprate: ',num2str(matParamStruct.matSave.frequency),'Hz)'])
            end

            %Match framespikes to btData
            if nSpikes == btData(end,5)-1 %-1 based on exactly N=1 testing
                disp(['Perfect match between number of detected framespikes (',num2str(nSpikes),') and reported last valid i value (',num2str(btData(end,5)-1),')'])
            
                temp = bwlabel( btData(:,9) == 2 ); %Assumption of having btData, but is this analysis possible without it?
                [~,flipOnsetIndices] = ismember( [1:nanmax(temp)], temp ); %First position of state switch in onOff
                    %Necessary to find first position in case of different duty cycles?
                %QA
                if length( flipOnsetIndices ) ~= nSpikes
                    ['## Alert: Disparity between detected framespikes and detected flip positions ##']
                    crash = yes
                end

                estPTBInferTimeStart = inferTimes( frameLOCS(1) ) - btData( flipOnsetIndices(1), 6); %Theoretical time (inferTimes reference) PTB started at
                    %Use framespike to find first flip position in inferTimes, then subtract known duration since PTB start from that
                [inferPTBStartDisp, inferPTBStartInd] = min( abs(inferTimes - estPTBInferTimeStart) );  %Find closest match in inferTimes for estimated PTB start time based on framespike
                    %Returns respectively the time error between the inferTimes point and PTB start, as well as the index of inferTimes (/TS) where PTB started (But not necessarily first stim)
                 
                estPTBInferTimeEnd = inferTimes( frameLOCS(end) ) + ( btData(end,6) - btData( flipOnsetIndices(end), 6) ); 
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
        %Little follow-on
        if isShortcutting
            hasDaqData = shortStruct.hasDaqData;
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
        disp(['Best guess commencement time (from ',guessIndex{guessMode},'): ', ...
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
    
        if doPlot
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
            else
                frameOnsetIndices = shortStruct.frameOnsetIndices;
                tsImFrameCount = shortStruct.tsImFrameCount;
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
    
            shortStruct.tsLastFrameTime =tsLastFrameTime;
            shortStruct.tsLastArduinoTime = tsLastArduinoTime;
            shortStruct.tsLastTime = tsLastTime;
    
            if hasPTB %&& hasData 
                shortStruct.lastImStimFrameInd = lastImStimFrameInd;
                shortStruct.firstImStimFrameInd = firstImStimFrameInd;
            end

            shortStruct.hasDaqData = hasDaqData;
    
            %save( [hFolder,filesep,'h5Shortcut.mat'], 'shortStruct' )
            %Check/Make folder to put shortcut files into
            shortFolder = strcat( dataFolder,filesep,'SHORT' );
            if exist(shortFolder) ~= 7
                mkdir( shortFolder )
                disp(['Shortcut folder made at ',shortFolder])
            end
            save( [shortFolder,filesep,shortStruct.h5Name,'_h5Shortcut.mat'], 'shortStruct' )
            %sword
    
        end
    
    
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
                ['-# Alert: Failure to identify imaging start OR end #-']
                if isShortcutting ~= 1
                ['Imaging started ', num2str( inferTimes( frameOnsetIndices(1) ) ),'s after TS start']
                ['Imaging ended apparently after ', num2str( inferTimes( frameOnsetIndices(end) ) / 60 ),'m']
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
            [~,adjImStimStartVol] = find( imFrameListZ == adjImStimStart );
            [~,adjImStimEndVol] = find( imFrameListZ == adjImStimEnd );
            disp(['Adjusted imaging volumes of stimulation interest: ', num2str(adjImStimStartVol),' : ', num2str(adjImStimEndVol),...
                ' (',num2str(adjImStimEndVol-adjImStimStartVol+1),' total)' ])
    
            %Interpolate sequence element across (usable) volume space
            %Old
            %{
            btSeqPosInterp = interp1([1:size(btData,1)], btData(:,5) , linspace(1,size(btData,1), tsImFrameCount ), 'previous' )'; %Interpolate btData 5th column to be imaging size
            btSeqInterp = randomSequence( btSeqPosInterp ); %Grab appropriate randomSequence elements for the timing
            btSeqInterpZ = reshape( btSeqInterp, (thisFlyBlock.Steps + thisFlyBlock.FlybackFrames), size( avg_z_green_aligned , 3 )  ); %Reshape as above
                %Basically represents for each volume the stimulus being presented (Changes mid-column represent the stim changing mid z-stack)
            %}
            %New
            btSeqPosInterpZ = nan( size(imFrameListZ,1), size(imFrameListZ,2) ); %Make empty array for randomSequence to go into, of *total* imaging size
            btSeqPosInterpZ( imStimStart: imStimEnd) = interp1([1:size(btData,1)], btData(:,5) , linspace(1,size(btData,1), imStimEnd-imStimStart+1 ), 'previous' )'; %Interpolate and transplant randomSequence into 'frame' list
            btSeqPosInterpZ = btSeqPosInterpZ( : , adjImStimStartVol:adjImStimEndVol );
    
            %Pick a method to identify preceding stimuli by
            disp(['Interpolating bt sequence'])
            switch btInterpolationMethod
                case 'median'
                    %medBTSeqInterpZ = nanmedian( btSeqPosInterpZ , 1 ); %Time-matched (Imaging reference) list of stimulus present during that frame
                        %E.g. pos 360 being 2 means frame 360 was being presented with stimulus #2
                    methodBTSeqInterpZ = nanmedian( btSeqPosInterpZ , 1 );
                    disp(['Interpolated with median'])
                case 'mode'
                    %modeBTSeqInterpZ = mode( btSeqPosInterpZ , 1 ); %Switch to mode, so that can never be non-integer
                    methodBTSeqInterpZ = mode( btSeqPosInterpZ , 1 );
                    disp(['Interpolated with mode'])
                case {'last' , 'intelligent'} %"aka, yaoi"
                    methodBTSeqInterpZ = btSeqPosInterpZ( end , : );
                    if isnan( methodBTSeqInterpZ(end) )
                        methodBTSeqInterpZ(end) = []; %May cause desync?
                    end
                    if any( isnan(methodBTSeqInterpZ) )
                        ['## Alert: Excessive NaN presence in last-interpolated bt sequence ##']
                        crash = yes
                    end
                    disp(['Interpolated with last OR intelligent'])
            end
            %QA for critical aliasing issues
            %if numel(unique(medBTSeqInterpZ)) < numel( [medBTSeqInterpZ(1):medBTSeqInterpZ(end)] )
            %if numel(unique(modeBTSeqInterpZ)) < numel( [modeBTSeqInterpZ(1):modeBTSeqInterpZ(end)] )
            if numel(unique(methodBTSeqInterpZ)) < numel( [methodBTSeqInterpZ(1):methodBTSeqInterpZ(end)] )
                    %In theory detects if more elements likely existed than present in imaging frames
                        %The most common cause for this is imaging not being (much) faster than stimulation
                ['-# Warning: Sequence elements apparently lost due to aliasing/imaging framerate #-']
            end
    
             %Find actual stimuli delivered (via randomSequence)   
             randomSeqOrig = matParamStruct.matSave.randomSequence; %Note that PTB script saves randomSequence as -1 (Blank), 1, (Stim 1), and 2 (Stim 2), rather than 5, 0, and 1, respectively
             if isequal( unique( randomSeqOrig ), [-1, 1, 2] ) || isequal( unique( randomSeqOrig ), [1, 2] )
                 randomSeqCorr = randomSeqOrig;
                 randomSeqCorr( randomSeqCorr == -1 ) = 5;
                 randomSeqCorr( randomSeqCorr == 1 ) = 0;
                 randomSeqCorr( randomSeqCorr == 2 ) = 1; %Note order as critical, lest all 2 -> 1 -> 0
             else
                 ['## Alert: Unrecognised elements present in randomSequence! ##']
                 unique( randomSeqOrig )
                 crash = yes
             end
    
             %imStimTerp = randomSeqCorr( medBTSeqInterpZ );
             %imStimTerp = randomSeqCorr( modeBTSeqInterpZ );
             imStimTerp = randomSeqCorr( methodBTSeqInterpZ );
                %Note: Due to volume/timing inefficiencies this may skip the first actual presented element/etc
                    %Also, number of imaging frames/volumes per element/block may be inconsistent
    
              %Trim data to relevant (i.e. Stimulated) portion
              %dataStimTrim = avg_z_green_aligned( :,:, adjImStimStartVol:adjImStimEndVol );
              dataStimTrim = thisImData( :,:, adjImStimStartVol:adjImStimEndVol );
              disp(['Data trimmed to stimulated portion only'])
    
              %Siphon out blanks
                    %Start by extracting pre-experiment blank baseline (If existing)
              if find( imStimTerp ~= 5, 1, 'first' ) ~= 1 && imStimTerp(1) == 5
                  blankBaseline = [];
                  baseInds = [ 1:find( imStimTerp ~= 5, 1, 'first' )-1 ]; %Note: Assumption that there is no pre-baseline of stimuli/etc
                  blankBaseline = dataStimTrim( :,:,baseInds );
                  blankSeq = imStimTerp( baseInds );
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
    
                  disp(['Blank pre-experiment baseline of ',num2str(numel(baseInds)),' elements collected; Data/sequence trimmed'])
              end
    
              %if thisBlock.blockNum == 3
              %prequel
              %end
    
              %Perform bendy rolling specific processing
              if bendyBlockDesign == 0
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
                  disp( ['Bendy rolling data reindicised to length ',num2str( size(dataStimTrim,3) ),' (',num2str( size(dataStimTrim,3) / nVol ),' "blocks") for analysis (',...
                      num2str( ( 1 - ( size(dataStimTrim,3) / max(newFrameInds) ) ) * 100 ),'% Loss)'] )
                  %Loss in this context is how much of the original imaging data could not be associated with 'blocks' 
                    %It is likely to be highest when imaging rates are low and stimulus frequencies are high (i.e. 1 in 3 is 33%, but 1 in 8 is only 12.5%)
    
              end
    
              hasSiphoned = 0; %Flag to indicate whether blanks siphoned off
              stillRollable = 1; %Flag to indicate whether data can still be analysed rolling-style after blank removal
              if any( imStimTerp == 5 ) || bendyBlockDesign == 1
                  if bendyBlockDesign == 0 %"Blanks indicate true blank periods, intended for alternative analysis/etc"
                    %Note: This isn't really main rolling analysis here, just blank removal, unlike below for block design bendy
                        disp(['Rolling design; Siphoning blanks for alternative analysis'])
                        not coded for randomSeqActual yet
                        especially in context of reindicising
        
                        %blankStack = dataStimTrim( imStimTerp == 5 );
                        [inds] = find( imStimTerp == 5 );
                            %Note: For rolling bendy design, this can be thought of as an actual 'blank' period
                            %          For bendy block design, blanks are actually the analagous recording period
                        blankStack = dataStimTrim( :,:, inds );
                        blankInds = inds;
            
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
    
                        %Prepare a labelled form of randomSeqCorr for intelligent stimuls collection, if applicable
                        if isequal( btInterpolationMethod, 'intelligent' ) 
                            temp = randomSeqCorr;
                            temp( temp == 0 ) = 1;
                            temp( temp == 5 ) = 0;
                            randomSeqLabel = bwlabel( temp );
                        end
    
                        %melting
    
                        %Assemble stimulus periods for collection
                            %Note: This relies on the assumption that there were no huge issues in theoretical->actual stimulus presentation
                        stimInds = collInds(:,1) - repmat( fliplr([1:blockLength]), size(collInds,1), 1 ); %Identify stimulated frame periods (Not for collection, but to derive sequence)
                            %Use collInds, not intInds, because cleaned
                            %Note: Due to frame attribution inaccuracies/phase differences between imaging and display, there will be variability in number of 'stimulus' frames preceding an inter-stimulus period
                        %Quick QA
                        if nanmax( nanmax( stimInds ) ) > methodBTSeqInterpZ(end)
                                %Not 100% sure of correct reference frame here
                            ['## Alert: Requested post-stimulus imaging frame beyond apparent imaging duration ##']
                            crash = yes
                                %This is marginally more likely with last due to lastNaN removal
                        end
                        %stimSeq = randomSeqCorr( medBTSeqInterpZ( stimInds(:,:) ) ); %Pulls duplicates
                        %stimSeq = randomSeqCorr( [-blockLength+1:0] +  medBTSeqInterpZ( stimInds(:,end) ) );
                        stimSeq = nan( size(stimInds) );
                        stimSeqCorrInds = nan( size(stimInds) );
                        for row = 1:size( stimInds,1 )
                           % burnthem
                           switch btInterpolationMethod
                               case {'median' , 'mode', 'last'}
                                    %stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  medBTSeqInterpZ( stimInds(row,end) ); %Which elements of randomSeqCorr were pulled
                                    %stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  modeBTSeqInterpZ( stimInds(row,end) );
                                    stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  methodBTSeqInterpZ( stimInds(row,end) );
                                        %pull medBTSeqInterpZ values at stimInds, go 5 back, pull those elements of randomSeqCorr
                                    stimSeq( row, : ) = randomSeqCorr( stimSeqCorrInds( row, : ) ); %Stimuli
                               case 'intelligent'
                                   thisStimEnd =  methodBTSeqInterpZ( stimInds(row,end) );
                                   %Check for accidentally landing in inter-stimulus period
                                   if randomSeqLabel( thisStimEnd ) == 0
                                       if randomSeqLabel( thisStimEnd-1 ) ~= 0
                                           thisStimEnd = thisStimEnd - 1; %Adjust                                   
                                       else
                                            ['## Fatal failure to intelligently identify stimulus block identity ##']
                                            crash = yes
                                            %Being 1 stimulus before the end of the block = okay, being 1 after = acceptable, but any more probably indicates a large error has occurred
                                            %Likely some grand loss of phase
                                       end
                                   end
                                   randomSeqBlock = randomSeqLabel( thisStimEnd );
                                   thisStimInds = [ find( randomSeqLabel == randomSeqBlock,1, 'first' ) : find( randomSeqLabel == randomSeqBlock,1, 'last' ) ];
                                   %QA
                                   if numel(thisStimInds) ~= blockLength
                                       ['## Fatal mismatch between identified stimulus block length and theoretical block length ##']
                                       crash = yes
                                   end
                                   %If no issues, use labelled region to get stimuli
                                   stimSeq( row, : ) = randomSeqCorr( thisStimInds );
                           end
                        end
                        %QA
                        if any( stimSeq == 5, 'all' )
                            ['## Alert: Blanks collected in stimSequence; Probable phase failure ##'] %phailure
                            crash = yes
                                %If using interpolation other than last, this may be result of median/mode 'skipping' stimulus elements 
                                    %i.e. >1 stimuli occurred in the span of a single volume
                        end
    
                        %Reshape derived stimulus sequence
                        deRandomSeq = reshape(stimSeq', 1, size(stimSeq,1)*blockLength); %Transposition v important here for correct rowwise-ness
                         
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
                        end
    
                        %Report
                        disp(['Final number of imaging events: ',num2str(size( postStimData,3 )/nomInter)]) %Add potential max # imaging events here
                        disp(['Final number of stimulus events: ',num2str(size( deRandomSeq,2 )/blockLength)]) %Add potential total (btData) here
    
                  end
              end
    
    
              %Overwrite data (if function source)
              if dataSource == -1
                  %thisBlock.greenChannel = dataStimTrim; %Mostly unnecessary
                  %thisBlock.randomSequence = imStimTerp; %Mostly unnecessary
                  %BLOCKS(thisFlyRowInd,:) = thisBlock;
                  if bendyBlockDesign == 0 %Rolling
                      BLOCKS( thisFlyRowInd ).greenChannel = dataStimTrim;
                      %BLOCKS( thisFlyRowInd ).randomSequence = imStimTerp;
                      BLOCKS( thisFlyRowInd ).randomSequence = randomSeqActual; %New, not interpolated
                        %Note: Using randomSequence not of exact same length as imaging may cause issues with rolling implementation in analyseBlock
                      BLOCKS( thisFlyRowInd ).nVol = nVol;
                      BLOCKS( thisFlyRowInd ).nStimuli = nStimuli;
                      %clear dataStimTrim imStimTerp randomSeqActual nVol nStimuli %Just in case
                  else %Block
                      BLOCKS( thisFlyRowInd ).greenChannel = postStimData;
                      BLOCKS( thisFlyRowInd ).randomSequence = deRandomSeq;
                      BLOCKS( thisFlyRowInd ).nVol = nomInter;
                      BLOCKS( thisFlyRowInd ).blankBlocks = 0; %Need to add support later for blank blocks
                      BLOCKS( thisFlyRowInd ).fauxBlockDesign = 1; %Just to keep track
                      disp(['Faux-block design created'])
                      %clear postStimData deRandomSeq nomInter
                  end
                  if isfield( matParamStruct.matSave, 'blockDesign' )
                     BLOCKS( thisFlyRowInd ).bendyBlockDesign =  matParamStruct.matSave.blockDesign;
                  end
                  disp(['Modified data and randomSequence inserted into BLOCKS'])
                  %{
                    %Blank data not currently valid for insertion into BLOCKS
                  if hasSiphoned == 1 && bendyBlockDesign == 0
                      BLOCKS( thisFlyRowInd ).blankImageStack = blankStack; %Note: Architecture ostensibly should be repeating, 'nVol' sized groups of frames                  
                      BLOCKS( thisFlyRowInd ).blankImageInds = blankInds;
                      disp(['Modified blank data (and inds) inserted into BLOCKS'])
                  end
                  %}
                  if stillRollable == 0
                      BLOCKS( thisFlyRowInd ).isRolling = 0;                
                      disp(['Rolling status revoked'])
                  end
                  BLOCKS( thisFlyRowInd ).syncModified = 1;
              end
    
              %QA for empty randomSequence data
                %In theory there might be reasons for this to be the case, but none are good for following analysis
              if isempty( BLOCKS( thisFlyRowInd ).randomSequence ) && allowRandomSequenceEmpty == 0
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
    if disregardRollingDesign == 1
        for b = size( BLOCKS ,2 ):-1:1
            if isfield( BLOCKS, 'bendyBlockDesign' ) && ~isempty(  BLOCKS( b ).bendyBlockDesign  ) && BLOCKS( b ).bendyBlockDesign == 0
                disp([ '-# Block ', num2str( BLOCKS( b ).blockNum ), ' removed from analysis due to disregardation of rolling design #-' ])
                BLOCKS( b ) = [];
            end
        end
    end

    if dataSource == -1
        FLIES(fly).BLOCKS = BLOCKS;
        disp(['-- Modified BLOCKS reinserted into FLIES --'])
    end
    
    disp(['-- Newtype data synchronised --'])

end

%function end
end