function BLOCKS = syncMaster(BLOCKS, flyRecord, options)
%Script/Function for synchronising newtype (2025+) 2p data synchronised with BT/ThorSync

%To do: Check support for follow-on analysis of true bendy rolling data

%To improve: Absolute best frame specificity for im/stim calcs
    %Specifically:
        %   - Bendy block design has some by-necessity inaccuracies with regards to stimulus aliasing (see stimCollInds/temp4)

    %Also, support for a bendy panel block design with blank periods (Currently all blanks are assumed to be inter-stim periods)

    %Also, ability to analyse LED data as if rolling?




arguments
    BLOCKS struct
    flyRecord table
    options.dataDirectory string = "I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\"
    options.doPlot double = 0
    options.doVid double = 0
    options.rollingAnalysis double = -1
    options.dataSource double = -1 %Whether to use data from function call (-1) or to look external for registered data (1)
    options.allowRandomSequenceEmpty double = 0 %Whether to allow randomSequence to end up empty
    options.disregardRollingDesign double = 0 %Whether to discard bendy rolling design blocks (e.g. When comparing LEDs to bendy block design)
end
%}
%{
BLOCKS = FLIES(fly).BLOCKS;
flyRecord = flyRecord;
options.dataDirectory = dataDirectory;
options.doPlot = 0
options.doVid = 0
options.rollingAnalysis = -1
options.dataSource = -1 %Whether to use data from function call (-1) or to look external for registered data (1)
options.allowRandomSequenceEmpty = 0
options.disregardRollingDesign = 1
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
rollingAnalysis = -1; %Make dynamic?
dataSource = options.dataSource;
allowRandomSequenceEmpty = options.allowRandomSequenceEmpty;
disregardRollingDesign = options.disregardRollingDesign;

for thisBlock = [BLOCKS]
    %thisFlyRowInd = find( flyRecord.Fly == thisBlock.flyNum & flyRecord.Block == thisBlock.blockNum ); %Wrong; Only valid with N=1 fly, due to flyRecord size
    thisFlyRowInd = find( [BLOCKS.flyNum] == thisBlock.flyNum & [BLOCKS.blockNum] == thisBlock.blockNum ); %Fixed; Will also be used later to overwrite data
    thisFlyRecord = flyRecord( thisFlyRowInd ,:)
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
    dataFolder = [strcat(options.dataDirectory, currentDate)];
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

    toc

    %Also attempt to acquire TS settings
    h5SetFile = dir([ hFiles(i).folder, filesep, 'ThorRealTimeDataSettings.xml' ]);
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

        %Read BT
        targetBTName = strcat( MATDate, '_', expNameMAT{1}, '_', expNameMAT{2},'_btData' );
        btFile =dir( strcat(dataFolder, filesep, '**/', targetBTName, '.csv') );
        %QA
        if isempty( btFile )
            ['## BT data not found! ##']
            crash = yes
        end
        btData = csvread( [btFile.folder,filesep,btFile.name] );

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

         estimatedPTBEndTime = imStartTime + inferTimes( bleachEndInd ); 
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


    %% Dynamic calcs
    guessIndex = [{'PTB'},{'Exp*.xml'}]; %List of ways to guess the start time
    guessMode = 0;
    if hasPTB
        bestGuessCommenceTime = ptbEndTime - inferTimes( bleachEndInd ); %Use PTB end as a known point in posix, match with inferred relative TS time
        guessMode = 1;
        %disp(['Best guess commencement time (from PTB): ',datestr(datetime( bestGuessCommenceTime , 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10'))])
    else
        bestGuessCommenceTime = imStartTime; %Use self-reported start time blindly
        guessMode = 1;
        %disp(['Best guess commencement time (from Exp*.xml): ',datestr(datetime( bestGuessCommenceTime , 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10'))])
    end
    disp(['Best guess commencement time (from ',guessIndex{guessMode},'): ',datestr(datetime( bestGuessCommenceTime , 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10'))])

    inferTimesPosix = inferTimes + bestGuessCommenceTime; %Posix form of inferTimes, for simplicity

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
    num2str( (inferTimes( frameEndInd )/60) ),' mins (',...
    num2str( (inferTimes( frameEndInd )) ),' s)'] )  
    %disp(['Estimated imaging end time (Exp self report -> Last frame high): ', ...
    %    datestr( datetime( imStartTime + (inferTimes( frameEndInd )), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])
    disp(['Estimated imaging end time (',guessIndex{guessMode},' -> Last frame high): ', ...
        datestr( datetime( bestGuessCommenceTime + (inferTimes( frameEndInd )), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])

    disp([ 'Estimated arduino duration: ', num2str((inferTimes( bleachEndInd ))),'s' ])
    %disp(['Estimated arduino end time (Exp self report -> Last bleach high): ', ...
    %    datestr( datetime( imStartTime + (inferTimes( bleachEndInd )), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])
    disp(['Estimated arduino end time (',guessIndex{guessMode},' -> Last bleach high): ', ...
        datestr( datetime( bestGuessCommenceTime + (inferTimes( bleachEndInd )), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])

    disp(['Estimated TS end time: ',...
        datestr( datetime( bestGuessCommenceTime + (inferTimes( end )), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC+10' ) ) ])

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
        bwFrameOut = bwlabel( syncStruct.DI.frameOutData );
        [~,frameOnsetIndices] = ismember( [1:nanmax(bwFrameOut)], bwFrameOut );
        tsImFrameCount = numel(frameOnsetIndices);
        disp(['ThorSync calculated number of imaging frames: ', num2str(tsImFrameCount)])
        %disp(['Expected number of frames: ',num2str(size( avg_z_green_aligned , 3 ) * (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames))])
        disp(['Expected number of frames: ',num2str(size( thisImData , 3 ) * (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames))])
        %QA
        if numel(frameOnsetIndices) ~= size( thisImData , 3 ) * (thisFlyRecord.Steps + thisFlyRecord.FlybackFrames) %Note: All of these vars need to be manually imported/defined
            ['-# Alert: TS framecount and imaging framecount differ #-']
            crash = yes %Eventually will probably have systems to handle this case
        end
    end

    %% Some more plots

    if hasData && doPlot
        %Imaging frame timing instability (courtesy of TS)
            %If there are many unique values (i.e. >2) on this graph, timing was unstable 
        figure
        hist( diff(frameOnsetIndices)/sampRate, 128 )
        title(['Inter-frame interval hist'])
        xlabel(['Time (s)'])
        ylabel(['Count'])
    end


    %% Tether image to known reference frame

    %To do: Non-linear frame placement based on TS frame timings
        %Probably not a huge issue, but likely to be useful
    if hasData && hasPTB
        %lastImStimFrameInd = find( frameOnsetIndices - bleachEndInd >= 0 , 1, 'first' ) - 1; %The last frame to receive a full stimulation; Arduino low inference
        lastImStimFrameInd = find( inferTimesPosix - ptbEndTime >= 0 , 1, 'first' ) - 1; %The last frame to receive a full stimulation; PTB posix inference
            %Note: May be at any point in a volume
        firstImStimFrameInd = find( inferTimesPosix - ptbStartTime >= 0 , 1, 'first' ); %The first imaging frame to receive stimulation
            %Again, can be any point
            %Note: Both of these in TS reference frame
        imStimStart = find( frameOnsetIndices - firstImStimFrameInd > 0 , 1 , 'first' ); %Might theoretically be better to find min, rather than first after?
        imStimEnd = find( frameOnsetIndices - lastImStimFrameInd > 0 , 1 , 'first' );
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

        %medBTSeqInterpZ = nanmedian( btSeqPosInterpZ , 1 ); %Time-matched (Imaging reference) list of stimulus present during that frame
        modeBTSeqInterpZ = mode( btSeqPosInterpZ , 1 ); %Switch to mode, so that can never be non-integer
            %E.g. pos 360 being 2 means frame 360 was being presented with stimulus #2
        %QA for critical aliasing issues
        %if numel(unique(medBTSeqInterpZ)) < numel( [medBTSeqInterpZ(1):medBTSeqInterpZ(end)] )
        if numel(unique(modeBTSeqInterpZ)) < numel( [modeBTSeqInterpZ(1):modeBTSeqInterpZ(end)] )
                %In theory detects if more elements likely existed than present in imaging frames
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
         imStimTerp = randomSeqCorr( modeBTSeqInterpZ );
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
              modeBTSeqInterpZ( baseInds ) = [];

              disp(['Blank pre-experiment baseline of ',num2str(numel(baseInds)),' elements collected; Data/sequence trimmed'])
          end
            

          hasSiphoned = 0; %Flag to indicate whether blanks siphoned off
          stillRollable = 1; %Flag to indicate whether data can still be analysed rolling-style after blank removal
          if any( imStimTerp == 5 )
              
              if bendyBlockDesign == 0 %"Blanks indicate true blank periods, intended for alternative analysis/etc"
                    disp(['Rolling design; Siphoning blanks for alternative analysis'])
    
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

                    %melting

                    %Assemble stimulus periods for collection
                        %Note: This relies on the assumption that there were no huge issues in theoretical->actual stimulus presentation
                    stimInds = collInds(:,1) - repmat( fliplr([1:blockLength]), size(collInds,1), 1 ); %Identify stimulated frame periods (Not for collection, but to derive sequence)
                        %Use collInds, not intInds, because cleaned
                        %Note: Due to frame attribution inaccuracies/phase differences between imaging and display, there will be variability in number of 'stimulus' frames preceding an inter-stimulus period
                    %stimSeq = randomSeqCorr( medBTSeqInterpZ( stimInds(:,:) ) ); %Pulls duplicates
                    %stimSeq = randomSeqCorr( [-blockLength+1:0] +  medBTSeqInterpZ( stimInds(:,end) ) );
                    stimSeq = nan( size(stimInds) );
                    stimSeqCorrInds = nan( size(stimInds) );
                    for row = 1:size( stimInds,1 )
                        %stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  medBTSeqInterpZ( stimInds(row,end) ); %Which elements of randomSeqCorr were pulled
                        stimSeqCorrInds( row, : ) =  [-blockLength+1:0] +  modeBTSeqInterpZ( stimInds(row,end) );
                            %pull medBTSeqInterpZ values at stimInds, go 5 back, pull those elements of randomSeqCorr
                        stimSeq( row, : ) = randomSeqCorr( stimSeqCorrInds( row, : ) ); %Stimuli
                    end
                    %QA
                    if any( stimSeq == 5 )
                        ['## Alert: Blanks collected in stimSequence; Probable phase failure ##'] %phailure
                        crash = yes
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
                  BLOCKS( thisFlyRowInd ).randomSequence = imStimTerp;
                  clear dataStimTrim imStimTerp %Just in case
              else %Block
                  BLOCKS( thisFlyRowInd ).greenChannel = postStimData;
                  BLOCKS( thisFlyRowInd ).randomSequence = deRandomSeq;
                  BLOCKS( thisFlyRowInd ).nVol = nomInter;
                  BLOCKS( thisFlyRowInd ).blankBlocks = 0; %Need to add support later for blank blocks
                  BLOCKS( thisFlyRowInd ).fauxBlockDesign = 1; %Just to keep track
                  disp(['Faux-block design created'])
                  clear postStimData deRandomSeq nomInter
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
        %fullStimTerp = nan( 1, size(avg_z_green_aligned,3) );
        fullStimTerp = nan( 1, size(thisImData,3) );
        fullStimTerp( adjImStimStartVol:adjImStimEndVol  ) = imStimTerp;

        %Make a copy of the frames
        %imCopy = avg_z_green_aligned;
        imCopy = thisImData;
        %Use the sequence data to burn in some sequence representors
        %for unI = unique(medBTSeqInterpZ) %Old
        for unI = unique(imStimTerp) %New
            %thisSeqFrameCoords = find( medBTSeqInterpZ == unI ); %Old, assumption of stim==im
            thisSeqFrameCoords = find( fullStimTerp == unI ); %New, accounts for pre/post non-stim time
            imCopy( 1:12, 1 + (unI-1)*12:unI*12, thisSeqFrameCoords ) = 255*ones( 12,12, size( thisSeqFrameCoords,2 ) ); %Make little white boxes
        end
        imCopy = repmat( imCopy, 1, 1, 1, 3 ); %Add a fourth 'colour' dimension (Fake)
        imCopy = reshape( imCopy, size( imCopy,1 ), size( imCopy,2 ), size( imCopy,4 ), size( imCopy,3 ) ); %Reshape
        %Make a sacrificial array
        vidFrames2 = uint8( repmat( [0] , [size(imCopy, [1:2] ) , 3]) ); %Make a single frame of zeroes as an unsigned 8-bit
        vidFrames2 = repmat( vidFrames2 , 1 , 1 , 1 , size(imCopy,4) ); %Replicate to be full frame number
        vidFrames2( 1:size(imCopy,1), 1:size(imCopy,2), 1:3, : ) = imCopy; %Paste in original vid data
        %Save
        %vidOutObj = VideoWriter([strcat(dataFolder, filesep, expName, filesep, 'rollVid','.mp4')],'MPEG-4')
        vidOutObj = VideoWriter([strcat(dataFolder, filesep, expName, filesep, 'rollVid','.avi')], 'Motion JPEG AVI')
        %vidOutObj.FrameRate = size(avg_z_green_aligned,3) / inferTimes( frameEndInd ); %Estimate framerate
        vidOutObj.FrameRate = size(thisImData,3) / inferTimes( frameEndInd ); %Estimate framerate
        tic
        open(vidOutObj)
        writeVideo(vidOutObj,vidFrames2);
        close(vidOutObj)
        clear imCopy vidFrames2
        disp(['Written in ',num2str(toc),'s'])
    end
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

disp(['-- Newtype data synchronised --'])

end
     