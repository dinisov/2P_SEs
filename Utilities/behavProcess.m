%Utility for 

close all
clear

vidDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\Matt\Videos'; %Note lack of terminal backslash
    %Based on vidDirectory
    
blocks = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");

chosenFlies = [222:224];
    %Note: Corresponding vids/CSVs must be terminated with "mov" and start with "fly<number>" where <number> is the fly number on the day 

%%

syncWithSequence = 1; %Whether to synchronise this data with available sequence data
if syncWithSequence == 1
    mainDirectory = '\\uq.edu.au\uq-inst-gateway1\\RFDG2021-Q4413\2P_Data\Gcamp7s_CC';
        %Formerly dataDirectory
    sequenceDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\Matt\Data_LEDs';
    durationMethod = 2; %Whether to use 1 - Experiment frames/flyback/etc parameters (Least accurate), 2 - tif file DateModified (Better accurate),
                                        %or 3 - RaspberryPi theoretical experiment calculations (Superior accuracy)[Not implemented yet], to calculate recording duration
end

doAcInac = 1;
if doAcInac == 1
    acMeanThresh = +2; %How far relative to mean activity must be (in SD) to be classified as "active"
    minAcTime = 20; %Min. time in seconds that activity must be separated by for it to be classified as two 'pockets' of activity
        %Based on empirical data showing majority of ac pockets separated by <18s of inactivity (at mean+2*sd)
        %Note: This format (s) is a departure from current state of SOFAS use
        %Set to 0 to effectively disable
    minSleepBoutTime = 1*60; %How many seconds minimum for an inactivity region to be classified as 'sleep'
        %Note: This is applied after minAcTime
end

saveOutput = 1; %Whether to actually save

%%

failCount = 0;

flagParamSaveList = who;
flagParamSaveList = [flagParamSaveList; 'flagParamSaveList'; 'fly'];

for fly = 1:length(chosenFlies)
    clearvars('-except', flagParamSaveList{:})
    
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    flyOnDay = thisFlyBlocks.FlyOnDay(1);
    
    thisVidFolder = [vidDirectory,filesep,currentDate];

    disp([char(10),'-- Now processing: #',num2str(chosenFlies(fly)),' - ',currentDate,' - fly ',num2str(flyOnDay),' --'])
    disp(thisVidFolder)
    tic
    
    thisCSVs = dir( [thisVidFolder,filesep,'fly',num2str(flyOnDay),'*mov.csv'] );
    %QA
    if isempty(thisCSVs)
        ['-# No data found for ',num2str(chosenFlies(fly)),' #-']
        failCount = failCount + 1;
        continue
    end
    
    %Collate
    allData = [];
    sepData = [];
    for IIDN = 1:size(thisCSVs,1)
        %allData = [allData; readtable([thisCSVs(IIDN).folder,filesep,thisCSVs(IIDN).name]) ];
        sepData{IIDN} = readtable([thisCSVs(IIDN).folder,filesep,thisCSVs(IIDN).name]);
        allData = [allData;...
            readtable([thisCSVs(IIDN).folder,filesep,thisCSVs(IIDN).name]), ...
            array2table( repmat(IIDN, size(sepData{IIDN},1), 1) ,'VariableNames',{'CSVnum'})];
    end

    %Convert to nice datettimes
    allDates = num2str([allData.Year,allData.Month,allData.Date,allData.Hour,allData.Mins,allData.Seconds,allData.usec]);
    %Clear up potentially missing proximal zeroes
    for wU = 5:-1:0 %uWu
        zeroInds = find(allDates(:,end-wU) == ' ');
        allDates(zeroInds,end-wU) = '0';
    end
    allDatetimes = datetime(allDates,'Format', 'yyyy      MM      dd      HH      mm      ss  SSSSSS');

    behavFrameRate = 1 / nanmedian(diff([allData.Seconds+(allData.usec/1000000)]));
        %Bootleg derive the framerate by finding the median time difference between frames

    if doAcInac == 1
        %Calculate sleep/wake
        disp(['-- Beginning activity/inactivity separation calculations --'])
        movData = allData.avCntrSize;
        movMean = nanmean(movData);
        movSD = nanstd(movData);
        
        acThreshLevel = movMean + acMeanThresh*movSD;
        
        acUpper = movData > acThreshLevel;
        
        %New, BWLabel based method
        %minActivityTime = minAcTime * BaseFrameRate;
        minActivityTime = minAcTime * behavFrameRate;
        tic
        %tempUpperBinary = isnan( acUpper ) ~= 1;
        %tempUpperBW = bwlabel( tempUpperBinary );
        tempUpperBW = bwlabel( ~acUpper ); %Label gaps between ac pockets
        %invTempUpperBW = bwlabel( acUpper ); %Label ac pockets
        
        temp = nansum( tempUpperBW == [1:nanmax(tempUpperBW)] , 1); %Find sizes of all gaps
        %temp2 = nansum( invTempUpperBW == [1:nanmax(invTempUpperBW)] , 1); %Find sizes of all pockets
        
        %temp2 = nansum(tempUpperBW == find(temp <= (minAcTime*behavFrameRate)),2);
        
        acUpperProc = acUpper;
        acUpperProc( nansum(tempUpperBW == find(temp <= (minAcTime*behavFrameRate)),2) == 1 ) = 1;
        %Find all gaps separated by <minAcTime, set to 1 in binary data
        
        inacBinaryProc = ~acUpperProc;
        inacBinaryBW = bwlabel( inacBinaryProc );
        temp = nansum( inacBinaryBW == [1:nanmax(inacBinaryBW)] , 1); %Find sizes of all inacBouts
        disp(['There are ',num2str(nanmax(inacBinaryBW)),' prospective inactivity bouts (Av: ',num2str(nanmean(temp/behavFrameRate)),'s)'])
        
        inacBinaryProc( nansum(inacBinaryBW == find(temp < (minSleepBoutTime*behavFrameRate)),2) == 1 ) = 0; 
        %Find all prospective sleep bouts < minimum size, flatten
        inacBinaryBW = bwlabel( inacBinaryProc ); %Recalculate
        
        %tempUpperBinary = tempUpperBW ~= 0; %Recalculate binary ac/inac
        %tempUpperBW = bwlabel( tempUpperBinary );
        %invTempUpperBW = bwlabel(~tempUpperBinary);  
        
        %Plot raw data/thresh
        figure
        plot( movData, 'k' )
        hold on
        line([0,length(movData)],[movMean + acMeanThresh*movSD,movMean + acMeanThresh*movSD], 'Color', 'r', 'LineStyle', '--')
        plot( normalize(allData.DeltaProp, 'range', [nanmin(movData),nanmean(movData)] ), 'color', [0.8,0.8,0.8] )
        title(['Activity (k) and calculated threshold (r) [mean+',num2str(acMeanThresh),'*SD]'])
        set(gcf,'Name',['Fly #',num2str(chosenFlies(fly)),' - ',currentDate,' - fly ',num2str(flyOnDay),' ac/inac separation'])
        
        %Plot binaries, for usefulness
        figure
        subplot(3,1,1)
        plot(acUpper,'k')
        title(['Original activity binary'])
        subplot(3,1,2)
        plot(acUpperProc, 'r')
        title(['Processed activity binary'])
        subplot(3,1,3)
        plot( inacBinaryProc ,'b')
        title(['Final sleep binary'])
        set(gcf,'Name',['Fly #',num2str(chosenFlies(fly)),' - ',currentDate,' - fly ',num2str(flyOnDay),' behav binaries'])
        
        %Save
        inStruct = struct;
        inStruct.holeSizes = nansum( inacBinaryBW == [1:nanmax(inacBinaryBW)] , 1);
        inStruct.holeStarts = [find( diff(inacBinaryProc) == 1 )+1]';
        inStruct.holeEnds =  find( diff(inacBinaryProc) == -1 )';
        %And cover some 'edge' cases (Note: Minor risk of introducing asynchrony between starts and ends)
        if inacBinaryProc(1) == 1
            inStruct.holeStarts = [1,inStruct.holeStarts];
        end
        if inacBinaryProc(end) == 1
            inStruct.holeEnds = [inStruct.holeEnds,length(inacBinaryProc)];
        end
        
        disp(['There were ',num2str(length(inStruct.holeSizes)),' sleep bouts meeting criteria (Min ',num2str(minSleepBoutTime),'s)'])
        
        toc
        
        %quetz
    end


    if syncWithSequence == 1
        %------------
        disp('Synchronising behavioural data with 2p sequence')
        
        nBlocks = height(thisFlyBlocks);
        overSave = struct;
        for b = 1:nBlocks
            currentBlock = thisFlyBlocks(b,:);
            flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate]; 
            currentDataDirectory = fullfile(mainDirectory,currentDate,flyID);
            %QA
            if isempty( dir(currentDataDirectory) ) == 1
                ['## Data not found! ##']
                %crash = yes
                continue
            end
            disp(flyID)
            
            %Synchronise with experiment params
            %(Stolen from Transient analysis)
            filename = [currentDataDirectory,filesep,'Experiment.xml'];
            %Acquire imaging parameters
            expParams = readstruct(filename); %Note: Only works with MATLAB >=2021b
            imagingStartTime = expParams.Date.dateAttribute;
            %Timing notes:
            %   - Python time can be calculated by multiplying enforced delays (3.02s not including overhead) by number of blocks (e.g. 1000) to give a relatively good estimate
            %   - However, subtracting the dateModified time of the ChanA_Preview.tif file from the experiment start time (ExpParams) seems to actually be more accurate
            %       (If overhead is assumed to exist)
            %   - Testing shows the former value to yield an example duration of 50m20s, with the latter guessing 50m31s
            %       (For comparison, purely calculating based on ExpParams attributes yields 51m28s)
            %   - Note: DateModified for all files *except* ChanA_Preview seem to be slaved to when the operator returned to the station, not the true end of recording
            %   - Empirical testing shows true duration to be ~50m20s for a standard 800/200 experiment

            if durationMethod == 1
                %Old, frame-based system
                %expDurSecs = ( expParams.Streaming.framesAttribute / expParams.LSM.frameRateAttribute ) +...
                %    ( expParams.LSM.flybackCyclesAttribute/expParams.ZStage.stepsAttribute*expParams.Streaming.framesAttribute / expParams.LSM.frameRateAttribute );
                expDurSecs = ( ( expParams.Streaming.framesAttribute / expParams.LSM.frameRateAttribute ) +...
                ( expParams.LSM.flybackCyclesAttribute/expParams.ZStage.stepsAttribute*expParams.Streaming.framesAttribute / expParams.LSM.frameRateAttribute ) )...
                - (expParams.Streaming.dmaFramesAttribute/expParams.LSM.frameRateAttribute);
                %This extraordinarily laborious (not-yet-generalisability-tested) line calculates the imaging duration by dividing the number of frames by the (apparent) framerate...
                %       then adding the flyback overhead as a factor of the actual frames duration
                %           It is probably accurate within +- 1 minute (Although improvements based on empirical discoveries are feasible)
                %               %Added: Subtraction of dmaFrames, because reasons (It seemed to bring calculations closer to 'real' qualitative durations)
                %           Note: It is probably not affected by very late experiment termination/etc
                %}
                imagingEndTime =  imagingStartTime + seconds(expDurSecs);
                disp([strcat( 'Derived block duration: ',string( seconds(expDurSecs) ) )])
            elseif durationMethod == 2
            %New, dateModified-based system
                tifDeets = dir( [currentDataDirectory,filesep,'ChanA_Preview.tif'] );
                if isempty(tifDeets)
                    ['## TIF file could not be found! ##']
                    crash = yes
                end

                imagingEndTime = tifDeets.date;

            elseif durationMethod == 3
                 %Pi-based method
                 ['Not implemented yet']
                 crash = yes
            end
            
            disp([ strcat(string(imagingStartTime),'=>',string(imagingEndTime)) ])
            disp([ strcat('Calculated recording duration: ',string(imagingEndTime-imagingStartTime),'[Method ',num2str(durationMethod),']' )])

            %Place in behavioural coordinates
            behavStartIndex = find(time(between(imagingStartTime,allDatetimes)) >= 0, 1); %Index of behav data where 2p data starts
            %(If behav started after 2p this will probably be 0)
            behavEndIndex = find(time(between(imagingEndTime,allDatetimes)) >= 0, 1); %Index of behav data where 2p data (~)ends
            %QA
            if isempty( behavStartIndex ) %|| isempty( behavEndIndex )
                ['## Error in finding 2p data start in behavioural reference frame ##']
                crash = yes
            elseif isempty( behavEndIndex ) && seconds(time(between( allDatetimes(end),imagingEndTime ))) <= expDurSecs*0.05
                disp(['-# (Apparent) 2p imaging end time postcedes behav timescale by ',num2str(seconds(time(between( allDatetimes(end),imagingEndTime ))) ),'s (<5%); Using last behav timepoint #-'])
                behavEndIndex = size(allDatetimes,1);
            elseif isempty( behavEndIndex ) && seconds(time(between( allDatetimes(end),imagingEndTime ))) > expDurSecs*0.05
                ['-# (Apparent) 2p imaging end time postcedes behav timescale by ',num2str(seconds(time(between( allDatetimes(end),imagingEndTime ))) ),'s (>5%) #-']
                crash = yes
            end
            thisAcData = allData.avCntrSize( behavStartIndex:behavEndIndex ); %Acquire behav data of this block timeframe (Note: Behav timescale)
            
            %thisAcData = ...
            %   interp1([1:size(thisAcData,1)], thisAcData , linspace(1,size(thisAcData,1), size(thisData,ndims(thisData)) ) )';
            %Interpolate behavioural data to 2p timescale
            
            %----
            
            %Derive sequence data (Stolen from 2p processing functions)
            % get sequences
            %if exist(fullfile(sequenceDirectory,[flyID '.mat']),'file')
            %    randomSequence = load(fullfile(sequenceDirectory,[flyID '.mat']),'randomSequence').';
            %    BLOCKS(b).randomSequence = randomSequence.randomSequence-1;
            %else
            %BLOCKS(b).randomSequence = csvread(fullfile(sequenceDirectory,'Data_LEDs',[flyID '.csv'])).';
            randomSequence = csvread(fullfile(sequenceDirectory,[flyID '.csv'])).';
            %end
            
            nBadTrials = 0;
            nBadBlankTrials = 0;
            
            %remove bad trials and associated frames (this should be put inside a function)
            if ~isempty(currentBlock.removeFrames{1})
                %This segment untested in this script context
                removeFrames = eval(currentBlock.removeFrames{1});
                badTrials = eval(currentBlock.badTrials{1});
                nBadTrials = length(badTrials);
                
                %BLOCKS(b).greenChannel(:,:,removeFrames) = [];
                
                %auxRandomSequence = reshape(BLOCKS(b).randomSequence,[currentBlock.nStimuli length(BLOCKS(b).randomSequence)/currentBlock.nStimuli]);
                auxRandomSequence = reshape(randomSequence,[currentBlock.nStimuli length(randomSequence)/currentBlock.nStimuli]);
                
                nBadBlankTrials = sum(auxRandomSequence(1,badTrials) == 5);% 5's were used for a blank trial
                
                auxRandomSequence(:,badTrials) = [];
                
                %BLOCKS(b).randomSequence = auxRandomSequence(:).';
                randomSequence = auxRandomSequence(:).';
            
            end
            
            % calculate number of volumes per stimulus train
            %nVolTotal = size(BLOCKS(b).greenChannel,3);
            nVolTotal = currentBlock.realFrames; %No idea if correct transliteration
            %Should in fact be ~6000?
            
            %---
            savStruct = struct;
            savStruct.ancillary.flyID = flyID;
            savStruct.ancillary.dataDirectory = currentDataDirectory;
            savStruct.ancillary.behavIndsStartEnd = [behavStartIndex,behavEndIndex]; %Start and end (behav ref frame) of this block's data
            savStruct.ancillary.behavOriginalSize = size(thisAcData,1);
            savStruct.ancillary.behavFrameRate = behavFrameRate;
            savStruct.ancillary.derivedExpParams = expParams;
            savStruct.ancillary.imagingStartTime = imagingStartTime;
            savStruct.ancillary.durationMethod = durationMethod;
            savStruct.ancillary.imagingEndTime = imagingEndTime;

            savStruct.randomSequence = randomSequence; %Mostly unnecessary, but potentially handy to (re)save
            
            thisAcDataInterp = ...
            interp1([1:size(thisAcData,1)], thisAcData , linspace(1,size(thisAcData,1), length(randomSequence) ) )';
                            %Interpolate behavioural data to 2p timescale
                            %This variable represents raw activity data, on the same timescale as randomSequence is recorded
            savStruct.thisAcDataInterp = thisAcDataInterp;
            if doAcInac == 1
                %Interpolate binary
                thisInacBinary = inacBinaryProc( behavStartIndex:behavEndIndex ); %Subselect based on this experiment coords (as above)
                thisInacBinaryInterp = ...
                interp1([1:size(thisInacBinary,1)], double(thisInacBinary) , linspace(1,size(thisInacBinary,1), length(randomSequence) ), 'previous' )'; %Note conversion to double for interp
                %Note: Use of 'previous' (basically ZOH) method here probably critical to prevent fractional instead of binary values in interpolated data
                savStruct.acInac.thisInacBinaryInterp = thisInacBinaryInterp;
                savStruct.acInac.ancillary.acMeanThresh = acMeanThresh;
                savStruct.acInac.ancillary.minAcTime = minAcTime;
                savStruct.acInac.ancillary.minSleepBoutTime = minSleepBoutTime;

                %savStruct.acInac.inStruct = inStruct; %This saves the entire inStruct, which is less than useful
                %Subselect inStruct
                thisInStruct = struct; %Cutdown copy of the total inStruct, for just this recording block
                temp2 = [behavStartIndex:behavEndIndex];
                c = 1;
                for bout = 1:size(inStruct.holeStarts,2)
                    temp = [inStruct.holeStarts(bout):inStruct.holeEnds(bout)];
                    temp3 = intersect( temp2, temp ); %Allows for sleep bouts that cross recording boundaries/etc
                        %Note: Technically means you can have sleep bouts of like, 10 frames/etc, even if that is below the standard minimum requirements
                    if ~isempty(temp3)
                        thisInStruct.holeSizes(c) = length(temp3);
                        thisInStruct.holeStarts(c) = temp3(1);
                        thisInStruct.holeEnds(c) = temp3(end);

                        c = c + 1;
                    end
                end
                savStruct.acInac.thisInStruct = thisInStruct; %Cutdown version
                disp(['(',num2str(c - 1),' sleep bouts occurred in this recording block specifically)'])
                    %Note: May add to more than total number of sleep bouts if crossing boundaries/etc

            end
            
            %And save
            if saveOutput == 1
                save([currentDataDirectory,filesep,'behavSequence'],"savStruct"); %Save to 2p data directory
                disp(['Saved behavioural data (synced to sequence) to ',char(10),currentDataDirectory])
            else
                disp(['(Behavioural data not saved by request)'])
            end

            overSave(b).savStruct = savStruct;
            
        end
        
        
        %------------
    end
    
    %And save
    if saveOutput == 1
        save( [thisVidFolder,filesep,'fly',num2str(flyOnDay),'_mov'], 'allData' ,'allDatetimes','thisCSVs') %Save to vid directory
    else
        disp(['(Not saved)'])
    end
    disp(['Processed in ',num2str(toc),'s'])
    
end %chosenFlies end

disp([char(10),'Processing of ',num2str(length(chosenFlies)),' datasets complete'])
if failCount > 0
    ['-# There were ',num2str(failCount),' failures #-']
end
