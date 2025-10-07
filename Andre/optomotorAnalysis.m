% Optomotor Analysis

clear
close all

%reading fly record
blocks = readtable("I:\RFDG2021-Q4413\Andre\2p_Record\Andre_2P_record.xlsx");
blocks = blocks(~logical(blocks.Exclude),:);
noDataList = [];
nonBatteryList = [];

%directories
dataDirectory = 'I:\RFDG2021-Q4413\Andre\2p_Data';

% chosenFlies = [85, 86, 87, 89, 91, 104, 106, 108:110, 138:140];
% chosenBlocks = repmat({1}, 1, size(chosenFlies, 2));

% chosenFlies = 96:103;
% chosenBlocks = {[1], [1], [1], [1], [1], [2], [1], [1]};

%back-of-heads
chosenFlies = [126:128, 141:145]; %56 pre-dates batSeq with all conditions

chosenBlocks = {[1], [1], [1], [1], [1], [1], [1], [1]}; %come back to later, maybe not, could do every battery case using nStimuli == -2

%options
figPlots = [0, 0, 1, 0, 1, 1]; %[veloc, dis, dXs, ficTrack, optoLefts, optoRights]
filterInactive = 0; %option to remove directional displacements based on linearity (based on assumption that fly doesn't move ball in a linear fashion and so must be caused by gas-induced ball rotation)
shadedError = 1;

%Quick Plot options
immediateAnalysis = 0; % whether to use flyRecord or date and info passed below
quickFlyDirectory = 'I:\RFDG2021-Q4413\Andre\2p_Data\5Jun25\MAT'; %must target a single folder
quickFlyInfo = '050625_f2_b1'; %fly details as expected for file name (ddmmyy_f#_b#)

%%

%predefining important variables and tests
xtest = [];
x2test = [];

%QA
if immediateAnalysis == 1
    disp(['\/\/\/ WARNING: quick fly analysis chosen, fly record info disregarded \/\/\/'])
end

if immediateAnalysis == 1
    chosenFlies = [-1];
    nBlocks = 1;
end

%% collection and individual plots
optoStruct = struct;
structCount = 1;
for fly = 1:length(chosenFlies)
    if immediateAnalysis ~= 1
        thisFly = chosenFlies(fly);
        thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);

        if ~isempty(chosenBlocks) && ~isempty(chosenBlocks{fly})
            thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
        end

        nBlocks = height(thisFlyBlocks);

        if nBlocks == 0
            ['-# Alert: No apparent blocks (for fly ',num2str(chosenFlies(fly)),') #-']
            %crash = yes
            noDataList = [noDataList,thisFly];
            continue
        end
    end

    for b = 1:nBlocks
        if immediateAnalysis ~= 1
            currentBlock = thisFlyBlocks(b,:);
            currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
            %         flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
            currentDate_2 = char(datetime(thisFlyBlocks.Date(1),'Format','ddMMyy')); %2nd date uses different format
            flyID = [currentDate_2, '_f', num2str(currentBlock.FlyOnDay), '_b', num2str(currentBlock.Block)];
            %         guacamole

            if currentBlock.nStimuli == -2
                btDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btData.csv'];
                btAncillaryDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btAncillary.mat'];
                matDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_MAT.mat'];
            else
%             ['-| Alert: non-battery data detected for ' (num2str(chosenFlies(fly))), '|-']
            nonBatteryList = [nonBatteryList, flyID];
            continue
            end
        else % quick analysis (immediateAnalysis == 1) - intended only for use 
            btDataDirectory = [quickFlyDirectory, filesep, quickFlyInfo, '_btData.csv'];
            btAncillaryDirectory = [quickFlyDirectory, filesep, quickFlyInfo, '_btAncillary.mat'];
            matDataDirectory = [quickFlyDirectory, filesep, quickFlyInfo, '_MAT.mat'];
        end
            %%%%%%%%%% actual stuff happens %%%%%%%%%%%%%%

            %Loading files
            BTData = csvread(btDataDirectory);
            load(btAncillaryDirectory);
            load(matDataDirectory);

            %defining variables
            numOpto = matSave.numOpto;
            optoElements = matSave.optoElements;
            panelFrequency = matSave.panelFrequency;
            optoDur = matSave.optoCycleDuration;
            
            if isfield(matSave, 'batSequence')
                batSequence = matSave.batSequence;
            else
                %rebuilding batSequence 
                batConditions = [optoElements, matSave.numLR*2, matSave.numLum*2, matSave.freqElements];
                batSequence = [];
                for i = [1,2,3,4]
                    batSequence = [batSequence, i*ones(1, batConditions(i))];
                end
                batSequence = [batSequence, -1]; %has an added element at end to match sequenceLength
            end

            optoStartInd = find(batSequence == 1, 1, 'first');
            optoEndInd = find(batSequence == 1, 1, 'last');           
           
%             flamenate

            %ensuring phaseTrack clean
            if size(btAncillary.phaseTrack, 2) > optoElements
                btAncillary.phaseTrack(optoElements + 1:size(btAncillary.phaseTrack, 2)) = [];
            end
%             checkengine

            %setting up phaseTrack
            phaseTrackUnwrap = unwrap(deg2rad(btAncillary.phaseTrack)); 
            phaseTrackDiff = diff(phaseTrackUnwrap);
%             for i = 1:length(phaseTrackDiff)
%                 if phaseTrackDiff(i) == 0
%                     phaseTrackDiff(i) = mean([phaseTrackDiff(i-1), phaseTrackDiff(i+1)]);
%                 end
%             end
            %cleaning up phaseTrackDiff
            %0.3 currently arbitrary - comes from rate of phase change (using 0.3 could get buggy if temporal frequency < 3Hz)
            phaseThreshold = max(phaseTrackDiff) - 0.01; %0.01 arbitrary but small (should only get funky if temporal frequency ridiculously slow)
            phaseTrackDiff(phaseTrackDiff > phaseThreshold) = 1;
            phaseTrackDiff(phaseTrackDiff < -phaseThreshold) = -1;
%             phaseTrackDiff(phaseTrackDiff > -phaseThreshold & phaseTrackDiff < phaseThreshold) = 0;
            phaseTrackDiff(phaseTrackDiff > -1 & phaseTrackDiff < 1) = 0;

            phaseTrack = btAncillary.phaseTrack;

            for i = 1:length(phaseTrackDiff)
                if phaseTrackDiff(i) == 0
                    if i == 1
                        phaseTrackDiff(i) = phaseTrackDiff(i + 1);
                    elseif i == length(phaseTrackDiff)
                        phaseTrackDiff(i) = phaseTrackDiff(i - 1);
%                     elseif phaseTrackDiff(i+1) ~= 0 && phaseTrackDiff(i-1) ~= 0 
%                         phaseTrackDiff(i) = nanmean([phaseTrackDiff(i-1), phaseTrackDiff(i+1)]);
                    elseif phaseTrack(i) ~= phaseTrack(i-1) && phaseTrack(i) ~= phaseTrack(i+1)
                        if phaseTrackDiff(i+1) == 0
                            phaseTrackDiff(i) = phaseTrackDiff(i-1);
                        elseif phaseTrackDiff(i-1) == 0
                            phaseTrackDiff(i) = phaseTrackDiff(i+1);
                        else
                            phaseTrackDiff(i) = nanmean([phaseTrackDiff(i-1), phaseTrackDiff(i+1)]);
                        end 
                    end
                end
            end   

            %cheating to get phase to work for fast freqs
%             for i = 1:size(phaseTrackDiff, 2)
%                 if phaseTrackDiff(i) == 0
%                     phaseTrackDiff(i) = phaseTrackDiff(i+1);
%                 end
%             end

%             phaseTrackDiff (phaseTrackDiff < 1) = 0;
%             phase = bwlabel(phaseTrackDiff);

            % if no stationary periods, setting one direction (-1) to 0 so that bwlabel works
            if phaseTrackDiff(:) ~= 0
                modPhaseTrackDiff = phaseTrackDiff;
                modPhaseTrackDiff(modPhaseTrackDiff < 1) = 0;
                phase = bwlabel(modPhaseTrackDiff);
            else
                phase = bwlabel(phaseTrackDiff);
            end
%             kaplunk

%             uniCoords = []; %for filling in alternating phases later
%             for i = 1:nanmax(phase)
%                 uniCoords(i, 1) = find(phase == i, 1, 'first');
%                 uniCoords(i, 2) = find(phase == i, 1, 'last');
%             end
% 
%             otherCoords = [];
%             for i = 1:nanmax(phase)
%                 if i == nanmax(phase)
%                     otherCoords(i, 1) = uniCoords(i, 2) + 1;
%                     otherCoords(i, 2) = size(btAncillary.phaseTrack, 2);
%                 else
%                     otherCoords(i, 1) = uniCoords(i, 2) + 1;
%                     otherCoords(i, 2) = uniCoords(i + 1, 1) - 1;
%                 end
%             end
% 
%             phaseCoords = [];
%             for i = 1:nanmax(phase)
%                 phaseCoords = [phaseCoords; uniCoords(i,[1:2]); otherCoords(i,[1:2])];
%             end

%             phaseCoords = [];
%             for i = 1:nanmax(phase)
%                 phaseCoords(i, 1) = find(phase == i, 1, 'first');
%                 phaseCoords(i, 2) = find(phase == i, 1, 'last');
%             end
% 
%             phaseCoords = [];
%             for i = 1:2*nanmax(phase)
%                 if i - 1 == 0 || phaseCoords(i-1, 1) == 0
%                     phaseCoords(i, 1) = uniPhaseCoords(i - 1, 1);
%                     phaseCoords(i, 2) = uniPhaseCoords(i - 1, 2);
%                 else
%                     phaseCoords(i, 1) = 0;
%                     phaseCoords(i, 2) = 0;
%                 end
%             end


%             phaseCoords = []; %coords for start and end of each optomotor phase (left/right/stationary)
%             for i = 1:2:2*nanmax(uniPhaseCoords, 1) 
%                 phaseCoords(i, 1) = uniPhaseCoords(i, 1);
%                 phaseCoords(i + 1, 1)= uniPhaseCoords(i, 1) + 1;
%                 phaseCoords(i, 2) = uniPhaseCoords(i, 2);
%                 phaseCoords(i + 1, 2)= uniPhaseCoords(i, 2) + 1;
%             end
%             bandicoot

%             phaseRight = [];
%             for i = 1:size(phase, 2)
%                 if phase(i) == 0
%                     phaseRight = [phaseRight, 1];
%                 else
%                     phaseRight = [phaseRight, 0];
%                 end
%             end
%             phaseRight = bwlabel(phaseRight);
%             sputnik

            %QA
            if size(btAncillary.btAncillaryArray, 1) < optoElements
                disp(['-_- WARNING: optomotor did not complete -_-'])
                disp(['Using maximum number of elements instead'])
                optoElements = size(btAncillary.btAncillaryArray, 1);
                numOpto = max(phase) - 1;
            end

            % if blank before opto included, skipping to next i due to frequency error in stimlus code
            if isfield(matSave, 'blank')
                disp(['-_- WARNING: first opto direction is 1/60th shorter -_-'])
                disp(['-_- note: currently assumes frequency error in all instances of code containing blanks'])
%                 phaseCoords(1, 1) = phaseCoords(1, 1) + 1;
                optoStartInd = optoStartInd + 1;
                optoEndInd = optoEndInd + 1;
            end


            leftCount = 1;
            rightCount = 1;
            % phase splicing
            leftIts = [nan*ones(ceil(size(phase, 2)/numOpto), 1)];
            rightIts = [nan*ones(ceil(size(phase, 2)/numOpto), 1)];
            lefts = [];
            rights = [];
            for i = 1:size(phase, 2)

                if phase(i)
                    lefts = [lefts; i];
                else
                    rights = [rights; i];
                end

                if phase(i) == 0 && ( i == 1 || phase(i-1) ~= 0 )
                    for j = 1:size(lefts, 1)
                        leftIts(j,leftCount) = lefts(j);
%                       rightIts{:, directionCount} = rights;
                    end
                    lefts = [];
                    leftCount = leftCount + 1;

%                 elseif (i == size(phase, 2) || phase(i)) && (i == 1 || phase(i-1) == 0)
                elseif phase(i) == 0 && (i ==size(phase, 2) || phase(i+1))
                    for f = 1:size(rights, 1)
                        rightIts(f, rightCount) = rights(f);
                    end
                    rights = [];
                    rightCount = rightCount + 1;
                end
% 
%                 if phase(i)
%                     lefts = [lefts; i];
%                 else
%                     rights = [rights; i];
%                 end
            end
%             leftIts = leftIts(:, 2:size(leftIts, 2)); %remove initialisation row of nan <--- FALSE
            leftIts(leftIts(:,:) == 0) = nan; %make all 0's nans (i should never be 0), ensures averaging works
            rightIts(rightIts(:,:) == 0) = nan;
%             crush

            %ensures correct i values stored if optomotor isn't first condition
            leftIts = leftIts + (optoStartInd - 1); 
            rightIts = rightIts + (optoStartInd - 1); 

            %generating phaseCoords using leftIts
            phaseCoords = [];
            %Assumes optomotor always starts with left direction and that always same number each direction
            for i = 1:size(leftIts, 2)
                leftCoords = [leftIts(1, i), leftIts(find(isnan(leftIts(:, i)), 1, 'first') - 1, i)];
                rightCoords = [rightIts(1, i), rightIts(find(isnan(rightIts(:, i)), 1, 'first') - 1, i)];
                combinedCoords = [leftCoords; rightCoords];
                phaseCoords = [phaseCoords; combinedCoords];
            end
%             splash
            
            % if blank before opto included, skipping to next i due to frequency error in stimlus code
%             if isfield(matSave, 'blank')
%                 disp(['-_- WARNING: first opto direction is 1/60th shorter -_-'])
%                 disp(['-_- note: currently assumes frequency error in all instances of code containing blanks'])
%                 phaseCoords(1, 1) = phaseCoords(1, 1) + 1;
%                 optoStartInd = optoStartInd + 1;
%             end


%             %QA
%             if size(btAncillary.btAncillaryArray, 1) < optoElements
%                 disp(['-_- WARNING: optomotor did not complete -_-'])
%                 disp(['Using maximum number of elements instead'])
%                 optoElements = size(btAncillary.btAncillaryArray, 1);
%             end

            % grabbing BT velocites
            rawBT = cell2mat(btAncillary.btAncillaryArray(optoStartInd:optoEndInd,2));
            vX1 = signedHex2dec( rawBT(:, 6:7));
            vX2 = signedHex2dec( rawBT(:, 12:13));
            vY1 = signedHex2dec( rawBT(:, 8:9));
            vY2 = signedHex2dec( rawBT(:, 14:15));

            %%generating integrals
            % distance/displacement per axis
            dX1 = 0;
            dX2 = 0;
            dY1 = 0;
            dY2 = 0;

            for i = 1:optoElements-1
                if i == 1
                    dX1(i) = 0 + vX1(i);
                    dX2(i) = 0 + vX2(i);
                    dY1(i) = 0 + vY1(i);
                    dY2(i) = 0 + vY2(i);
                else
                    dX1(i) = dX1(i-1) + vX1(i);
                    dX2(i) = dX2(i-1) + vX2(i);
                    dY1(i) = dY1(i-1) + vY1(i);
                    dY2(i) = dY2(i-1) + vY2(i);
                end
            end

            %dX1(1:540) = sin(dX1(1:540))*1e4;

            %mapping displacements to directions
%             leftRotate =[];
%             rightRotate = [];
%             for i = 1:numOpto
%                 for b = 1:size(leftIts, 1)
%                     if ~isnan(leftIts(b, i))
%                         leftRotate(b, i) = dXleftIts(b, i);
%                     end
%                     if ~isnan(rightIts(b, i))
%                         rightRotate(b, i)= rightIts(b, i);
%                     end
%                 end
%             end
%             leftRotate(leftRotate(:, :)==0) = nan;
%             rightRotate(rightRotate(:,:)==0) = nan;


            dX1PerDirection = [nan*ones(numOpto, size(leftIts, 1), 2)]; % [ rep, number of i elements, direction of grating ] 
            dX2PerDirection = [nan*ones(numOpto, size(leftIts, 1), 2)]; % [ ^^ ]

            for i = 1:numOpto
%                 leftRotate = [];
%                 rightRotate = [];
%                 for j = 1:size(leftIts, 1)
%                     leftRotate = dX
%                 end
                for j = 1:size(leftIts, 1)
                     if ~isnan(leftIts(j, i))
                        dX1PerDirection(i, j, 1) = dX1(leftIts(j, i) - (optoStartInd - 1));
                        dX2PerDirection(i, j, 1) = dX2(leftIts(j, i) - (optoStartInd - 1));
                    end

                    if ~isnan(rightIts(j, i))
                        dX1PerDirection(i, j, 2) = dX1(rightIts(j, i) - (optoStartInd - 1));
                        dX2PerDirection(i, j, 2) = dX2(rightIts(j, i) - (optoStartInd - 1));
                    end

                    %dX1PerDirection(i, :, 1) = leftIts(:, i)';
                    %dX1PerDirection(i, :, 2) = rightIts(:, i)';
%                     if ~isnan(leftIts(j, i))
%                         dX1PerDirection(i, j, 1) = dX1(leftIts(j, i));
%                         dX2PerDirection(i, j, 1) = dX2(leftIts(j, i));
%                     end
% 
%                     if ~isnan(rightIts(j, i))
%                         dX1PerDirection(i, j, 2) = dX1(rightIts(j, i));
%                         dX2PerDirection(i, j, 2) = dX2(rightIts(j, i));
%                     end

%                 for j = 1:size(rotatePerDirection, 2)
%                     if ~isnan(rotatePerDirection(i, j))
%                         
%                     end
                end
            end

%             for i = size(dX1PerDirection, 2)
%                 for b = 1:size(dX1PerDirection, 1)
%                     if ~isnan(dX1PerDirection(b, i, 1))
%                         dX1PerDirection(b, i, 1) = dX1(dX1PerDirection(b, i, 1));
%                     end
%                 end
%             end

%             smoothie

%             for i = 1:size(rotatePerDirection)
%                 if ~isnan(rotatePerDirection)
%                 end
%             end

            %generating averages for rotation in each direction



            %distance/displacement in single axis (y-axis?)
            fV = []; % forward velocity (math taken from xMMPerSec)
            hV = []; %  horizontal velocity (math taken from yMMPerSec)
            aMod = sin(pi/4) *( (6*pi)/2750 * (200/4) ); % formula taken from FSD script (maybe make 200 and 4 triggerRate and divider respectively eventually)

            %maybe add in a boolean for different sensorPositions eventually
            for i = 1:optoElements
                fV(i) = ( vY2(i) - vY1(i) ) * aMod;
                hV(i) = -( vY1(i) + vY2(i) ) * aMod;
            end
            % guava

            fD = 0;
            hD = 0;
            for i = 2:optoElements
                fD(i) = fD(i-1) + fV(i);
                hD(i) = hD(i-1) + hV(i);
            end


            %%plotting stuff
            optoFinish = BTData(find(BTData(:, 5)==optoEndInd, 1), 6);
            if optoStartInd == 1
                optoStart = BTData( find( (BTData(:, 5) == 2), 1)-1, 6 );
            else
                optoStart = BTData( find( (BTData(:, 5) == optoStartInd), 1), 6);
            end
%             optoTime = BTData(find(BTData(:, 5)==optoElements, 1), 6) - BTData( find( (BTData(:, 5) == 2), 1)-1, 6 ); %boils down to (time at end of opto) - (time at start of opto)
            optoTime = optoFinish - optoStart;
            % shanghai

            %velocities
            if figPlots(1)
                figure
                plot(linspace(0, optoTime, optoElements), vX1)
                hold on
                for velo = [vX2, vY1, vY2]
                    plot(linspace(0, optoTime, optoElements), velo)
                end
                %             title("(vX1, vX2, vY1, vY2) vs Time (s)") %change later to be better descritption
                if immediateAnalysis ~= 1
                    velocityTitle = ['f', num2str(thisFly), ' b', num2str(currentBlock.Block), ': velocity vs time(s)'];
                else
                    velocityTitle = [quickFlyInfo, ': velocity vs time(s)'];
                end
                title(velocityTitle)
                legend([{'vX1'},{'vX2'},{'vY1'},{'vY2'},])
            end

            %distance/displacement
            % tic
            if figPlots(2)
                figure
                plot(linspace(0, optoTime, optoElements), dX1)
                hold on
                for dist = [dX2', dY1', dY2']
                    plot(linspace(0, optoTime, optoElements), dist)
                    %     toc
                end
                %             title("(dX1, dX2, dY1, dY2) vs Time (s)")
                if immediateAnalysis ~= 1
                    displacementTitle = ['f', num2str(thisFly), ' b', num2str(currentBlock.Block), ': displacement vs time(s)'];
                else
                    displacementTitle = [quickFlyInfo, ': displacement vs time(s)'];
                end
                title(displacementTitle)
                legend([{'dX1'},{'dX2'},{'dY1'},{'dY2'},])
            end

            %with phase of optomotor
            if figPlots(3)
                figure
                plot(linspace(0, optoTime, size(dX1, 2)), dX1, 'LineWidth', 2, 'Color', 'k')
                hold on
                plot(linspace(0, optoTime, size(dX2, 2)), dX2, 'LineWidth', 2, 'Color', 'k')

                if immediateAnalysis ~= 1
                    rotateTitle = ['f', num2str(thisFly), ' b', num2str(currentBlock.Block), ': rotation vs time(s)'];
                else
                    rotateTitle = [quickFlyInfo, ': rotation vs time(s)'];
                end
                title(rotateTitle)

                %             fill(optoDur, max([max(dX1), max(dX2)]), [0 0 1])
                maxDX = max([max(dX1), max(dX2)]);
                minDX = min([min(dX1), min(dX2)]);

                % ##--------------Theoretical calc --------------##
                %             altBright = 1;
                %             for x = [optoDur: optoDur: optoTime]
                %                 initX = x - optoDur;
                %                 if altBright == 1
                %                     fill( [initX, initX, x, x], [-20, maxDX, maxDX, -20], 'b', 'FaceAlpha', 0.25) %-20 arbitrary, just want fill to be larger than figure limits
                %                     altBright = 0;
                %                 else
                % %                     fill( [initX, initX, x, x], [0, maxDX, 0, maxDX], 'b')
                %                     altBright = 1;
                %                 end
                %             end
                %             axis([0, optoTime, 0, maxDX])

                % ##---------using BT for time--------##
                for i = 1:size(phaseCoords, 1)
                    %                 initX = x - optoDur;
                    phaseStart = phaseCoords(i, 1);
                    phaseEnd = phaseCoords(i, 2);
                    startTime = optoStart;
                    %                 phaseIteration = phaseC

                    x = BTData(find(BTData(:, 5)== phaseEnd, 1), 6) - startTime;

                    if phaseStart == 1
                        initX = x - (x - BTData(find(BTData(:, 5)== 2, 1)-1 , 6)) - startTime; %x minus the difference between currentTime(x) and time at previous i
                    else
                        initX = x - (x - BTData(find(BTData(:, 5)== phaseStart, 1), 6)) - startTime;
                    end

                    if phaseTrackDiff(phaseStart) > 0
                        fill( [initX, initX, x, x], [minDX, maxDX, maxDX, minDX], 'b', 'FaceAlpha', 0.15)
                    elseif phaseTrackDiff(phaseStart) < 0
                        fill( [initX, initX, x, x], [minDX, maxDX, maxDX, minDX], 'r', 'FaceAlpha', 0.15)
                        %                     fragile
                    end
                end
                axis([0, optoTime, minDX, maxDX])
                set(gcf, 'Name', ['fly ', num2str(chosenFlies(fly)), ' block ', num2str(chosenBlocks{fly}), ' rotational displacement'])
            end
            % %saving
            % optoStruct(structCount).phaseCoords = phaseCoords;
            % optoStruct(structCount).dX1Times = 
            %^^ DISCONTINUED - will take too long to implement cross fly for full opto

            %%fictive Track
            %using single axis
            if figPlots(4)
                figure
                plot(hD, fD)
                hold on
                for i = 1:floor(length(hD)/100):length(hD)
                    scatter( hD(i), fD(i) )
                end

                if immediateAnalysis ~= 1
                    ficTrackTitle = ['f', num2str(thisFly), ' b', num2str(currentBlock.Block), ': ficTrack'];
                else
                    ficTrackTitle = [quickFlyInfo, ': ficTrack'];
                end
                title(ficTrackTitle)
            end

            %rotation lefts
            if figPlots(5)
                firstNanX1 = find(isnan(mean( dX1PerDirection(:,:,1), 1 )), 1, 'first');
                firstNanX2 = find(isnan(mean( dX2PerDirection(:,:,1), 1 )), 1, 'first');
                figure
                hold on
                %             leftLegend = [];
                for i = 1:numOpto
                    %                 t1 = BTData(BTData(:, 5) == dX1PerDirection(i, size(dX1PerDirection(i), 2), 1), 6);
                    %                     t1 = BTData(BTData(:, 5) == leftITs(find(isnan(leftIts(:, i), 1) - 1 i) )
                    iVal = leftIts(find(isnan(leftIts(:, i)), 1) - 1, i);
                    t1 = BTData(find(BTData(:, 5) == iVal, 1, 'last'), 6);
                    if leftIts(1, i) == 1
                        cycleTime = t1 - optoStart;
                    else
                        cycleTime = t1 - BTData(find(BTData(:, 5) == leftIts(1, i), 1), 6);
                    end
                    x2test = [x2test, cycleTime];

                    %                     firstNanX1 = find(isnan(mean( dX1PerDirection(:,:,1), 1 )), 1, 'first');
                    %                     firstNanX2 = find(isnan(mean( dX2PerDirection(:,:,1), 1 )), 1, 'first');

                    %                     plot(linspace(0, cycleTime, find(isnan(dX1PerDirection(i, :, 1)), 1) - 1), dX1PerDirection(i, :, 1) - dX1PerDirection(i, 1, 1))
                    plot(linspace(0, cycleTime, firstNanX1), dX1PerDirection(i, [1:firstNanX1], 1) - dX1PerDirection(i, 1, 1))
                    plot(linspace(0, cycleTime, firstNanX2), dX2PerDirection(i, [1:firstNanX2], 1) - dX2PerDirection(i, 1, 1))
                    %                     xlim([0,cycleTime])
                    %                     plot(dX2PerDirection(i, :, 1) - dX2PerDirection(i, 1, 1))
                end
                %plotting averages
                plot(linspace(0, cycleTime, firstNanX1), nanmean(dX1PerDirection(:, [1:firstNanX1], 1) - dX1PerDirection(:, 1, 1)), 'Color', [0 0 0], 'LineWidth', 3)
                plot(linspace(0, cycleTime, firstNanX2), nanmean(dX2PerDirection(:, [1:firstNanX2], 1) - dX2PerDirection(:, 1, 1)), 'Color', [0 0 0], 'LineWidth', 3)

                %calculating SEM
                SEMdX1Left = nanstd(dX1PerDirection(:, [1:firstNanX1], 1) - dX1PerDirection(:, 1, 1), 1)/sqrt(size(dX1PerDirection, 1));
                SEMdX2Left = nanstd(dX2PerDirection(:, [1:firstNanX2], 1) - dX2PerDirection(:, 1, 1), 1)/sqrt(size(dX2PerDirection, 1));
                %adding error bars
                errorbar(linspace(0, cycleTime, firstNanX1), nanmean(dX1PerDirection(:, [1:firstNanX1], 1) - dX1PerDirection(:, 1, 1)), SEMdX1Left)
                errorbar(linspace(0, cycleTime, firstNanX2), nanmean(dX2PerDirection(:, [1:firstNanX2], 1) - dX2PerDirection(:, 1, 1)), SEMdX2Left)
                title("optomotor lefts")
                set(gcf, 'Name', ['fly ', num2str(chosenFlies(fly)), ' block ', num2str(chosenBlocks{fly}), ' optoLefts'])
            end

            % rotation rights
            if figPlots(6)
                firstNanX1 = find(isnan(mean( dX1PerDirection(:,:,1), 1 )), 1, 'first');
                firstNanX2 = find(isnan(mean( dX2PerDirection(:,:,1), 1 )), 1, 'first');

                figure
                hold on
                for i = 1:numOpto
                    iVal = rightIts(find(isnan(rightIts(:, i)), 1) - 1, i);
                    t1 = BTData(find(BTData(:, 5) == iVal, 1, 'last'), 6);
                    if rightIts(1, i) == 1
                        cycleTime = t1 - optoStart;
                    else
                        cycleTime = t1 - BTData(find(BTData(:, 5) == rightIts(1, i), 1), 6);
                    end

                    plot(linspace(0, cycleTime, firstNanX1), dX1PerDirection(i, [1:firstNanX1], 2) - dX1PerDirection(i, 1, 2))
                    plot(linspace(0, cycleTime, firstNanX2), dX2PerDirection(i, [1:firstNanX2], 2) - dX2PerDirection(i, 1, 2))
                end
                %plotting averages
                plot(linspace(0, cycleTime, firstNanX1), nanmean(dX1PerDirection(:, [1:firstNanX1], 2) - dX1PerDirection(:, 1, 2), 1), 'Color', [0 0 0], 'LineWidth', 3)
                plot(linspace(0, cycleTime, firstNanX2), nanmean(dX2PerDirection(:, [1:firstNanX2], 2) - dX2PerDirection(:, 1, 2), 1), 'Color', [0 0 0], 'LineWidth', 3)

                %calculating SEM
                SEMdX1Right = nanstd(dX1PerDirection(:, [1:firstNanX1], 2) - dX1PerDirection(:, 1, 2), 1)/sqrt(size(dX1PerDirection, 1));
                SEMdX2Right = nanstd(dX2PerDirection(:, [1:firstNanX2], 2) - dX2PerDirection(:, 1, 2), 1)/sqrt(size(dX2PerDirection, 1));
                %adding error bars
                errorbar(linspace(0, cycleTime, firstNanX1), nanmean(dX1PerDirection(:, [1:firstNanX1], 2) - dX1PerDirection(:, 1, 2)), SEMdX1Right)
                errorbar(linspace(0, cycleTime, firstNanX2), nanmean(dX2PerDirection(:, [1:firstNanX2], 2) - dX2PerDirection(:, 1, 2)), SEMdX2Right)

                title("optomotor rights")
                set(gcf, 'Name', ['fly ', num2str(chosenFlies(fly)), ' block ', num2str(chosenBlocks{fly}), ' optoRights'])
            end

%         else
% %             ['-| Alert: non-battery data detected for ' (num2str(chosenFlies(fly))), '|-']
%             nonBatteryList = [nonBatteryList, flyID];
%             continue
%         end

        %saving 
        optoStruct(structCount).fly = thisFly;
        optoStruct(structCount).block = table2array(thisFlyBlocks(b, "Block"));
        optoStruct(structCount).cycleTime = cycleTime;
        optoStruct(structCount).dX1Direction = dX1PerDirection;
        optoStruct(structCount).dX2Direction = dX2PerDirection;
        structCount = structCount + 1;
    end

%     BTData = csvread("I:\RFDG2021-Q4413\Andre\2p_Data\22Apr25\MAT\220425_f1_b2_btData.csv");
%     load("I:\RFDG2021-Q4413\Andre\2p_Data\22Apr25\MAT\220425_f1_b2_btAncillary.mat");
%     load("I:\RFDG2021-Q4413\Andre\2p_Data\22Apr25\MAT\220425_f1_b2_MAT.mat");

end

%% cross-fly analyses
if ~immediateAnalysis
    blockNums = [optoStruct.block];
    crossFlydX1 = nan(size(chosenFlies, 2), 720, max(blockNums), 2); %720 arbitrary (maybe optoduration*panelfrquency?) (is there a way to find largest matrix in a field of a struct?)
    crossFlydX2 = nan(size(chosenFlies, 2), 720, max(blockNums), 2); % [ fly, displacement values, block, direction ] ^^same for above^^
    currentFly = optoStruct(1).fly;
    flyCount = 1;
    for block = 1:size(blockNums, 2)
        if currentFly ~= optoStruct(block).fly
            currentFly = optoStruct.fly;
            flyCount = flyCount + 1;
        end

        currentdX1 = optoStruct(block).dX1Direction - optoStruct(block).dX1Direction(:,1,:);
        currentdX2 = optoStruct(block).dX2Direction - optoStruct(block).dX2Direction(:,1,:);

        %potentially filtering here
        if filterInactive
            for direction = [1,2]
                linies = nan(2, size(currentdX1, 1)); % [ BT1/2, length of cycle (in i) ]
                for o = 1:size(linies, 2)
                    %fitting ball movement to linear model
                    lin1 = fitlm(linspace(0, optoStruct(block).cycleTime, size(currentdX1, 2)), currentdX1(o, :, direction));
                    lin2 = fitlm(linspace(0, optoStruct(block).cycleTime, size(currentdX2, 2)), currentdX2(o, :, direction));
                    %storing quality of linear model
                    linies(1,o) = lin1.Rsquared.Ordinary;
                    linies(2,o) = lin2.Rsquared.Ordinary;
                end
                %removing reps which are linear (r-squared threhsold set to 0.9)
                currentdX1(find(linies(1,:)<0.9), :, direction) = nan;
                currentdX2(find(linies(2,:)<0.9), :, direction) = nan;
            end
        end

        crossFlydX1(flyCount, :, blockNums(block), :) = nanmean(currentdX1, 1);
        crossFlydX2(flyCount, :, blockNums(block), :) = nanmean(currentdX2, 1);
    end

    % transforming opto movements into statistically analysable data
    directionSums = nan(size(chosenFlies, 2), 2, 2, max(blockNums)); % [ fly, sum statistic for left(1) and right(1), BT# , block]
    for fly = 1:size(directionSums, 1)
        for block = 1:size(directionSums, 4)
            for direction = [1,2]
                if isnan(crossFlydX1(fly, 1, block, direction)) % if first element is nan - block doesn't exist - currently assumes correct regardless of which BT used
                    continue
                end

                directionSums(fly, direction, 1, block) = nansum(crossFlydX1(fly, :, block, direction));
                directionSums(fly, direction, 2, block) = nansum(crossFlydX2(fly, :, block, direction));
            end
        end
    end
    % ^^NOT SURE ABOUT FUNCTIONALITY/APPROPIATENESS

    directions = {'left', 'right'};
    avgCycleTime = nanmean([optoStruct.cycleTime]);
    %plotting cross-fly
    for block = 1:max(blockNums)
        for direction = [1,2]
            % avgCycleTime = nanmean([optoStruct.cycleTime]);

            %BT1
            figure
            hold on
            currentSize = find(isnan(crossFlydX1( 1,:,block, direction ) ), 1, 'first'); %may be out slightly
            % plot(linspace(0, avgCycleTime, find(isnan(crossFlydX1( 1,:,block, 1 ) ), 1, 'first')), crossFlydX1(:,1:361,block,1)')
            plot(linspace(0, avgCycleTime, currentSize), crossFlydX1(1:7,1:currentSize,block,direction)')
            plot(linspace(0, avgCycleTime, currentSize), nanmean(crossFlydX1(1:7,1:currentSize,block,direction)), 'LineWidth', 3, 'Color', 'k')

            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'Box', 'off')
            ylabel('rotational displacement (mm)')
            xlabel('time(s)')
            title(['dX1 ' , directions{direction}, 's, cross-fly block ' num2str(block)])

            %BT2
            figure
            hold on
            currentSize = find(isnan(crossFlydX2( 1,:,block, direction ) ), 1, 'first'); %may be out slightly
            % plot(linspace(0, avgCycleTime, find(isnan(crossFlydX1( 1,:,block, 1 ) ), 1, 'first')), crossFlydX1(:,1:361,block,1)')
            plot(linspace(0, avgCycleTime, currentSize), crossFlydX2(:,1:currentSize,block,direction)')
            plot(linspace(0, avgCycleTime, currentSize), nanmean(crossFlydX2(:,1:currentSize,block,direction)), 'LineWidth', 3, 'Color', 'k')

            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'Box', 'off')
            ylabel('rotational displacement (mm)')
            xlabel('time(s)')
            title(['dX2 ' directions{direction}, 's, cross-fly block ' num2str(block)])
        end
    end

    directionColour = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};
    dX1SEM = nanstd(crossFlydX1, [], 1)/sqrt(size(crossFlydX1, 1));
    %plotting dX1 average displacements
    for block  = 1:size(crossFlydX1, 3)
        figure 
        hold on
        for direction = [1,2] %1:size(crossFlydX1, 4)
            currentSize = find(~isnan(crossFlydX1( 1,:,block, direction ) ), 1, 'last') - 1; %may be out slightly
            if shadedError
                slimShadeX = [linspace(0, avgCycleTime, currentSize),linspace(avgCycleTime, 0, currentSize)];
                slimShadeY = [ nanmean(crossFlydX1(:,1:currentSize,block,direction)) + dX1SEM(:,1:currentSize,block,direction) , ...
                    fliplr( nanmean(crossFlydX1(:,1:currentSize,block,direction)) - dX1SEM(:,1:currentSize,block,direction)) ];
                fill(slimShadeX,slimShadeY,directionColour{direction})
                alpha(0.15)
            else
                errorbar(linspace(0, avgCycleTime, currentSize), nanmean(crossFlydX1(:,1:currentSize,block,direction), 1), dX1SEM(:,1:currentSize,block,direction), 'LineStyle', 'none', 'Color', directionColour{direction})
            end
            plot(linspace(0, avgCycleTime, currentSize), nanmean(crossFlydX1(:,1:currentSize,block,direction)), 'LineWidth', 3, 'Color', directionColour{direction})
        end
            set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'Box', 'off')
            ylabel('rotational displacement (mm)')
            xlabel('time(s)')
            % legend(directions)
            title(['dX1, cross-fly block ' num2str(block)])
    end
end

% first code done!