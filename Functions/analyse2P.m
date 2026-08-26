function R = analyse2P(FLIES, chosenFlies, outputDirectory, groupedBlocks, saveFull, skipIso)
%analyse2P Summary of this function goes here
%   Detailed explanation goes here

    %structure to put results in
    R = struct;
    
    disp('Analysing SEs...');
    tic;
    for fly = 1:length(FLIES)
        disp(['Fly ' num2str(fly)]);
        thisFly = FLIES(fly);
        %for b = [thisFly.BLOCKS.blockNum]
        for b = 1:size( thisFly.BLOCKS,2 )
            %disp(['Block ' num2str(b)]);
            blockNum = thisFly.BLOCKS(b).blockNum;
            disp(['Block ' num2str(blockNum)])
            thisBlock = thisFly.BLOCKS(b);
            disp(thisBlock)
            %earthquke
            %R(fly).BLOCK(b) = analyse2PBlock(thisBlock);
            R(fly).BLOCK(b) = analyse2PBlock(thisBlock,'skipIso',skipIso);
            
            %thisBlockDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]);
            thisBlockDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(blockNum)]);
            if ~exist(thisBlockDirectory,'dir')
                mkdir(thisBlockDirectory); 
            end
            meanDataSeq = R(fly).BLOCK(b).meanDataSeq;
            meanBlankTransient = R(fly).BLOCK(b).meanBlankTransient;
            meanTransient = R(fly).BLOCK(b).meanTransient;
            if isfield(thisFly.BLOCKS(b),'ancillary')
                ancillary = thisFly.BLOCKS(b).ancillary;
                if isfield(thisFly.BLOCKS(b),'approxVolRate')
                    ancillary.approxVolRate = thisFly.BLOCKS(b).approxVolRate;
                end
            else
                ancillary = struct; %Make empty, put stuff in momentarily, if applicable
            end
            if isfield( thisBlock, 'singularZ' )
                ancillary.singularZ = thisBlock.singularZ;
            end
            if isfield(thisFly.BLOCKS(b),'secondStageSmooth')
                ancillary.secondStageSmooth = thisFly.BLOCKS(b).secondStageSmooth;
                if ancillary.secondStageSmooth == 1
                    ancillary.smooth.smoothWindowTime = thisFly.BLOCKS(b).smoothWindowTime;
                    ancillary.smooth.smoothWindowSize = thisFly.BLOCKS(b).smoothWindowSize;
                    ancillary.smooth.smoothWindowMethod = thisFly.BLOCKS(b).smoothWindowMethod;
                end
            end

            %save(fullfile(thisBlockDirectory,'results'),'meanDataSeq','meanBlankTransient','meanTransient');
            save(fullfile(thisBlockDirectory,'results'),'meanDataSeq','meanBlankTransient','meanTransient', 'ancillary');
            
            %Save isomer data (if requested)
            if ~skipIso
                dataSeqIso = R(fly).BLOCK(b).dataSeqIso;
                %Note: This is the full form
                auxSeq = R(fly).BLOCK(b).auxSeq;

                %%meanDataSeq = sum(dataSeq,5)./sum(dataSeq~=0,5); %This is the math (in analyse2PBlock) that turns dataSeq into meanDataSeq, for reference
                %%meanDataSeq(isnan(meanDataSeq)) = 0;
                    %Note: Empirical testing seems to confirm that under normal circumstances sum and mean yield same shape, but not absolute values

                meanSeqIso = sum(dataSeqIso,5)./sum(dataSeqIso~=0,5);
                meanSeqIso(isnan(meanSeqIso)) = 0;
                    %Note: I do not personally agree with this system to find the mean, but it is Dinis' canonical method

                save(fullfile(thisBlockDirectory,'resultsIso'),'meanSeqIso','auxSeq','meanBlankTransient','ancillary');

                if saveFull
                    save(fullfile(thisBlockDirectory,'resultsIsoFull'),'dataSeqIso','auxSeq','meanBlankTransient','ancillary','-v7.3');
                end

            end

            %Save full (per-event) sequence data, if requested
            if ~isempty(saveFull) && saveFull
                dataSeq = R(fly).BLOCK(b).dataSeq; %Inefficient memory, but necessary for save?

                temp = nansum( squeeze( dataSeq(1,:,1,1,:) ) ~= 0 , 2 ); %Instance count
                    %Note: Leans heavily on nothing being wierd with T:1, pixel 1,1
                if saveFull == 1
                    disp(['Saving (reduced) full sequence data, as requested'])
                elseif saveFull == 2
                    disp(['Saving (actual) full sequence data, as requested'])
                end
                disp(['# of instances of seqs:'])
                disp( temp )

                dataSeqReduced = nan( [size( dataSeq, [1:4] ), nanmax(temp) ] );
                    %Save memory, at the cost of obliterating true time
                    %Architecture: T*Seq*X*Y*Event#
                dataSeqReducedInd = nan( size( dataSeq, 2 ) , nanmax(temp) );
                    %Might as well store original event # alongside
                for seq = 1:size( dataSeq, 2 )
                    temp = find( dataSeq(1,seq,1,1,:) ~= 0 );
                    dataSeqReduced( :, seq, :, :, [1:numel(temp)] ) = ...
                        dataSeq(:,seq,:,:, temp );
                    dataSeqReducedInd( seq, [1:numel(temp)] ) = temp;
                end

                %And save
                if saveFull == 1
                    save(fullfile(thisBlockDirectory,'resultsFull'),'dataSeqReduced','dataSeqReducedInd',...
                        '-v7.3');
                        %Since dataSeq is very big, we will not save it unless absolutely necessary
                            %Like, this takes multiple minutes per individual block
                elseif saveFull == 2
                    tic
                    save(fullfile(thisBlockDirectory,'resultsFull'),'dataSeq','dataSeqReduced','dataSeqReducedInd',...
                        '-v7.3');
                    toc
                end


            end
            
            %Check for behav state data
            if isfield( R(fly).BLOCK(b), 'dataSeqBehav' ) && ~isempty(R(fly).BLOCK(b).dataSeqBehav)
                %disp(R(fly).BLOCK(b).dataSeqBehav)
                for stat = 1:size(R(fly).BLOCK(b).dataSeqBehav,2)
                    %disp( R(fly).BLOCK(b).dataSeqBehav(stat) )
                    %R(fly).BLOCK(b).dataSeqBehav(stat).state
                    thisState = R(fly).BLOCK(b).dataSeqBehav(stat).state;
                    %dataSeq = R(fly).BLOCK(b).dataSeqBehav(stat).dataSeq; %Not in use currently cos filled with unnecessary zeroes
                    dataSeqBehavReduced = R(fly).BLOCK(b).dataSeqBehav(stat).dataSeqReduced;
                    %meanDataSeqBehavReduced = R(fly).BLOCK(b).dataSeqBehav(stat).meanDataSeqReduced;
                    meanDataSeqBehav = R(fly).BLOCK(b).dataSeqBehav(stat).meanDataSeqReduced;
                    %disp( num2str(size(dataSeq)) )
                    thisDir = [thisBlockDirectory,filesep,'State_',num2str(R(fly).BLOCK(b).dataSeqBehav(stat).state)];
                    if exist(thisDir) ~= 7
                        mkdir(thisDir);                        
                    end
                    %disp(dataSeq)
                    tic
                    %save( [thisDir,filesep,'stateResults'],...
                    %    'thisState','dataSeqBehavReduced', 'meanDataSeqBehavReduced', '-v7.3'); %MATLAB version specification necessary in case dataSeq is large
                        %Note: meanDataSeqReduced is *vastly* smaller than dataSeqBehavReduced, so it might be good to save as separate
                    %save( [thisDir,filesep,'stateResults'],...
                    %    'meanBlankTransient','thisState','meanDataSeqBehav'); %Non-trial
                    %save( [thisDir,filesep,'stateResultsExtendedReduced'],...
                    %    'meanBlankTransient','thisState','dataSeqBehavReduced', '-v7.3'); %W/ trial
                    save( [thisDir,filesep,'stateResults'],...
                        'meanBlankTransient','thisState','meanDataSeqBehav', 'ancillary'); %Non-trial
                    save( [thisDir,filesep,'stateResultsExtendedReduced'],...
                        'meanBlankTransient','thisState','dataSeqBehavReduced', 'ancillary', '-v7.3'); %W/ trial
                    disp(['Saved data for behav state ',num2str(thisState),' to ',thisDir,' in ',num2str(toc),'s'])
                    %Make a little figure to see where this state occurred in time
                    if isfield( R(fly).BLOCK(b).dataSeqBehav(stat), 'behavSeq' )
                        behavSeq = R(fly).BLOCK(b).dataSeqBehav(stat).behavSeq;
                        figure
                        plot( behavSeq, 'Color', 'k' )
                        if length( unique(behavSeq) ) > 1
                            ylim([nanmin(behavSeq)*1.1, nanmax(behavSeq)*1.1])
                        else
                            disp(['-# Only one state present in plot data #-'])
                        end
                        hold on
                        behavSeq( behavSeq ~= thisState ) = NaN;
                        plot( behavSeq , 'Color', 'g' )
                        xlabel('Time')
                        ylabel('State')
                        title(['State\_',num2str(thisState)])
                        saveas(gcf,fullfile(thisDir,['StateSeqFigure.png']));
                    end
                end
            end
            
        end
        
        % add brain images to results structure
        %(need to do this separately to avoid dissimilar structures)
        %for b = [thisFly.BLOCKS.blockNum]
        for b = 1:size( thisFly.BLOCKS,2 )
            R(fly).BLOCK(b).brainImage = FLIES(fly).BLOCKS(b).brainImage;
        end
        
        %THIS NEEDS WORK
        %{
        if groupedBlocks
            disp(['Fly ' num2str(fly) ' grouped blocks']);
            R(fly).ALL = analyseSequentialEffectsTwoPhoton(thisFly.greenChannel,thisFly.randomSequence,thisFly.nVol);
            R(fly).ALL.brainImage = FLIES(1).BLOCKS(b).brainImage;
            save(fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],'All','results'),'results');%UNTESTED
        end
        %}
    end
    toc;

end