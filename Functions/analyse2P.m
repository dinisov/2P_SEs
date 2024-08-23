function R = analyse2P(FLIES, chosenFlies, outputDirectory, groupedBlocks)
%analyse2P Summary of this function goes here
%   Detailed explanation goes here

    %structure to put results in
    R = struct;
    
    disp('Analysing SEs...');
    tic;
    for fly = 1:length(FLIES)
        disp(['Fly ' num2str(fly)]);
        thisFly = FLIES(fly);
        for b = [thisFly.BLOCKS.blockNum]
            disp(['Block ' num2str(b)]);
            thisBlock = thisFly.BLOCKS(b);
            disp(thisBlock)
            %earthquke
            R(fly).BLOCK(b) = analyse2PBlock(thisBlock);
            
            thisBlockDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]);
            if ~exist(thisBlockDirectory,'dir')
                mkdir(thisBlockDirectory); 
            end
            meanDataSeq = R(fly).BLOCK(b).meanDataSeq;
            meanBlankTransient = R(fly).BLOCK(b).meanBlankTransient;
            meanTransient = R(fly).BLOCK(b).meanTransient;

            save(fullfile(thisBlockDirectory,'results'),'meanDataSeq','meanBlankTransient','meanTransient');
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
                    save( [thisDir,filesep,'stateResults'],...
                        'meanBlankTransient','thisState','meanDataSeqBehav'); %Non-trial
                    save( [thisDir,filesep,'stateResultsExtendedReduced'],...
                        'meanBlankTransient','thisState','dataSeqBehavReduced', '-v7.3'); %W/ trial
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
        for b = [thisFly.BLOCKS.blockNum]
            R(fly).BLOCK(b).brainImage = FLIES(fly).BLOCKS(b).brainImage;
        end
        
        %THIS NEEDS WORK
        if groupedBlocks
            disp(['Fly ' num2str(fly) ' grouped blocks']);
            R(fly).ALL = analyseSequentialEffectsTwoPhoton(thisFly.greenChannel,thisFly.randomSequence,thisFly.nVol);
            R(fly).ALL.brainImage = FLIES(1).BLOCKS(b).brainImage;
            save(fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],'All','results'),'results');%UNTESTED
        end
    end
    toc;

end