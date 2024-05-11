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
            R(fly).BLOCK(b) = analyse2PBlock(thisBlock);
            
            thisBlockDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]);
            if ~exist(thisBlockDirectory,'dir')
                mkdir(thisBlockDirectory); 
            end
            meanDataSeq = R(fly).BLOCK(b).meanDataSeq;
            meanBlankTransient = R(fly).BLOCK(b).meanBlankTransient;
            meanTransient = R(fly).BLOCK(b).meanTransient;

            save(fullfile(thisBlockDirectory,'results'),'meanDataSeq','meanBlankTransient','meanTransient');
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