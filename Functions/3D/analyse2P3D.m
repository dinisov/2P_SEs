function R = analyse2P3D(FLIES, chosenFlies, outputDirectory, groupedBlocks)
%analyse2P Summary of this function goes here
%   Detailed explanation goes here

    %structure to put results in
    R = struct;
    
    disp('Analysing SEs...');
    tic;
    for fly = 1:length(FLIES)
        disp(['Fly ' num2str(fly)]);
        thisFly = FLIES(fly);
        for b = 1:length(thisFly.BLOCKS)
            thisBlock = thisFly.BLOCKS(b);
            thisBlockID = thisFly.BLOCKS(b).thisBlockID;
            disp(['Analysing block ' num2str(b),' of ',num2str(length(thisFly.BLOCKS)),' (ID: ',num2str(thisBlockID),')']);
            R(fly).BLOCK(b) = analyse2PBlock3D(thisBlock);
            %R(fly).BLOCK(b).thisBlockID = [thisBlockID]; %Causes crashes due to dissimilar sizes with above function output
            
            %thisBlockDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]); %Dinis system, does not account for gaps
            thisBlockDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(thisBlockID)]); %Matt new, skips excluded block IDs
            if ~exist(thisBlockDirectory,'dir')
                mkdir(thisBlockDirectory); 
            end
            
            meanDataSeq = R(fly).BLOCK(b).meanDataSeq;
            meanBlankTransient = R(fly).BLOCK(b).meanBlankTransient;
            meanTransient = R(fly).BLOCK(b).meanTransient;

            save(fullfile(thisBlockDirectory,'results'),'meanDataSeq','meanBlankTransient','meanTransient');
            
            dataSeq = R(fly).BLOCK(b).dataSeq;
            inds = find( sum( dataSeq(:,:,:,:,:,:), [1,2,3,4,5] ) ~= 0 ); %Find only non-zero trials
            dataSeqReduced = dataSeq(:,:,:,:,:,inds);
            %save(fullfile(thisBlockDirectory,'results_extended_reduced'),'dataSeqReduced','-v7.3');
            dataStruct = struct;
            dataStruct.dataSeqReduced = dataSeqReduced;
            dataStruct.flyName = chosenFlies(fly);
            %if isfield(thisFly.BLOCKS,'thisBlock') == 1
            dataStruct.blockName = thisBlockID;
            %end
            save(fullfile(thisBlockDirectory,'results_extended_reduced'),'dataStruct','-v7.3');
            disp(['Extended transient results saved'])
        end
        
        % add brain images to results structure
        for b = 1:length(thisFly.BLOCKS)
            R(fly).BLOCK(b).brainImage = FLIES(fly).BLOCKS(b).brainImage;
            thisBlockID = thisFly.BLOCKS(b).thisBlockID;
            R(fly).BLOCK(b).thisBlockID = thisBlockID; %And append this too
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