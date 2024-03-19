% register blocks of SE trials within and across blocks
% currently does the Z averaged stacks but can be extended to align the
% original green/red channels (i.e. align each slice in the Z stack)

close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

% scratchDirectory = '../../2P Data';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded flies
% blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
imageSize = [128 128];

chosenFlies = [20];

chosenBlocks = [1];

% chosenFlies = [4 5 6 7 13 20 22 23 38 50 54];
% 
% chosenBlocks = {[1 3],1,2,[1 2],2,1,3,2,2,2,[2 3]};

%%

for fly = chosenFlies

    thisFlyBlocks = blocks(blocks.Fly == fly,:);

    if ~isempty(chosenBlocks)
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks),:);
    end
    
    nBlocks = height(thisFlyBlocks);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    %align inside each block
    for b = 1:nBlocks
        
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        currentDirectory = fullfile(mainDirectory,currentDate,flyID);
        
        nSlices = currentBlock.Steps + currentBlock.FlybackFrames;
        
        realFrames = currentBlock.realFrames;
        idealFrames = ceil(realFrames/100)*100;% frame number to nearest 100
        
        missingVols = (idealFrames-realFrames)/nSlices;
        
        if mod(missingVols,1) ~= 0
           disp('NON-INTEGER NUMBER OF VOLUMES MISSING!!');
        end
        
        nTrials = currentBlock.BlockLength + currentBlock.BlankBlocks;
        
        volPerTrial = idealFrames/(nTrials*nSlices);
        
%         load(fullfile(currentDirectory, 'avg_z_green_aligned'),'avg_z_green_aligned');

        load(fullfile(currentDirectory, 'green_channel_aligned'),'green_channel_aligned'); 
        
        green_channel_aligned = squeeze(green_channel_aligned(:,:,1,:));

        trace = squeeze(mean(mean(green_channel_aligned,1),2));
        
        volsMissing = [1];

        badTrialList = [309];
        
        %WARNING: this assumes only one volume is missing per trial
        for i = 1:length(volsMissing)
           badTrial = findBadTrial(trace,nTrials-(i-1),volPerTrial);
%            trace = removeBadTrial(badTrialList(i)-(i-1),trace,volPerTrial,volsMissing(i));
        end
        
    end

end

% this uses the first 3 trials of each block because it is a strong motif
% due to a dip in brightness during the first volume of each block
% there is a small chance the error was during the first 3 volumes in which
% case badTrial is not the current but the next one; inspect on a
% case-by-case basis
function badTrial = findBadTrial(trace,nTrials,volPerTrial)
    badTrial = NaN;
    for i = 0:nTrials-2
        thisTrial = trace(i*volPerTrial + 1:(i+1)*volPerTrial);
%         nextTrial = trace((i+1)*volPerTrial + 1:(i+2)*volPerTrial);

% correlation coefficient method
%         cc = corrcoef(thisTrial(1:2),nextTrial(1:2));
%         if cc(2) < 0.9
%            disp(['trial ' num2str(i+1) ' is bad']);
%            figure; plot([thisTrial nextTrial]); legend({'This trial','Next trial'});
%            badTrial = i+1;
%            break;
%         end

        % index of minimum method
        [~,ind] = min(thisTrial);
        if ind ~= 1
           disp(['trial ' num2str(i+1) ' is bad']);
           figure; plot([thisTrial]); %legend({'This trial','Next trial'});
           badTrial = i+1;
           break;
        end

    end
end

function trace = removeBadTrial(badTrial,trace,volPerTrial,volsMissing)
    badTrial = badTrial - 1;
    disp((badTrial*volPerTrial) + 1:(badTrial*volPerTrial) + volPerTrial - volsMissing);
    trace((badTrial*volPerTrial) + 1:(badTrial*volPerTrial) + volPerTrial - volsMissing) = [];
end