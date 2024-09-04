%Reduce Files wrapper
%   Designed to loop reduceFilesFunc for as many fly/blocks requested

clear
close all

rdmDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");

%get rid of excluded flies
% blocks = blocks(~logical(blocks.Exclude),:);

% the numbers here should be the original size divided by some power of 2
finalSize = [128 128];

chosenFlies = [224];
chosenBlocks = {[1,2,3,4]}; % leave empty if reducing all blocks for one fly
    %MUST BE IN FORMAT {[blocks]}
    
%%

for fly = 1:length(chosenFlies)
    thisFly = chosenFlies(fly);
    
   % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == thisFly,:);

    if ~isempty(chosenBlocks) && ~isempty(chosenBlocks{fly})
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
    end

    nBlocks = height(thisFlyBlocks);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    parfor b = 1:nBlocks
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate]; %Borrowed from pre_process
        currentRDMDirectory = fullfile(rdmDirectory,currentDate,flyID);
        try
            reduceFilesFunc( currentRDMDirectory, currentBlock, flyID, finalSize )
        catch
            disp(['-# Failure to process ',flyID,' #-'])
        end
        
    end
    
    
end