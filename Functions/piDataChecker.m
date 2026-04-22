function piDataChecker(FLIES,flyRecord,chosenFlies,sequenceDirectory)
%Designed to check for obvious errors in CSV usage for analysis

for fly = 1:length(chosenFlies)
    
    thisFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(fly),:);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    nBlocks = height(thisFlyBlocks);

    for b = 1:nBlocks
        
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        
        if exist(fullfile(sequenceDirectory,[flyID '.mat']),'file')
            randomSequence = load(fullfile(sequenceDirectory,[flyID '.mat']),'randomSequence').';
            randomSequence = randomSequence.randomSequence-1;
        else
            randomSequence = csvread(fullfile(sequenceDirectory,'Data_LEDs',[flyID '.csv'])).';
        end

    end

end

end