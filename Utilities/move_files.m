%

close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

scratchDirectory = '../../2P Data';

flies = readtable('../../2P Record/2P_record');

flies = flies(~logical(flies.Exclude),:);

%%

% this level is the list of dates
for block = 51:height(flies)
    
    currentBlock = flies(block,:);
    
    currentDate = char(datetime(currentBlock.Date,'Format','dMMMyy'));

    currentFlyDirectory = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];

    currentDirectory = fullfile(mainDirectory,currentDate,currentFlyDirectory);

    matFiles = dir([currentDirectory '\*.mat']);
    
    if ~exist(fullfile(scratchDirectory,currentDate,currentFlyDirectory),'dir')    
        mkdir(fullfile(scratchDirectory,currentDate,currentFlyDirectory));
    end
    
    tic;
    for i = 1:length(matFiles) 
        currentFile = matFiles(i).name;
        copyfile(fullfile(mainDirectory,currentDate,currentFlyDirectory,currentFile),fullfile(scratchDirectory,currentDate,currentFlyDirectory,currentFile));
    end
    toc;
    
end