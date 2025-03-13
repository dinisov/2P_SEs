% plot many fly things (per block, per volume, collapsed, etc)
function plotFly(R, groupedBlocks, outputDirectory, visibility)

%% per block
%for b = [R.BLOCK.blockNum]
for b = 1:size( R.BLOCK,2 )
    blockNum = R.BLOCK(b).blockNum;
    %subDirectory = fullfile(outputDirectory,['Block' num2str(b)]);
    subDirectory = fullfile(outputDirectory,['Block' num2str(blockNum)]);
    if ~exist(subDirectory,'dir')
       mkdir(subDirectory); 
    end
    plotStuff(R.BLOCK(b), subDirectory,visibility);
end

%% all blocks concatenated
if groupedBlocks
    subDirectory = fullfile(outputDirectory,'All');
    if ~exist(subDirectory,'dir')
       mkdir(subDirectory); 
    end
    plotStuff(R.ALL, subDirectory,visibility);
end