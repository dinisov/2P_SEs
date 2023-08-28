% plot many fly things (per block, per volume, collapsed, etc)
function plotFlyBatch(R, iter, outputDirectory)

%% per block
for b = 1:length(R.BLOCK)
    subDirectory = fullfile(outputDirectory,['Block' num2str(b)]);
    if ~exist(subDirectory,'dir')
       mkdir(subDirectory); 
    end
    plotStuffBatch(R.BLOCK(b), iter, subDirectory);
end

%% all blocks concatenated
% if groupedBlocks
%     subDirectory = fullfile(outputDirectory,'All');
%     if ~exist(subDirectory,'dir')
%        mkdir(subDirectory); 
%     end
%     plotStuff(R.ALL, subDirectory);
% end