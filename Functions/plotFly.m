% plot many fly things (per block, per volume, collapsed, etc)
function plotFly(R, brainImage, outputDirectory)

%% per block
for b = 1:length(R.BLOCK)
    subDirectory = fullfile(outputDirectory,['Block' num2str(b)]);
    if ~exist(subDirectory,'dir')
       mkdir(subDirectory); 
    end
    plotStuff(R.BLOCK(b), brainImage, subDirectory);
end

%% all blocks concatenated
subDirectory = fullfile(outputDirectory,'All');
if ~exist(subDirectory,'dir')
   mkdir(subDirectory); 
end
plotStuff(R.ALL, brainImage, subDirectory);