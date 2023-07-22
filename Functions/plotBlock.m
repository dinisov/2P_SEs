%BLOCK is a streucture containing
%     meanDataSeq: [nVol×16×sizeX×sizeY double]
%           nData: [1×16 double]
%           r2Vol: [1×nVol struct]
%            rVol: [1×nVol struct]
%              r2: [1×1 struct]
%               r: [1×1 struct]
function plotBlock(BLOCK, brainImage, outputDirectory)

% meanDataSeq = BLOCK.meanDataSeq;
% nVol = size(meanDataSeq,1);
% 
% %% per volume
% subDirectory = fullfile(outputDirectory,'Per volume');
% if ~exist(subDirectory,'dir')
%    mkdir(subDirectory); 
% end
% 
% for vol = 1:nVol
%     thisVolDirectory = fullfile(subDirectory,num2str(vol));
%     if ~exist(thisVolDirectory,'dir')
%        mkdir(thisVolDirectory); 
%     end
% %     thisVolData = permute(squeeze(BLOCK.meanDataSeq(vol,:,:,:)),[2,3,1]);
%     plotStuff(BLOCK, brainImage,thisVolDirectory);
% end

%% collapsed across volumes
subDirectory = fullfile(outputDirectory,'Collapsed');
if ~exist(subDirectory,'dir')
   mkdir(subDirectory); 
end
% thisBlockData = permute(squeeze(sum(BLOCK.meanDataSeq,1)),[2,3,1]);
plotStuff(BLOCK, brainImage, subDirectory);
