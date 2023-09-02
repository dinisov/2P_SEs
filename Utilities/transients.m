close all; clear;

addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_swinderen\Dinis\2P SEs\Functions\');

close all; clear;

resultsDirectory = '../../2P Results';

outputDirectory = '../../2P Clusters';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [32 32];

flyList = unique(blocks.Fly);

chosenFlies = 1:41;

%%
for fly = chosenFlies
   
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == fly,:);
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
        
        thisBlockDirectory = fullfile(resultsDirectory,['Fly' num2str(thisFlyBlocks(1,:).Fly)],['Block' num2str(b)]);
    
        results = load(fullfile(thisBlockDirectory,'results.mat'));
        
        subDirectory = fullfile(thisBlockDirectory,'Transients');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        
        % global response transient
        allSeq = permute(squeeze(sum(results.meanDataSeq,2)),[2 3 1]);
        allSeq = allSeq/max(allSeq,[],'all');
        makeMovie(allSeq,fullfile(subDirectory,'global.avi'));
        
        % response transient per sequence
        for s = 1:16     
            seq = permute(squeeze(results.meanDataSeq(:,s,:,:)),[2 3 1]);
            seq = seq/max(seq,[],'all');
            makeMovie(seq,fullfile(subDirectory,['seq' num2str(s) '.avi']));
        end
        
    end
        
end

function makeMovie(data, filename)

    VidObj = VideoWriter(filename, 'Uncompressed AVI');

    VidObj.FrameRate = 2.5; %set your frame rate
    open(VidObj);
    for f = 1:size(data, 3)
        writeVideo(VidObj, ind2rgb(uint8(255 * data(:,:,f)), parula(256)));
    end
    close(VidObj);
        
end