close all; clear;
% 
addpath('..\..\Scripts\Global functions\');
addpath('..\..\Scripts\Indexes and legends\');
addpath('..\..\2P SEs\Functions\');
addpath('../../Extracted datasets/');

resultsDirectory = '../../2P Results 3D 2';
% dataDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded blocks
blocks = blocks(~logical(blocks.Exclude),:);

% chosenFlies = 25:35; %cholinergic CC LED
% chosenFlies = [19:24 57:63];% pan-neuronal LED
chosenFlies = [95:102 104:106];

chosenBlocks = {1,1,1,1,1,1,1,1,1,1,1};
% chosenBlocks = {2,2,2,2,2,2,2,2,2,2,2};
% chosenBlocks = {3,3,3,3,3,3,3,3,3,3,3};

totalBlocks = length(cell2mat(chosenBlocks));

imageSize = [32 32];

trim = 3;

imageSize = imageSize-2*trim;

n_pixels = prod(imageSize);

startVol = 1;

slices = 1:4;

%%

FLIES = struct;

for fly = 1:length(chosenFlies)

    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);
    
    % the blocks corresponding to this fly
    thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
        currentBlock = thisFlyBlocks(b,:);
        results = load(fullfile(resultsDirectory,['Fly' num2str(thisFlyBlocks(1,:).Fly)],['Block' num2str(currentBlock.Block)],'results.mat'));

        n_z = thisFlyBlocks(b,:).Steps;

        %maybe remove the first volume or some slices
        results.meanDataSeq = results.meanDataSeq(startVol:end,:,:,:,slices);

        % if blank blocks were use as pedestal, otherwise use mean
        % transient
        if thisFlyBlocks(b,:).BlankBlocks
            results.meanBlankTransient = results.meanBlankTransient(:,:,slices,startVol:end);
            sizeAux = size(results.meanBlankTransient); sizeAux = sizeAux([4 1 2 3]); sizeAux = [sizeAux(1) 1 sizeAux(2:4)];
    
            % normalise each sequence transient by the mean blank transient (i.e. make dF/F)
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanBlankTransient,[4 1 2 3]),sizeAux),[1 16 1 1 1]);
        else
            sizeAux = size(results.meanTransient); sizeAux = sizeAux([4 1 2 3]); sizeAux = [sizeAux(1) 1 sizeAux(2:4)];
            results.meanTransient = results.meanTransient(:,:,slices,startVol:end);
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanTransient,[4 1 2 3]),sizeAux),[1 16 1 1 1]);
        end

        % construct a matrix of activity over time
        transients = permute(squeeze(mean(results.meanDataSeq,2)),[2 3 4 1]);
        
        % remove sides from images
        transients([1:trim end-(trim-1):end],:,:,:) = [];
        transients(:, [1:trim end-(trim-1):end],:,:) = [];

        FLIES(fly).BLOCK(b).transients = transients;

    end
    
end

%% group and average transients

% close all


% this assumes all transient matrices are the same size and were trimmed
% the same way
allTransients = zeros([size(FLIES(1).BLOCK(1).transients) totalBlocks]);

block = 1;

for fly = 1:length(chosenFlies)

    for b = 1:length(FLIES(fly).BLOCK)
        
        allTransients(:,:,:,:,block) = FLIES(fly).BLOCK(b).transients;

        block = block + 1;

    end

end

% calculate mean and plot
meanTransients = squeeze(mean(mean(mean(allTransients,1),2),3));

meanTransient = mean(meanTransients,2);%singular

stdMeanTransient = std(meanTransients,[],2);

semMeanTransient = stdMeanTransient./sqrt(totalBlocks);

errorbar(meanTransient,semMeanTransient); xlim([0 size(meanTransients,1)+1]);

save('before_transient','meanTransient','semMeanTransient');
% save('during_transient','meanTransient','semMeanTransient');
% save('after_transient','meanTransient','semMeanTransient');