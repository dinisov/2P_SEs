close all; clear;
% 
addpath('..\..\Scripts\Global functions\');
addpath('..\..\Scripts\Indexes and legends\');
addpath('..\..\2P SEs\Functions\');
addpath('../../Extracted datasets/');

resultsDirectory = '../../2P Results';
% dataDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded blocks
blocks = blocks(~logical(blocks.Exclude),:);

% chosenFlies = 25:35; %cholinergic CC LED
% chosenFlies = [19:24 57:63];% pan-neuronal LED
chosenFlies = [111 112];

chosenBlocks = {1,1};

totalBlocks = length(cell2mat(chosenBlocks));

imageSize = [32 32];

trim = 3;

imageSize = imageSize-2*trim;

n_pixels = prod(imageSize);

% load six_hertz.mat

%%

FLIES = struct;

for fly = 1:length(chosenFlies)
    
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);
    
    % the blocks corresponding to this fly
    thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
    
    % at some point block index should migrate to real block number, right
    % now it is an index
    for b = 1:nBlocks
    
        results = load(fullfile(resultsDirectory,['Fly' num2str(thisFlyBlocks(1,:).Fly)],['Block' num2str(b)],'results.mat'));

        %remove first time point (often it is anomalous)
        results.meanDataSeq = results.meanDataSeq(1:end,:,:,:);

        % if blank blocks were use as pedestal, otherwise use mean
        % transient
        if thisFlyBlocks(b,:).BlankBlocks
            sizeAux = size(results.meanBlankTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
    
            % normalise each sequence transient by the mean blank transient (i.e. make dF/F)
            results.meanDataSeq = results.meanDataSeq./repmat(reshape(permute(results.meanBlankTransient,[3 1 2]),sizeAux),[1 16 1 1]);
        else
            sizeAux = size(results.meanTransient); sizeAux = sizeAux([3 1 2]); sizeAux = [sizeAux(1) 1 sizeAux(2:3)];
            results.meanTransient = results.meanTransient(:,:,1:end);
        end

        transients = permute(squeeze(mean(results.meanDataSeq,2)),[2 3 1]);
            
        % remove sides from images
        transients([1:trim end-(trim-1):end],:,:) = [];
        transients(:, [1:trim end-(trim-1):end],:) = [];

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
        
        allTransients(:,:,block) = FLIES(fly).BLOCK(b).transients;

        block = block + 1;

    end

end