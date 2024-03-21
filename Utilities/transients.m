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

chosenFlies = 80;

trim = 3;

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
        
        seqStruct = struct;
        
        % global response transient
        allSeq = permute(squeeze(sum(results.meanDataSeq,2)),[2 3 1]);
        allSeq = allSeq(trim+1:end-trim,trim+1:end-trim,:);

        makeMovie(allSeq,fullfile(subDirectory,'global.avi'));
        
        % response transient per sequence
        for s = 1:16     
            seq = permute(squeeze(results.meanDataSeq(:,s,:,:)),[2 3 1]);
            seq = seq(trim+1:end-trim,trim+1:end-trim,:);
            
            makeMovie(seq,fullfile(subDirectory,['seq' num2str(s) '.avi']));
            seqStruct(s).data = seq;
        end
        
        %Analyse and plot
        seqMeans = nan(1,size(seqStruct(1).data,3));
        for s = 1:16
            seqMeans(s,:) = squeeze( nanmean(nanmean(seqStruct(s).data,1),2) )';
        end
        
        figure
        hold on
        for s = 1:16
           plot(seqMeans(s,:)) 
        end
        legend
        
    end
        
end

function makeMovie(data, filename)

    VidObj = VideoWriter(filename, 'Uncompressed AVI');
    
    data = data(:,:,2:end); %get rid of first time point
    data = data-mean(data,3); % normalise by mean along time
%     data = data-data(:,:,end); % normalise by time point

    %map to interval [0,1]
    data = data+abs(min(data,[],'all')); % make minimum 0
    data = data/max(data,[],'all'); %make maximum 1

    VidObj.FrameRate = 2.5; %set your frame rate
    open(VidObj);
    for f = 1:size(data, 3)
        writeVideo(VidObj, ind2rgb(uint8(255 * data(:,:,f)), jet(256)));
    end
    close(VidObj);
        
end