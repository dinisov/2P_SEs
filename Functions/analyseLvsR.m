function R = analyseLvsR(R, FLIES, chosenFlies, outputDirectory, trim)
%analyseLvsR Summary of this function goes here
%   Detailed explanation goes here

% newGridSize = gridSize-2*trim;

disp('Calculating L vs R');
tic;
% for each block
for fly = 1:length(FLIES)
    for b = 1:length(FLIES(fly).BLOCKS)
        subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'LvsR');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        thisBlockData = squeeze(mean(R(fly).BLOCK(b).dataSeqIso,1));
        R(fly).BLOCK(b).LvsR  = calculateTtestsLR(thisBlockData);
        figure; imagesc(R(fly).BLOCK(b).LvsR.h); saveas(gcf,fullfile(subDirectory,'h.png'));
        figure; imagesc(R(fly).BLOCK(b).LvsR.p); colormap(hot(256));saveas(gcf,fullfile(subDirectory,'p.png'));
        
        LRDiff = R(fly).BLOCK(b).LvsR.LRDiff(trim+1:end-trim,trim+1:end-trim); %LRDiff = normalize(LRDiff,'range',[0 1]);
        figure; imagesc(LRDiff); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'LminusR.png'));
        
        meanL = R(fly).BLOCK(b).LvsR.meanL(trim+1:end-trim,trim+1:end-trim);% meanL = normalize(meanL,'range',[0 1]);
        figure; imagesc(meanL); colormap(jet(256));colorbar; saveas(gcf,fullfile(subDirectory,'meanL.png'));
        
        meanR = R(fly).BLOCK(b).LvsR.meanR(trim+1:end-trim,trim+1:end-trim);% meanR = normalize(meanR,'range',[0 1]);
        figure; imagesc(meanR); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'meanR.png'));
        close all;
        
        %make movies for L and R transients and difference
        left = permute(squeeze(mean(mean(R(fly).BLOCK(b).dataSeqIso(:,1:2:31,trim+1:end-trim,trim+1:end-trim,:),5),2)),[2 3 1]);
        makeMovie(prepareMovieData(left),fullfile(subDirectory,'left_transient.avi'),false);
        
        right = permute(squeeze(mean(mean(R(fly).BLOCK(b).dataSeqIso(:,2:2:32,trim+1:end-trim,trim+1:end-trim,:),5),2)),[2 3 1]);
        makeMovie(prepareMovieData(right),fullfile(subDirectory,'right_transient.avi'),false);
        
        diffLR = left-right;
%         diffLR = normalize(diffLR,3,'range',[0 1]);
        makeMovie(prepareMovieData(diffLR),fullfile(subDirectory,'difference_transient.avi'),false);
         
    end
end
toc;

function data = prepareMovieData(data)

    data = data(:,:,2:end); %get rid of first time point
    data = data-mean(data,3); % normalise by mean along time
%     data = data-data(:,:,end); % normalise by time point

    %map to interval [0,1]
    data = data+abs(min(data,[],'all')); % make minimum 0
    data = data/max(data,[],'all'); %make maximum 1
    
end

end

