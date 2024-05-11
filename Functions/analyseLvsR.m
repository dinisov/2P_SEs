function R = analyseLvsR(R, chosenR, outputDirectory)
%analyseLvsR Summary of this function goes here
%   Detailed explanation goes here

disp('Calculating L vs R');
tic;
% for each block
for fly = 1:length(R)
    for b = [R(fly).BLOCK.blockNum]  

        subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenR(fly))],['Block' num2str(b)],'LvsR');

        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        
        trim = R(fly).BLOCK(b).Trim;
        
        thisBlockData = squeeze(sum(R(fly).BLOCK(b).dataSeqIso,1)./size(R(fly).BLOCK(b).dataSeqIso,1));

        % trim data
        thisBlockData = thisBlockData(:,trim+1:end-trim,trim+1:end-trim,:);

        if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
            blankTrials = R(fly).BLOCK(b).meanBlankTransient;
            blankTrials = blankTrials(trim+1:end-trim,trim+1:end-trim,:);
            meanBlankTrials = mean(blankTrials,3);
        else
            blankTrials = 0;
            meanBlankTrials = 0;
        end

        %% plot means over time
        
        R(fly).BLOCK(b).LvsR  = calculateTtestsLR(thisBlockData);
        
        figure; imagesc(R(fly).BLOCK(b).LvsR.h); saveas(gcf,fullfile(subDirectory,'h.png'));
        figure; imagesc(R(fly).BLOCK(b).LvsR.p); colormap(hot(256));saveas(gcf,fullfile(subDirectory,'p.png'));
        
        LRDiff = R(fly).BLOCK(b).LvsR.LRDiff-meanBlankTrials./(meanBlankTrials+~any(meanBlankTrials(:)));
        figure; imagesc(LRDiff); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'LminusR.png'));
        
        meanL = (R(fly).BLOCK(b).LvsR.meanL-meanBlankTrials)./(meanBlankTrials+~any(meanBlankTrials(:)));
        figure; imagesc(meanL); colormap(jet(256));colorbar; saveas(gcf,fullfile(subDirectory,'meanL.png'));
        
        meanR = (R(fly).BLOCK(b).LvsR.meanR-meanBlankTrials)./(meanBlankTrials+~any(meanBlankTrials(:)));
        figure; imagesc(meanR); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'meanR.png'));
        close all;
        
        %% make movies over time

        left = R(fly).BLOCK(b).dataSeqIso(:,1:16,trim+1:end-trim,trim+1:end-trim,:);
        right = R(fly).BLOCK(b).dataSeqIso(:,17:32,trim+1:end-trim,trim+1:end-trim,:);
        
        meanLeft = sum(left,5)./sum(left~=0,5);
        meanRight = sum(right,5)./sum(right~=0,5);
        
        meanLeft(isnan(meanLeft)) = 0;
        meanRight(isnan(meanRight)) = 0;
        
        meanLeft = permute(squeeze(mean(meanLeft,2)),[2 3 1]);
        meanRight = permute(squeeze(mean(meanRight,2)),[2 3 1]);
        
        makeMovie(prepareMovieData((meanLeft-blankTrials)./(blankTrials+~any(blankTrials(:)))),fullfile(subDirectory,'left_transient.avi'),false);
        makeMovie(prepareMovieData((meanRight-blankTrials)./(blankTrials+~any(blankTrials(:)))),fullfile(subDirectory,'right_transient.avi'),false);

        makeMovie(prepareMovieData((meanRight-meanLeft)./(blankTrials+~any(blankTrials(:)))),fullfile(subDirectory,'difference_transient.avi'),false);
         
    end
end
toc;

function data = prepareMovieData(data)

%     data = data(:,:,2:end); %get rid of first time point
%     data = data-mean(data,3); % normalise by mean along time
%     data = data-data(:,:,end); % normalise by time point

    %map to interval [0,1]
    data = data+abs(min(data,[],'all')); % make minimum 0
    data = data/max(data,[],'all'); %make maximum 1
    
end

end

