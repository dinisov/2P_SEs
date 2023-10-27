function patternPlots(R, FLIES, chosenFlies, outputDirectory)
%patternPlots Summary of this function goes here
%   Detailed explanation goes here
    trim = 3;

    % newGridSize = gridSize-2*trim;

    for fly = 1:length(FLIES)
        for b = 1:length(FLIES(fly).BLOCKS)
            subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'Oddballs');
            if ~exist(subDirectory,'dir')
               mkdir(subDirectory); 
            end
            
            if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
                blankTrials = R(fly).BLOCK(b).meanBlankTransient;
                blankTrials = blankTrials(trim+1:end-trim,trim+1:end-trim,:);
                meanBlankTrials = mean(blankTrials,3);
            else
                blankTrials = 0;
                meanBlankTrials = 0;
            end
            
           thisBlockData = squeeze(mean(R(fly).BLOCK(b).meanDataSeq,1));
           
           AAAA = squeeze(thisBlockData(16,trim+1:end-trim,trim+1:end-trim)); AAAR = squeeze(thisBlockData(8,trim+1:end-trim,trim+1:end-trim)); 
           RRRR = squeeze(thisBlockData(1,trim+1:end-trim,trim+1:end-trim)); RRRA = squeeze(thisBlockData(9,trim+1:end-trim,trim+1:end-trim));

           figure; imagesc(AAAA-meanBlankTrials); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAA.png')); 
           figure; imagesc(RRRR-meanBlankTrials); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRR.png'));
           figure; imagesc(RRRA-meanBlankTrials); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRA.png'));
           figure; imagesc(AAAR-meanBlankTrials); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAR.png'));
           figure; imagesc((AAAA-AAAR)); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAAminusAAAR.png'));
           figure; imagesc((RRRR-RRRA)); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRRminusRRRA.png'));
           close all;

           %videos of differences over time

           %AAAA-AAAR
           AAAAminusAAAR_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,16,trim+1:end-trim,trim+1:end-trim)),[2 3 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,8,trim+1:end-trim,trim+1:end-trim)),[2 3 1]);
           RRRRminusRRRA_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,1,trim+1:end-trim,trim+1:end-trim)),[2 3 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,9,trim+1:end-trim,trim+1:end-trim)),[2 3 1]);

%            AAAAminusAAAR_t = normalize(AAAAminusAAAR_t,3,'range',[0 1]);
           makeMovie(prepareMovieData(AAAAminusAAAR_t),fullfile(subDirectory,'AAAAminusAAAR_transient.avi'),false);

%            RRRRminusRRRA_t = normalize(RRRRminusRRRA_t,3,'range',[0 1]);
           makeMovie(prepareMovieData(RRRRminusRRRA_t),fullfile(subDirectory,'RRRRminusRRRA_transient.avi'),false);

        end
    end
end

function data = prepareMovieData(data)

%     data = data(:,:,2:end); %get rid of first time point
%     data = data-mean(data,3); % normalise by mean along time
%     data = data-data(:,:,end); % normalise by time point

    %map to interval [0,1]
    data = data+abs(min(data,[],'all')); % make minimum 0
    data = data/max(data,[],'all'); %make maximum 1
    
end