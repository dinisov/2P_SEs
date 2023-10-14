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
           thisBlockData = squeeze(mean(R(fly).BLOCK(b).meanDataSeq,1));
           AAAA = squeeze(thisBlockData(16,trim+1:end-trim,trim+1:end-trim)); AAAR = squeeze(thisBlockData(8,trim+1:end-trim,trim+1:end-trim)); 
           RRRR = squeeze(thisBlockData(1,trim+1:end-trim,trim+1:end-trim)); RRRA = squeeze(thisBlockData(9,trim+1:end-trim,trim+1:end-trim));

           figure; imagesc(AAAA); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAA.png')); 
           figure; imagesc(RRRR); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRR.png'));
           figure; imagesc(RRRA); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRA.png'));
           figure; imagesc(AAAR); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAR.png'));
           figure; imagesc((AAAA-AAAR)); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAAminusAAAR.png'));
           figure; imagesc((RRRR-RRRA)); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRRminusRRRA.png'));
           close all;

           %videos of differences over time

           %AAAA-AAAR
           AAAAminusAAAR_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,16,trim+1:end-trim,trim+1:end-trim)),[2 3 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,8,trim+1:end-trim,trim+1:end-trim)),[2 3 1]);
           RRRRminusRRRA_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,1,trim+1:end-trim,trim+1:end-trim)),[2 3 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,9,trim+1:end-trim,trim+1:end-trim)),[2 3 1]);

           AAAAminusAAAR_t = normalize(AAAAminusAAAR_t,3,'range',[0 1]);
           makeMovie(AAAAminusAAAR_t,fullfile(subDirectory,'AAAAminusAAAR_transient.avi'),false);

           RRRRminusRRRA_t = normalize(RRRRminusRRRA_t,3,'range',[0 1]);
           makeMovie(RRRRminusRRRA_t,fullfile(subDirectory,'RRRRminusRRRA_transient.avi'),false);

        end
    end
end

