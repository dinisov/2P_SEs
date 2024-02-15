function patternPlots3D(R, FLIES, chosenFlies, outputDirectory)
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
            
           %data averaged over time
           thisBlockData = squeeze(mean(R(fly).BLOCK(b).meanDataSeq,1));
           
           % trim data
           thisBlockData = thisBlockData(:,trim+1:end-trim,trim+1:end-trim,:);
            
            if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
              
                blankTrials = R(fly).BLOCK(b).meanBlankTransient;
                blankTrials = blankTrials(trim+1:end-trim,trim+1:end-trim,:,:);
                meanBlankTrials = mean(blankTrials,4);
                plot3D(meanBlankTrials,'off');
                saveas(gcf,fullfile(subDirectory,'average_blanks.png')); 
                
                % subtract and divide by blank trial average (dF/F)
                aux = permute(repmat(meanBlankTrials,[1 1 1 16]),[4 1 2 3]);
                thisBlockData = (thisBlockData-aux)./aux;
                 
            end
           
           AAAA = squeeze(thisBlockData(16,:,:,:)); AAAR = squeeze(thisBlockData(8,:,:,:)); 
           RRRR = squeeze(thisBlockData(1,:,:,:)); RRRA = squeeze(thisBlockData(9,:,:,:));
           
           load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
           
           %this is just to help turn horizontal sequences into vertical ones
           ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);
           
           for i = 1:16
                plot3D(squeeze(thisBlockData(i,:,:,:)),'off');
                saveas(gcf,fullfile(subDirectory,[binomial_x_labels_latex{i}(ind_horiz) '.png'])); 
           end

           plot3D((AAAA-AAAR),'off'); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'AAAAminusAAAR.png'));
           plot3D((RRRR-RRRA),'off'); colormap(jet(256)); colorbar; saveas(gcf,fullfile(subDirectory,'RRRRminusRRRA.png'));
           close all;

           %videos of differences over time

           %AAAA-AAAR and RRRR-RRRA over time
           AAAAminusAAAR_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,16,trim+1:end-trim,trim+1:end-trim,:)),[2 3 4 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,8,trim+1:end-trim,trim+1:end-trim,:)),[2 3 4 1]);
           RRRRminusRRRA_t = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,1,trim+1:end-trim,trim+1:end-trim,:)),[2 3 4 1]) - permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,9,trim+1:end-trim,trim+1:end-trim,:)),[2 3 4 1]);

           makeMovie3D(prepareMovieData(AAAAminusAAAR_t),fullfile(subDirectory,'AAAAminusAAAR_transient.avi'),false);
           makeMovie3D(prepareMovieData(RRRRminusRRRA_t),fullfile(subDirectory,'RRRRminusRRRA_transient.avi'),false);

        end
    end
end

function data = prepareMovieData(data)

%     data = data(:,:,2:end); %get rid of first time point
%     data = data-mean(data,3); % normalise by mean along time
%     data = data-data(:,:,end); % normalise by time point

    %map to interval [0,1]
    data = data - min(data,[],'all'); % make minimum 0
    data = data / max(data,[],'all'); %make maximum 1
    
end