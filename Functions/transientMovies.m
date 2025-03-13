function transientMovies(R, chosenFlies, resultsDirectory)
%transientMovies Summary of this function goes here
%   Detailed explanation goes here

% load six_hertz.mat

for fly = 1:length(R)
    
    %for b = [R(fly).BLOCK.blockNum]
    for b = 1:size( R(fly).BLOCK,2 )
        blockNum = R(fly).BLOCK(b).blockNum;

        %thisBlockDirectory = fullfile(resultsDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]);
        thisBlockDirectory = fullfile(resultsDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(blockNum)]);
        
        subDirectory = fullfile(thisBlockDirectory,'Transients');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        
        %m_
        if ~isnan(R(fly).BLOCK(b).Trim) && ( ~isfield(R(fly).BLOCK(b),'TrimCoords') || isempty(R(fly).BLOCK(b).TrimCoords) )
            %trim = R(fly).BLOCK(b).Trim;
            trim = repmat( R(fly).BLOCK(b).Trim , 1 , 4 ); %Make coord-like, to simplify later
            disp(['Using unitary trim for transient movie creation'])
        else
            trim = R(fly).BLOCK(b).TrimCoords; %Note matrix, not singular
            disp(['Using coordinate trim for transient movie creation'])
        end
        
        % blank trial transient movie
        if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
            
            blankTrials = R(fly).BLOCK(b).meanBlankTransient;
            %blankTrials = blankTrials(trim+1:end-trim,trim+1:end-trim,:);
            blankTrials = blankTrials(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:);
            
            makeMovie(prepareMovieData(blankTrials),fullfile(subDirectory,'blankTrials.avi'),false);
            %makeMovie(prepareMovieData(blankTrials),fullfile(subDirectory,'blankTrials.mp4'),false);
            
            % global response transient
            allSeq = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,2)),[2 3 1]);
            %allSeq = (allSeq(trim+1:end-trim,trim+1:end-trim,:)-blankTrials)./blankTrials;
            allSeq = (allSeq(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:)-blankTrials)./blankTrials;
        else
            % global response transient
            allSeq = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,2)),[2 3 1]);
            %allSeq = allSeq(trim+1:end-trim,trim+1:end-trim,:);
            allSeq = allSeq(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:);
        end
        
        save(fullfile(subDirectory,'global_transient'),'allSeq');
                
%         sizeSeq = size(allSeq);
% 
%         seqAux = reshape(allSeq,[prod(sizeSeq(1:2)) sizeSeq(3)]);
% 
%         figure; plot(mean(seqAux));

        makeMovie(prepareMovieData(allSeq),fullfile(subDirectory,'global.avi'), false);
        %makeMovie(prepareMovieData(allSeq),fullfile(subDirectory,'global.mp4'), false);
        
%         blah = cell(1,11);
        
        % response transient per sequence
        for s = 1:16     
            seq = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,s,:,:)),[2 3 1]);
            
            if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
                %seq = (seq(trim+1:end-trim,trim+1:end-trim,:)-blankTrials)./blankTrials;
                seq = (seq(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:)-blankTrials)./blankTrials;
            else
                %seq = seq(trim+1:end-trim,trim+1:end-trim,:);
                seq = seq(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:);
            end
            
%             for i = 1:11
%                 blah{i}(s) = max(seq(:,:,i),[],'all');%-min(seq(:,:,i),[],'all');
%             end
%             blah(s) = mean(seq(:,:,6),'all');
            
            seq = prepareMovieData(seq);
            
            makeMovie(seq,fullfile(subDirectory,['seq' num2str(s) '.avi']),false);
            %makeMovie(seq,fullfile(subDirectory,['seq' num2str(s) '.mp4']),false);
        end
        
        % make movie of difference between global transient and no stimulus
        % transient
%         makeMovie(prepareMovieData(allSeq-blankTrials),fullfile(subDirectory,'difference.avi'),false);
       
%         for i = 1:11
%             figure; create_seq_eff_plot(normalize(blah{i}).',normalize(-six_hertz));
%         end

        %Repeat process, for behav
            %Note: Could interrupt mainline above, but choosing not to currently
        if isfield( R(fly).BLOCK(b) , 'dataSeqBehav' )
            disp(['-- Fly #',num2str(chosenFlies(fly)),', Block ',num2str(b),': Calculating transient movies for behav separated data --'])
            
            for statInd = 1:size( R(fly).BLOCK(b).dataSeqBehav,2 )
                thisBlockStateDirectory = [thisBlockDirectory,filesep,'State_',num2str(R(fly).BLOCK(b).dataSeqBehav(statInd).state),filesep,'Transients'];
                if ~exist(thisBlockStateDirectory,'dir')
                   mkdir(thisBlockStateDirectory); 
                end

                trim = trim; %Inherit
                blankTrials = blankTrials;
                
                thisData = R(fly).BLOCK(b).dataSeqBehav(statInd).meanDataSeqReduced;
                
                if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
                    % global response transient
                    allSeqBehav = permute(squeeze(mean( thisData ,2)),[2 3 1]);
                    %allSeqBehav = (allSeqBehav(trim+1:end-trim,trim+1:end-trim,:)-blankTrials)./blankTrials;
                    allSeqBehav = (allSeqBehav(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:)-blankTrials)./blankTrials;
                else
                    % global response transient
                    allSeqBehav = permute(squeeze(mean( thisData ,2)),[2 3 1]);
                    %allSeqBehav = allSeqBehav(trim+1:end-trim,trim+1:end-trim,:);
                    allSeqBehav = allSeqBehav(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:);
                end
                
                save(fullfile(thisBlockStateDirectory,'global_behav_transient'),'allSeqBehav');
                makeMovie(prepareMovieData(allSeqBehav),fullfile(thisBlockStateDirectory,'global_behav.avi'), false);
                
                for s = 1:16     
                    seqBehav = permute(squeeze(thisData(:,s,:,:)),[2 3 1]);
                    if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
                        %seqBehav = (seqBehav(trim+1:end-trim,trim+1:end-trim,:)-blankTrials)./blankTrials;
                        seqBehav = (seqBehav(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:)-blankTrials)./blankTrials;
                    else
                        %seqBehav = seqBehav(trim+1:end-trim,trim+1:end-trim,:);
                        seqBehav = seqBehav(trim(1)+1:end-trim(3),trim(4)+1:end-trim(2),:);
                    end
                    seqBehav = prepareMovieData(seqBehav);
                    makeMovie(seqBehav,fullfile(thisBlockStateDirectory,['seqBehav' num2str(s) '.avi']),false);
                    %makeMovie(seqBehav,fullfile(thisBlockStateDirectory,['seqBehav' num2str(s) '.mp4']),false);
                end
                disp(['State ',num2str(R(fly).BLOCK(b).dataSeqBehav(statInd).state),' sequence data/movies saved'])
                
            end
            
        end

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