function transientMovies3D(R, chosenFlies, resultsDirectory)
%transientMovies Summary of this function goes here
%   Detailed explanation goes here

% load six_hertz.mat

for fly = 1:length(R)
    
    for b = [R(fly).BLOCK.blockNum]%1:length(R(fly).BLOCK)
        
        thisBlockDirectory = fullfile(resultsDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]);
        
        subDirectory = fullfile(thisBlockDirectory,'Transients');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end

        trim = R(fly).BLOCK(b).Trim;
        
        % blank trial transient movie
        if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
            
            blankTrials = R(fly).BLOCK(b).meanBlankTransient;
            blankTrials = blankTrials(trim+1:end-trim,trim+1:end-trim,:,:);
            
            makeMovie3D(prepareMovieData(blankTrials),fullfile(subDirectory,'blankTrials.avi'),false);
            
            % global response transient
            allSeq = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,2)),[2 3 4 1]);
            allSeq = (allSeq(trim+1:end-trim,trim+1:end-trim,:,:) - blankTrials)./blankTrials;
        else
            
            meanTransient = R(fly).BLOCK(b).meanTransient;
            meanTransient = meanTransient(trim+1:end-trim,trim+1:end-trim,:,:);
            
            % global response transient
            allSeq = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,2)),[2 3 4 1]);
            allSeq = (allSeq(trim+1:end-trim,trim+1:end-trim,:,:) - meanTransient)./meanTransient;
        end
        
        save(fullfile(subDirectory,'global_transient'),'allSeq');
        
%         sizeSeq = size(allSeq);
% 
%         seqAux = reshape(allSeq,[prod(sizeSeq(1:2)) sizeSeq(3)]);
% 
%         figure; plot(mean(seqAux));

        makeMovie3D(prepareMovieData(allSeq),fullfile(subDirectory,'global.avi'),false);
        
%         blah = cell(1,11);
        
        % response transient per sequence
        sepSeq = [];
        for s = 1:16     
            seq = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,s,:,:,:)),[2 3 4 1]);
            
            if ~isempty(R(fly).BLOCK(b).meanBlankTransient)
                seq = (seq(trim+1:end-trim,trim+1:end-trim,:,:) - blankTrials)./blankTrials;
            else
                seq = (seq(trim+1:end-trim,trim+1:end-trim,:,:) - meanTransient)./meanTransient;
            end
            
%             for i = 1:11
%                 blah{i}(s) = max(seq(:,:,i),[],'all');%-min(seq(:,:,i),[],'all');
%             end
%             blah(s) = mean(seq(:,:,6),'all');
            sepSeq(:,:,:,:,s) = seq;
            
            seq = prepareMovieData(seq);
            
            makeMovie3D(seq,fullfile(subDirectory,['seq' num2str(s) '.avi']),false);
        end
        
        save(fullfile(subDirectory,'separated_transients'),'sepSeq');
        
        % make movie of difference between global transient and no stimulus
        % transient
%         makeMovie(prepareMovieData(allSeq-blankTrials),fullfile(subDirectory,'difference.avi'),false);
       
%         for i = 1:11
%             figure; create_seq_eff_plot(normalize(blah{i}).',normalize(-six_hertz));
%         end

    end
        
end

end

function data = prepareMovieData(data)

    %map to interval [0,1]
    data = data - min(data,[],'all'); % make minimum 0
    data = data / max(data,[],'all'); %make maximum 1
    
end