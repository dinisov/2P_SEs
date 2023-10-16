function transientMovies(R, chosenFlies, resultsDirectory, trim)
%transientMovies Summary of this function goes here
%   Detailed explanation goes here

for fly = 1:length(R)
    
    for b = 1:length(R(fly).BLOCK)
        
        thisBlockDirectory = fullfile(resultsDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)]);
        
        subDirectory = fullfile(thisBlockDirectory,'Transients');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        
        % global response transient
        allSeq = permute(squeeze(mean(R(fly).BLOCK(b).meanDataSeq,2)),[2 3 1]);
        allSeq = allSeq(trim+1:end-trim,trim+1:end-trim,:);
        allSeq = prepareMovieData(allSeq);
        
        makeMovie(allSeq,fullfile(subDirectory,'global.avi'), false);
        
        % response transient per sequence
        for s = 1:16     
            seq = permute(squeeze(R(fly).BLOCK(b).meanDataSeq(:,s,:,:)),[2 3 1]);
            seq = seq(trim+1:end-trim,trim+1:end-trim,:);
            seq = prepareMovieData(seq);
            
            makeMovie(seq,fullfile(subDirectory,['seq' num2str(s) '.avi']),false);
        end
        
    end
        
end

end

function data = prepareMovieData(data)

    data = data(:,:,2:end); %get rid of first time point
    data = data-mean(data,3); % normalise by mean along time
%     data = data-data(:,:,end); % normalise by time point

    %map to interval [0,1]
    data = data+abs(min(data,[],'all')); % make minimum 0
    data = data/max(data,[],'all'); %make maximum 1
    
end