function fitMovies(R, outputDirectory, gridSize, chosenFlies)
% fitMovies Makes movies of component fits over time
%   Detailed explanation goes here

    disp('Making fit movies...');

    tic;

    for fly = 1:length(R)
        for b = [R(fly).BLOCK.blockNum]
            subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'Movies');
            if ~exist(subDirectory,'dir')
               mkdir(subDirectory); 
            end
            
            trimSize = R(fly).BLOCK(b).Trim;
            newGridSize = gridSize-2*trimSize;

            %alternation component movie
            movieMatrix = zeros([newGridSize R(fly).BLOCK(b).nVol]);
            for vol = 1:R(fly).BLOCK(b).nVol
                movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_slrp(trimSize+1:end-trimSize,trimSize+1:end-trimSize);
            end
            makeMovie(movieMatrix,fullfile(subDirectory,'ALT_fit.avi'),true);

            %repetition component movie
            movieMatrix = zeros([newGridSize R(fly).BLOCK(b).nVol]);
            for vol = 1:R(fly).BLOCK(b).nVol
                movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_lrpr(trimSize+1:end-trimSize,trimSize+1:end-trimSize);
            end
            makeMovie(movieMatrix,fullfile(subDirectory,'REP_fit.avi'),true);

            %weird component movie
            movieMatrix = zeros([newGridSize R(fly).BLOCK(b).nVol]);
            for vol = 1:R(fly).BLOCK(b).nVol
                movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_weird(trimSize+1:end-trimSize,trimSize+1:end-trimSize);
            end
            makeMovie(movieMatrix,fullfile(subDirectory,'WEIRD_fit.avi'),true);
        end
    end
    toc;
    
end

