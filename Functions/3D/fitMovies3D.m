function fitMovies3D(R, FLIES, outputDirectory, gridSize, chosenFlies, trimSize)
% fitMovies Makes movies of component fits over time
%   Detailed explanation goes here

    disp('Making fit movies...');

    tic;
    newGridSize = gridSize-2*trimSize;

    for fly = 1:length(FLIES)
        
        for b = 1:length(FLIES(fly).BLOCKS)
            subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))],['Block' num2str(b)],'Movies');
            if ~exist(subDirectory,'dir')
               mkdir(subDirectory); 
            end
            
            n_z = size(R(fly).BLOCK(b).meanDataSeq,5);
            
            %TODO: make for loop of movies below

            %alternation component movie
            movieMatrix = zeros([newGridSize n_z FLIES(fly).BLOCKS(b).nVol]);
            for vol = 1:FLIES(fly).BLOCKS(b).nVol
                movieMatrix(:,:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_slrp(trimSize+1:end-trimSize,trimSize+1:end-trimSize,:);
            end
            makeMovie3D(movieMatrix,fullfile(subDirectory,'ALT_fit.avi'),true);

            %repetition component movie
            movieMatrix = zeros([newGridSize n_z FLIES(fly).BLOCKS(b).nVol]);
            for vol = 1:FLIES(fly).BLOCKS(b).nVol
                movieMatrix(:,:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_lrpr(trimSize+1:end-trimSize,trimSize+1:end-trimSize,:);
            end
            makeMovie3D(movieMatrix,fullfile(subDirectory,'REP_fit.avi'),true);

            %weird component movie
            movieMatrix = zeros([newGridSize n_z FLIES(fly).BLOCKS(b).nVol]);
            for vol = 1:FLIES(fly).BLOCKS(b).nVol
                movieMatrix(:,:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_weird(trimSize+1:end-trimSize,trimSize+1:end-trimSize,:);
            end
            makeMovie3D(movieMatrix,fullfile(subDirectory,'WEIRD_fit.avi'),true);
        end
    end
    toc;
    
end

