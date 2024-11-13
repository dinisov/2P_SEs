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
            
            %trimSize = R(fly).BLOCK(b).Trim;
            %newGridSize = gridSize-2*trimSize;
            if ~isnan(R(fly).BLOCK(b).Trim) && ( ~isfield(R(fly).BLOCK(b),'TrimCoords') || isempty(R(fly).BLOCK(b).TrimCoords) )
                trimSize = repmat( R(fly).BLOCK(b).Trim , 1 , 4 ); %Make coord-like, to simplify later
                disp(['Using unitary trim for fit movie creation'])
            else
                trimSize = R(fly).BLOCK(b).TrimCoords; %Note matrix, not singular
                disp(['Using coordinate trim for fit movie creation'])
            end
            
            newGridSize = [gridSize(1) - (trimSize(1)+trimSize(3)), gridSize(2) - (trimSize(2)+trimSize(4)) ];
                %Under unitary conditions, this should converge to the same
                %resulting size as before

            %alternation component movie
            movieMatrix = zeros([newGridSize R(fly).BLOCK(b).nVol]);
            for vol = 1:R(fly).BLOCK(b).nVol
                %movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_slrp(trimSize+1:end-trimSize,trimSize+1:end-trimSize);
                movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_slrp(trimSize(1)+1:end-trimSize(3),trimSize(4)+1:end-trimSize(2));
            end
            makeMovie(movieMatrix,fullfile(subDirectory,'ALT_fit.avi'),true);
            %makeMovie(movieMatrix,fullfile(subDirectory,'ALT_fit.mp4'),true);

            %repetition component movie
            movieMatrix = zeros([newGridSize R(fly).BLOCK(b).nVol]);
            for vol = 1:R(fly).BLOCK(b).nVol
                %movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_lrpr(trimSize+1:end-trimSize,trimSize+1:end-trimSize);
                movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_lrpr(trimSize(1)+1:end-trimSize(3),trimSize(4)+1:end-trimSize(2));
            end
            makeMovie(movieMatrix,fullfile(subDirectory,'REP_fit.avi'),true);
            %makeMovie(movieMatrix,fullfile(subDirectory,'REP_fit.mp4'),true);

            %weird component movie
            movieMatrix = zeros([newGridSize R(fly).BLOCK(b).nVol]);
            for vol = 1:R(fly).BLOCK(b).nVol
                %movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_weird(trimSize+1:end-trimSize,trimSize+1:end-trimSize);
                movieMatrix(:,:,vol) = R(fly).BLOCK(b).rVol(vol).r.r_weird(trimSize(1)+1:end-trimSize(3),trimSize(4)+1:end-trimSize(2));
            end
            makeMovie(movieMatrix,fullfile(subDirectory,'WEIRD_fit.avi'),true);
            %makeMovie(movieMatrix,fullfile(subDirectory,'WEIRD_fit.mp4'),true);
        end
    end
    toc;
    
end

