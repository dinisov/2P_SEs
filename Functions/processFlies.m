function processFlies(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% collate, reduce, filter and concatenate pre-aligned data

FLIES = collate2PData(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, groupedBlocks);

%% analyse SEs
% separates images according to preceding sequence of stimuli and
% calculates mean images as a function of the sequence

R = analyse2P(FLIES, chosenFlies, outputDirectory, groupedBlocks);

%% make movies of transients as differences to mean
if analysisToggle(1)
    transientMovies(R, chosenFlies, outputDirectory);
end

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS (per volume/time and collapsed across time)
if analysisToggle(2)
    R = componentFits2P(R, groupedBlocks);
end
%% make movies of fits over time
if analysisToggle(3)
    fitMovies(R, outputDirectory, gridSize, chosenFlies);
end
%% calculate mass t-tests
if analysisToggle(4)
    R = ttests2P(R, groupedBlocks);
end
%% RRRR-RRRA and AAAA-AAAR (collapsed and over time videos)
if analysisToggle(5)
    patternPlots(R, chosenFlies, outputDirectory);
end
%% L vs R analysis (t-tests, L-R, L and R, L and R movies)
if analysisToggle(6)
    analyseLvsR(R, chosenFlies, outputDirectory);
end

%% PCA analysis
if analysisToggle(7)
    flyPCA2D(R, chosenFlies, outputDirectory);
end

%% global transient analysis


%% plotting
if analysisToggle(7)
    disp('Plotting stuff');
    tic;
    % plot results per block
    for fly = 1:length(FLIES)
        subDirectory = fullfile(outputDirectory,['Fly' num2str(chosenFlies(fly))]);
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        plotFly(R(fly), groupedBlocks, subDirectory,'off');
    end
    toc;
end

end

