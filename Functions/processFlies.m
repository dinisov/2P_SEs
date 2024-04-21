function processFlies(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, outputDirectory, groupedBlocks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% collate, reduce, filter and concatenate pre-aligned data

FLIES = collate2PData(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, groupedBlocks);

%% analyse SEs
% separates images according to preceding sequence of stimuli and
% calculates mean images as a function of the sequence

R = analyse2P(FLIES, chosenFlies, outputDirectory, groupedBlocks);

%% make movies of transients as differences to mean

transientMovies(R, chosenFlies, outputDirectory);

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS (per volume/time and collapsed across time)

R = componentFits2P(R, groupedBlocks);

%% make movies of fits over time

fitMovies(R, outputDirectory, gridSize, chosenFlies);

%% calculate mass t-tests

R = ttests2P(R, groupedBlocks);

%% RRRR-RRRA and AAAA-AAAR (collapsed and over time videos)

patternPlots(R, chosenFlies, outputDirectory);

%% L vs R analysis (t-tests, L-R, L and R, L and R movies)

analyseLvsR(R, chosenFlies, outputDirectory);

%% plotting
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

