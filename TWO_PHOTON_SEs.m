close all; clear;

addpath('C:\Users\uqmvan13\2p\Dinis Scripts\Global functions\');
addpath('C:\Users\uqmvan13\2p\Dinis Scripts\Indexes and legends\');
%addpath('D:\group_vanswinderen\Matt\2p\2P SEs\Functions\');
addpath('.\Functions\');
addpath('C:\Users\uqmvan13\2p\Matt Scripts\');

close all; clear;

%Matt/Dinish
%{{
RDMDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\';

%where the sequence data is located (stimulus files)
%sequenceDirectory = 'I:\RFDG2021-Q4413\2P_Data\RPiData';
sequenceDirectory = 'I:\RFDG2021-Q4413\Matt';


%where the main data is found
dataDirectory = fullfile(RDMDirectory,'Gcamp7s_CC/');

outputDirectory = '../2P_RESULTS_4';

%flyRecord = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");
recordPath = "I:\RFDG2021-Q4413\2P Record\2P_record.xlsx";
flyRecord = readtable(recordPath);
%}


%Andre params
%{
RDMDirectory = 'I:\RFDG2021-Q4413\Andre\2p_Data\\'; %Andre
sequenceDirectory = 'I:\RFDG2021-Q4413\Andre';
dataDirectory = fullfile(RDMDirectory);
outputDirectory = '../2P Results Andre';
%flyRecord = readtable("I:\RFDG2021-Q4413\Andre\2p_Record\Andre_2P_record.xlsx");
recordPath = "I:\RFDG2021-Q4413\Andre\2p_Record\Andre_2P_record.xlsx";
flyRecord = readtable(recordPath);
%}

%get rid of excluded flies
flyRecord = flyRecord(~logical(flyRecord.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [64 64];

% side trim
% trim = 3;

flyList = unique(flyRecord.Fly);

%chosenFlies = [354:357]; %Matt
chosenFlies = [384:387]; %Matt
chosenBlocks = {}; %Leave empty if not using
    %Note: If using, block/s must be specified for *all* chosen flies
chosenZ = {}; %Same format as chosenBlocks (Specified for all)

%chosenFlies = [119]; %Andre
%chosenBlocks = {[2]};
%chosenZ = {};

flyRecord = flyRecord(ismember(flyRecord.Fly,chosenFlies),:);

% flyRecord = flyRecord(1,:);% for testing

% chosenBlocks = 1;

%whether to analyse grouped blocks
groupedBlocks = 0;

%%
%Pre QA for bad block specification
if ~isempty(chosenBlocks)
    for fly = 1:length(chosenFlies)
        thisFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(fly),:);
        if any( ~ismember(chosenBlocks{fly},thisFlyBlocks.Block.') )
            ['## Alert: One or more requested blocks not existing in flyRecord ##']
            crash = yes
        end
    end
end

%% process flies

% transient movies; component fits; fit movies; t-tests; oddballs; LvsR; PCA; global transient
%analysisToggle = [1 1 1 0 1 0 1 1];
%analysisToggle = [1 0 0 0 0 0 1 0];
analysisToggle = [1 0 0 0 0 0 1 0];
separateByState = 0; %Whether to use available behav data to repeat processing on sleep vs wake, etc 
doRolling = 0; %Whether to also do rolling analysis
    %Note: Requires MATLAB >=2021

startTime = datetime('now');
%for fly = chosenFlies %Does not work with chosenBlocks/etc specification for some reason
for flyInd = 1:length(chosenFlies)
    fly = chosenFlies(flyInd)
    if ~isempty(flyRecord(flyRecord.Fly == fly,:))
        %processFlies(flyRecord, fly, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks);
        %processFlies(flyRecord, fly, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, separateByState);
        %processFlies(flyRecord, fly, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, separateByState, doRolling);
        %processFlies(flyRecord, fly, chosenBlocks, chosenZ, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, separateByState, doRolling);
        if ~isempty(chosenBlocks)
            theseChosenBlocks = chosenBlocks(flyInd); %Ameliorates issues with irregular block requests
        else
            theseChosenBlocks = {};
        end
        if ~isempty(chosenZ)
            theseChosenZs = chosenZ(flyInd); %Ditto
        else
            theseChosenZs = {};
        end
        %processFlies(flyRecord, fly, theseChosenBlocks, theseChosenZs, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, separateByState, doRolling);
        %fake
        processFlies(flyRecord, fly, theseChosenBlocks, theseChosenZs, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, separateByState, doRolling, recordPath);
    end
end
endTime = datetime('now');
MET = seconds(endTime-startTime);
disp([char(10),'---------------------------------',char(10),'Total analysis time: ',num2str(MET/60),'m'])
