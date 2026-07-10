close all; clear;

addpath('..\Dinis Scripts\Global functions\');
addpath('..\Dinis Scripts\Indexes and legends\');
%addpath('D:\group_vanswinderen\Matt\2p\2P SEs\Functions\');
addpath('.\Functions\');
addpath('..\Matt Scripts\');

close all; clear;

%Matt/Dinish
%{{
RDMDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\';
%where the sequence data is located (stimulus files)
%sequenceDirectory = 'I:\RFDG2021-Q4413\2P_Data\RPiData';
sequenceDirectory = 'I:\RFDG2021-Q4413\Matt';
%where the main data is found
dataDirectory = fullfile(RDMDirectory,'Gcamp7s_CC/');
outputDirectory = '../2P_RESULTS_5';
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
recordPath = []; %Overwrite empty to prevent updating
%}

%Bhanu
%{
RDMDirectory = 'I:\bhanu2026-Q9370\2P\2P_Data\';
%where the sequence data is located (stimulus files)
sequenceDirectory = 'I:\bhanu2026-Q9370\2P';
%where the main data is found
dataDirectory = fullfile(RDMDirectory);
outputDirectory = '../2P_Results';
recordPath = "I:\bhanu2026-Q9370\2P\2p_record\2P_record_Bhanu_New.xlsx";
flyRecord = readtable(recordPath);
recordPath = [];
%}

%get rid of excluded flies
flyRecord = flyRecord(~logical(flyRecord.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [64 64];

% side trim
% trim = 3;

flyList = unique(flyRecord.Fly);

chosenFlies = [396,399,343,344,327,328,312,314,316,319,322,327,328]; %Matt
%chosenFlies = [318,319,321,322,323,343,338,344,399,396]; %Matt
chosenBlocks = {}; %Leave empty if not using
    %Note: If using, block/s must be specified for *all* chosen flies
chosenZ = {}; %Same format as chosenBlocks (Specified for all)
%[chosenFlies, chosenBlocks] = flyProvider( "I:\RFDG2021-Q4413\2P Record\2P_record_datasets.xlsx", 'Whisky', ...
%    'expandedArchitecture', 0 );

%chosenFlies = [15]; %Bhanu
%chosenBlocks = {[2]};
%chosenZ = {};

flyRecord = flyRecord(ismember(flyRecord.Fly,chosenFlies),:);

% flyRecord = flyRecord(1,:);% for testing

% chosenBlocks = 1;

%whether to analyse grouped blocks
groupedBlocks = 0;

%Other, mostly syncMaster-specific options
unsiphonedSEs = 0;
overwriteShortcut = 0;
separateByState = 0; %Whether to use available behav data to repeat processing on sleep vs wake, etc 
doRolling = 0; %Whether to also do rolling analysis (DEPRECATED)
    %Note: Requires MATLAB >=2021
    %Secondary note: Not tested for a very long time

%%
%Pre QA for bad block specification
if ~isempty(chosenBlocks)
    for fly = 1:length(chosenFlies)
        thisFlyBlocks = flyRecord(flyRecord.Fly == chosenFlies(fly),:);
        if any( ~ismember(chosenBlocks{fly},thisFlyBlocks.Block.') )
            ['## Alert: One or more requested blocks not existing in flyRecord ##']
            ['chosenFly ',num2str(fly),' - #',num2str(chosenFlies(fly))]
            ['requested blocks: ',num2str(chosenBlocks{fly})]
            ['available blocks: ',num2str(thisFlyBlocks.Block.')]
            crash = yes
        end
    end
end

%% process flies

% transient movies; component fits; fit movies; t-tests; oddballs; LvsR; PCA; global transient
%analysisToggle = [1 1 1 0 1 0 1 1];
%analysisToggle = [1 0 0 0 0 0 1 0];
analysisToggle = [1 0 0 0 0 0 1 0];

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
        %processFlies(flyRecord, fly, theseChosenBlocks, theseChosenZs, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, separateByState, doRolling, recordPath);
        processFlies(flyRecord, fly, theseChosenBlocks, theseChosenZs, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, ...
            'separateByState', separateByState, 'doRolling', doRolling, 'recordPath',recordPath,...
            'syncManUnsiphonedSEs',unsiphonedSEs,'syncManOverwriteShortcut',overwriteShortcut,...
            'syncMan2ndStageSmooth',1);
    end
end
endTime = datetime('now');
MET = seconds(endTime-startTime);
disp([char(10),'---------------------------------',char(10),'Total analysis time: ',num2str(MET/60),'m'])
close all
