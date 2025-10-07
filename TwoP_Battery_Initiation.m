%Fly processing non-wrapper for Andre
%Based on a hybridisation of TWO_PHOTON_SEs and processFlies

close all; clear;

% addpath('D:\group_vanswinderen\Dinis\Scripts\Global functions\');
% addpath('D:\group_vanswinderen\Dinis\Scripts\Indexes and legends\');
% %addpath('D:\group_vanswinderen\Matt\2p\2P SEs\Functions\');
addpath('.\Functions\');
addpath('C:\Users\uqmvan13\2p\Dinis Scripts\Global functions\');

%%

%RDMDirectory = 'I:\RFDG2021-Q4413\Andre\2p_Data\\'; %Andre
%sequenceDirectory = 'I:\RFDG2021-Q4413\Andre';
%dataDirectory = fullfile(RDMDirectory);
%outputDirectory = '../2P Results Andre';
%flyRecord = readtable("I:\RFDG2021-Q4413\Andre\2p_Record\Andre_2P_record.xlsx");

RDMDirectory = 'I:\RFDG2021-Q4413\2P_Data\'; %Matt
sequenceDirectory = 'I:\RFDG2021-Q4413\Matt';
dataDirectory = fullfile(RDMDirectory,'Gcamp7s_CC/');
outputDirectory = '../2P_RESULTS_4';
recordPath = "I:\RFDG2021-Q4413\2P Record\2P_record.xlsx";
flyRecord = readtable(recordPath);

%%

flyRecord = flyRecord(~logical(flyRecord.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [64 64];

% side trim
% trim = 3;

flyList = unique(flyRecord.Fly);

chosenFlies = [379,380,382,383]; %[96:103]; %[85, 86, 87, 89, 91, 104, 106:110, 138:140]; %[112:118]; %[119, 120, 121, 122, 123, 124, 125, 131, 132]; %[85, 86, 87, 89, 91, 92, 93, 94]; %[85, 86, 87, 89, 91, 92, 93, 94, 104, 106:110]; %[85, 86, 87, 89, 91, 92, 93, 94]; %[83, 85, 86, 88, 89, 91, 92]; %was 37, 76 [+85, 87]
chosenBlocks = {}; %{[1:3], [1:3], [1:3], [1:3], [1:3], [2:4], [1:3], [1:3]}; %repmat({1}, 1, size(chosenFlies, 2)); %{[2], [2], [2], [2], [2], [2], [2]}; %{[1],[1],[1],[1],[2],[1],[1],[1],[1]}; %{[2:4], [2:4], [3,4], [2:4], [3:5], [2,3], [2,3], [3,4], [2:4]}; %{[1], [1], [1], [1], [1], [1], [1], [1]}; %{[1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1]};
chosenZ = {};

%QA for(and from) Andre
if ~isempty(chosenBlocks) && size(chosenFlies, 2) ~= size(chosenBlocks, 2)
    disp("blocks don't match flies")
    crash = yes;
end

flyRecord = flyRecord(ismember(flyRecord.Fly,chosenFlies),:);

groupedBlocks = 0;
%analysisToggle = [1 0 0 0 0 0 1 0];
separateByState = 0; %Whether to use available behav data to repeat processing on sleep vs wake, etc 
doRolling = 0;

%%

fly = chosenFlies


%%

disp([char(10),'-------------------------------------'])

%FLIES = collate2PData(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, groupedBlocks, separateByState,doRolling); 
%FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory, groupedBlocks, separateByState,doRolling, 0);
%FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory, 'alternateUseCase', 0, 'reqZ', chosenZ);
FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory,...
    'alternateUseCase', 0, 'reqZ', chosenZ,...
    'dynamicSG',1,'dynamicSGTimeWidth',10);

close all
%
%fake
% [FLIES] = syncMaster_legacy( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', 0, 'rollingAnalysis', -1, 'disregardRollingDesign', 0, 'disregardBattery', 0 );
[FLIES] = syncMaster( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', 0, 'rollingAnalysis', -1,...
    'disregardRollingDesign', 0, 'disregardBattery', 0, 'overwriteShortcut', 0, 'unsiphonedSEs', 0 );


%
% Your image analysis script here
%