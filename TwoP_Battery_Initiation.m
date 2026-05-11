%Fly processing non-wrapper for Andre
%Based on a hybridisation of TWO_PHOTON_SEs and processFlies

close all; clear;

% addpath('D:\group_vanswinderen\Dinis\Scripts\Global functions\');
% addpath('D:\group_vanswinderen\Dinis\Scripts\Indexes and legends\');
% %addpath('D:\group_vanswinderen\Matt\2p\2P SEs\Functions\');
% addpath('.\Functions\');
addpath("C:\Users\uqbshre2\2p\2P SEs\Functions");
addpath("C:\Users\uqbshre2\2p\dinis scripts\Global functions");
addpath('C:\Users\uqbshre2\2p\Matt Scripts\');
% addpath('C:\Users\uqaspen5\2p\2p\dinis scripts\Global functions');

%%

RDMDirectory = 'I:\bhanu2026-Q9370\2P\2P_Data\'; %Bhanu
sequenceDirectory = 'I:\bhanu2026-Q9370\2P';
dataDirectory = fullfile(RDMDirectory);
outputDirectory = '../2P_Results';
flyRecord = readtable("I:\bhanu2026-Q9370\2P\2p_record\Bhanu_2P_record.xlsx");

% RDMDirectory = 'I:\RFDG2021-Q4413\Andre\2p_Data\\'; %Andre
% sequenceDirectory = 'I:\RFDG2021-Q4413\Andre';
% dataDirectory = fullfile(RDMDirectory);
% outputDirectory = '../2P Results Andre';
% flyRecord = readtable("I:\RFDG2021-Q4413\Andre\2p_Record\Andre_2P_record.xlsx");

% RDMDirectory = 'I:\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\'; %Matt
% sequenceDirectory = 'I:\RFDG2021-Q4413\Andre';
% dataDirectory = fullfile(RDMDirectory);
% outputDirectory = '../2P Results Andre';
% flyRecord = readtable("I:\RFDG2021-Q4413\2P Record\2P_record.xlsx");

%%

flyRecord = flyRecord(~logical(flyRecord.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [64 64];

% side trim
% trim = 3;

flyList = unique(flyRecord.Fly);

chosenFlies = 13; 
chosenBlocks = {1}; 
chosenZ = {};


%QA for(and from) Andre
if size(chosenFlies, 2) ~= size(chosenBlocks, 2)
    disp("blocks don't match flies")
    crash = yes;
end

flyRecord = flyRecord(ismember(flyRecord.Fly,chosenFlies),:);

groupedBlocks = 0;
analysisToggle = [1 0 0 0 0 0 1 0];
separateByState = 0; %Whether to use available behav data to repeat processing on sleep vs wake, etc 
doRolling = 0;

%%

fly = chosenFlies


%%

disp([char(10),'-------------------------------------'])

%FLIES = collate2PData(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, groupedBlocks, separateByState,doRolling); 
%FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory, groupedBlocks, separateByState,doRolling, 0);
FLIES = collate2PData(flyRecord, chosenFlies, chosenBlocks, gridSize, dataDirectory, sequenceDirectory, 'alternateUseCase', 0, 'reqZ', chosenZ);

close all
%
%fake
% [FLIES] = syncMaster_legacy( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', 0, 'rollingAnalysis', -1, 'disregardRollingDesign', 0, 'disregardBattery', 0 );
[FLIES] = syncMaster( FLIES , flyRecord, 'dataDirectory', dataDirectory, 'doPlot', 0, 'doVid', 0, 'rollingAnalysis', -1,...
    'disregardRollingDesign', 0, 'disregardBattery', 0, 'overwriteShortcut', 0, 'unsiphonedSEs', 1 );


%
% Your image analysis script here
%