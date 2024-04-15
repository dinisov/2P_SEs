close all; clear;

addpath('D:\group_vanswinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_vanswinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_vanswinderen\Dinis\2P SEs\Functions\');

close all; clear;

RDMDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\';

%where the sequence data is located (stimulus files)
sequenceDirectory = fullfile(RDMDirectory,'RPiData/');

%where the main data is found
dataDirectory = fullfile(RDMDirectory,'Gcamp7s_CC/');

outputDirectory = '../2P Results 2';

flyRecord = readtable('../2P Record/2P_record');

%get rid of excluded flies
flyRecord = flyRecord(~logical(flyRecord.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [32 32];

flyList = unique(flyRecord.Fly);

chosenFlies = 118;

% chosenFlies = flyList;%do not choose any flies

%whether to analyse grouped blocks
groupedBlocks = 0;

%% collate, reduce, filter and concatenate pre-aligned data

FLIES = collate2PData(flyRecord, chosenFlies, gridSize, dataDirectory, sequenceDirectory, groupedBlocks);

%% analyse SEs
% separates images according to preceding sequence of stimuli and
% calculates mean images as a function of the sequence

R = analyse2P(FLIES, chosenFlies, outputDirectory, groupedBlocks);

%% make movies of transients as differences to mean

transientMovies(R, chosenFlies, outputDirectory, 4);

%% calculate fit to SLRP, LRPR, SLRP+LRPR, and EPHYS (per volume/time and collapsed across time)

R = componentFits2P(R, FLIES, groupedBlocks);

%% make movies of fits over time

fitMovies(R, FLIES, outputDirectory, gridSize, chosenFlies, 3);

%% calculate mass t-tests

R = ttests2P(R, FLIES, groupedBlocks);

%% RRRR-RRRA and AAAA-AAAR (collapsed and over time videos)

patternPlots(R, FLIES, chosenFlies, outputDirectory);

%% L vs R analysis (t-tests, L-R, L and R, L and R movies)

analyseLvsR(R, FLIES, chosenFlies, outputDirectory, 3);

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

%% fit and plot some seq eff profiles of interest

% options = optimset('Algorithm','interior-point','FinDiffType','central');
%
% %10 Feb both
% % profiles = [23 19;28 14;23 5;23 19;26 4];
% 
% %10 Feb exp3
% % profiles = [13 18;23 28;21 12; 24 31; 26 3];
% 
% %10 Feb exp4
% profiles = [23 18; 28 17; 23 6; 28 14; 10 14];
% 
% type = [1 2 3 4 4 5];
% 
% names = {'slrp','lrpr','weird','combined','ephys'};
% 
% for i = 1:size(profiles,1)
%     
%     seq_eff_pattern = squeeze(meanVolData(profiles(i,1),profiles(i,2),:));
%    
%     switch type(i)
%         
%         case 1
%   
%            %fit only to slrp or lrpr
%            [x,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 0 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
%            seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%          
%         case 2
%             
%            [x,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 0 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
%            seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%             
%         case 3
%             
%            [x,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
%            seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%             
%         case 4
%             
%            %fit overall
%             [x,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf -inf 0 -inf],[inf inf 0 inf],[],options);
%             seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%             
%         case 5
%             
%             [x,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,six_hertz,seq_eff_pattern),[0 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
%             seq_eff_fit =  x(3)*six_hertz + x(4);
% 
%     end
%     
%     figure; create_seq_eff_plot(seq_eff_pattern,seq_eff_fit);
%     saveas(gcf,[ names{type(i)} '_best.png'])
%     
% end