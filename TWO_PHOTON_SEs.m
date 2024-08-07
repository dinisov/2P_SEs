close all; clear;

addpath('D:\group_vanswinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_vanswinderen\Dinis\Scripts\Indexes and legends\');
addpath('D:\group_vanswinderen\Matt\2p\2P SEs\Functions\');

close all; clear;

RDMDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\';

%where the sequence data is located (stimulus files)
%sequenceDirectory = 'I:\RFDG2021-Q4413\2P_Data\RPiData';
sequenceDirectory = 'I:\RFDG2021-Q4413\Matt';

%where the main data is found
dataDirectory = fullfile(RDMDirectory,'Gcamp7s_CC/');

outputDirectory = '../2P Results 3';

flyRecord = readtable("D:\group_vanswinderen\Dinis\2P Record\2P_record");

%get rid of excluded flies
flyRecord = flyRecord(~logical(flyRecord.Exclude),:);

imageSize = [128 128];

% final grid size
gridSize = [64 64];

% side trim
% trim = 3;

flyList = unique(flyRecord.Fly);

chosenFlies = [185];%174,175,176,177,178,179,180,181];

flyRecord = flyRecord(ismember(flyRecord.Fly,chosenFlies),:);

% flyRecord = flyRecord(1,:);% for testing

% chosenBlocks = 1;

%whether to analyse grouped blocks
groupedBlocks = 0;

%% process flies

% transient movies; component fits; fit movies; t-tests; oddballs; LvsR;
% PCA; global transient
analysisToggle = [1 1 1 0 1 0 1 1];
separateByState = 1; %Whether to use available behav data to repeat processing on sleep vs wake, etc 

for fly = chosenFlies
    if ~isempty(flyRecord(flyRecord.Fly == fly,:))
        %processFlies(flyRecord, fly, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks);
        processFlies(flyRecord, fly, gridSize, dataDirectory, sequenceDirectory, outputDirectory, analysisToggle, groupedBlocks, separateByState);
    end
end
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