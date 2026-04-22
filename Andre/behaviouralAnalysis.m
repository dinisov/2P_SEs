%%% behavioural analysis script
%uses BT data to analyse sleep characteristics/general fly behaviour
%will eventually integrate pixel subtraction results for comparison
clear
close all

%reading fly record
blocks = readtable("I:\RFDG2021-Q4413\Andre\2p_Record\Andre_2P_record.xlsx");
blocks = blocks(~logical(blocks.Exclude),:);
noDataList = [];
nonBatteryList = [];

%directories
dataDirectory = 'I:\RFDG2021-Q4413\Andre\2p_Data';
videoDirectory = 'I:\RFDG2021-Q4413\Andre\VIdeos';

%options
doDots = 1;
activityThreshold = 0.5; %used to threshold dWLE data according to active vs inactive

% chosenFlies = [119]; %[126:128, 141:145];
% chosenBlocks = {};
% %bat ATRs
% chosenFlies = [96:103];
% chosenBlocks = {};
% 
%bat ctrls
% chosenFlies = [104, 106:111, 153];
% 
% chosenBlocks = {};

% % SE ctrls
% chosenFlies = [129, 130, 133:137, 146:147];
% 
% chosenBlocks = repmat({[2:4]}, 1, size(chosenFlies, 2));

% chosenBlocks = {[2:4], [2:4], [2:4], [2:4], [2:4], [2:4], [2:4], [2:4], [2:4], [2:4]};

%SE ATRs
chosenFlies = [119:125, 131, 132];

chosenBlocks = {[2:4], [2:4], [2:4], [2:4], [3:5], [2:4], [2:4], [2:4], [2:4]};
% chosenBlocks = repmat({[2:4]}, 1, size(chosenFlies, 2));

activityStruct = struct;

for fly = 1:length(chosenFlies)

    disp(['commencing analysis for fly ', num2str(chosenFlies(fly))])

    thisFly = chosenFlies(fly);
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);

    if ~isempty(chosenBlocks) && ~isempty(chosenBlocks{fly})
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
    end

    nBlocks = height(thisFlyBlocks);

    if nBlocks == 0
        ['-# Alert: No apparent blocks (for fly ',num2str(chosenFlies(fly)),') #-']
        %crash = yes
        noDataList = [noDataList,thisFly];
        continue
    end

    % x titles used later in cross-block analysis
    xTitles = cell(nBlocks);
    blockTitle = 'Block ';

    finalBTTime = 0;

    for b = 1:nBlocks
        % stuff for cross block/fly analyses
        xTitles{b} = strcat(blockTitle, num2str(b));

        %loading fly data
        currentBlock = thisFlyBlocks(b,:);
        currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
        %         flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        currentDate_2 = char(datetime(thisFlyBlocks.Date(1),'Format','ddMMyy')); %2nd date uses different format
        flyID = [currentDate_2, '_f', num2str(currentBlock.FlyOnDay), '_b', num2str(currentBlock.Block)];

        flyVidID = ['fly', num2str(currentBlock.FlyOnDay), '_', currentDate_2([1,2]), '_', currentDate_2([3,4]), '_', currentDate_2([5,6])]; %, '_', '*', '_mov.csv']; %alternate fly name format need for videos


        btDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btData.csv'];
        btAncillaryDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_btAncillary.mat'];
        matDataDirectory = [dataDirectory, filesep, currentDate, filesep, 'MAT', filesep, flyID, '_MAT.mat'];
        pixelSubDirectory = [videoDirectory, filesep, currentDate, filesep, flyVidID, '_', '*', '_mov.csv'];

        BTData = csvread(btDataDirectory);
        load(btAncillaryDirectory);
        load(matDataDirectory);
        % 
        % videoFiles = dir(pixelSubDirectory);
        % videoData = nan(size(videoFiles, 1)*46000, ); %46000 because largest (from few I've investigated) was 45971
        % for vidNum = 1:size(videoFiles, 1)
        %     currentVideoData = csvread([videoFiles(vidNum).folder, videoFiles(vidNum).name]);
        % 
        % end
        % videoData = csvread(pixelSubDirectory);

        %using dWLE for avg activity (col 3 BT)
        dWLEraw = [BTData(:,6) + finalBTTime, abs(BTData(:,3))];
        finalBTTime = finalBTTime + BTData(end, 6);

        avg_dWLE = nanmean(abs(BTData(:,3)));
        % avg_dWLE = nanmean(BTData(:,3));
        SEMdWLE = nanstd(abs(BTData(:,3)), [], 1)/sqrt(size(BTData(:,3), 1));

        %thresholded dWLE
        dWLEactive = dWLEraw;
        dWLEactive(find(dWLEactive(:,2)<activityThreshold), 2) = 0;
        dWLEactive(find(dWLEactive(:,2)>activityThreshold), 2) = 1;


        activityStruct(fly).blocks(b).dWLEraw = dWLEraw;
        activityStruct(fly).blocks(b).dWLE = avg_dWLE;
        activityStruct(fly).blocks(b).dWLEerror = SEMdWLE;

        activityStruct(fly).blocks(b).activeMoments = dWLEactive;

        %using rawBT for avg activity
    end

    %% cross block analysis
    flyBlocks = activityStruct(fly).blocks;

    %raw dWLE
    figure
    hold on
    for block = 1:size(flyBlocks, 2)
        plot(flyBlocks(block).dWLEraw(:,1), flyBlocks(block).dWLEraw(:,2))
    end
    set(gcf, 'Name', ['fly ', num2str(chosenFlies(fly)), 'dWLE'])
    ylabel('absolute rotational velocity')
    xlabel('time (s)')

    %activity Measure
    figure
    hold on
    for block = 1:size(flyBlocks, 2)
        plot(flyBlocks(block).activeMoments(:,1), flyBlocks(block).activeMoments(:,2))
    end
    set(gcf, 'Name', ['fly ', num2str(chosenFlies(fly)), 'active moments - from dWLE'])
    ylabel('absolute rotational velocity')
    xlabel('time (s)')

    %dWLE
    figure
    hold on
    dWLEfig = bar([1:nBlocks], [flyBlocks(:).dWLE], 'FaceColor', 'flat', 'LineWidth', 2);
    % dWLEfig.CData
    errorbar([1:nBlocks], [flyBlocks(:).dWLE], [flyBlocks(:).dWLEerror], '.', 'Color', 'k')
    xticks([1:nBlocks])
    xticklabels(xTitles)
    ylabel('dWLE') %change to something more informative later
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'box', 'off')

    set(gcf, 'Name', ['fly ', num2str(chosenFlies(fly)), ' dWLE'], 'NumberTitle', 'off')
end

%% cross block/multi fly analysis
if size(chosenFlies, 2) > 1
    % %plotting dWLE across time
    % %stitching BTData together
    % crossFlyBTData = nan(size())
    % for fly = 1:size(chosenFlies, 2)
    %
    % end

    %dWLE avg
    disp('currently analyses first 3 blocks exclusively')
    % blockNums = [ActivityStruct]
    allFlydWLE = nan(3, size(chosenFlies, 2));
    for fly = 1:size(chosenFlies, 2)
        for block = 1:size(activityStruct(fly).blocks, 2)
            allFlydWLE(block, fly) = activityStruct(fly).blocks(block).dWLE;
        end
    end
    allFlydWLEerror = nanstd(allFlydWLE, [], 2)/sqrt(size(allFlydWLE, 2));

     %stats
     dWLEData = allFlydWLE';

     nonNormCheck = 0;
     for col = 1:size(dWLEData, 2)
         normTest = lillietest(dWLEData(:,col));
         if ~normTest
             nonNormCheck = 1;
         end
     end

     if nonNormCheck
         disp(['a group was not normally distributed, Friedman test being used for stats'])
         % [~,~,LRStats] = kruskalwallis(LRData, [] , 'off');
         [~,~,dWLEStats] = friedman(dWLEData, 1 , 'off');
         % close, close %closing obligatory figures that anova makes
         dWLESig = multcompare(dWLEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
         % LRDescription = ['kruskal-wallis test', newline];
         behavDescription = ['Friedman test', newline];
     else
         % [~,~,LRStats] = anova1(LRData, [] , 'off');
         [~,~,dWLEStats] = ranova(dWLEData, [] , 'off');
         % close, close %closing obligatory figures that anova makes
         dWLESig = multcompare(dWLEStats, 'Display', 'off', 'CriticalValueType', 'bonferroni');
         % LRDescription = ['one-way anova', newline];
         behavDescription = ['repeated-measures anova', newline];
     end

     % LRDescription = [];
     for group = 1:size(dWLESig, 1)
         behavDescription = [behavDescription, [num2str(dWLESig(group,1)), 'v', num2str(dWLESig(group,2)), ': ', num2str(dWLESig(group,6)), newline]];
     end

    figure
    hold on
    multiFlydWLE = bar([1:size(allFlydWLE, 1)], nanmean(allFlydWLE, 2), 'FaceColor', 'flat', 'LineWidth', 2);
    % multiFlydWLE.CDAta
    errorbar([1:size(allFlydWLE, 1)], nanmean(allFlydWLE, 2), allFlydWLEerror, '.', 'Color', 'k')
    if doDots
        scatter([1:size(allFlydWLE, 1)], allFlydWLE', 10, 'k', 'filled')
    end
    xticks([1:nBlocks])
    xticklabels(xTitles)
    ylabel('dWLE') %change to something more informative later
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'box', 'off')

    set(gcf, 'Name', 'cross-fly dWLE', 'NumberTitle', 'off', 'Renderer', 'painters')

    figMax = get(gca, 'ylim');
    figMax(1) = [];
    text(4, figMax/2, behavDescription)
else
    disp('only one fly supplied, cross-fly analyses not applicable')
end

%farewell intrepid traveller