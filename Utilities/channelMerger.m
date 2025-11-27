%Utility designed to merge green and red channel data (e.g. GCamP and reporter)

close all; clear;

rdmDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';
blocks = readtable("I:\RFDG2021-Q4413\2P Record\2P_record");

showDebugInfo = 1;

%%

chosenFlies = [394];
chosenBlocks = {[2]}; % leave empty if merging all blocks for fly/s

%%

flagParamSaveList = who;
flagParamSaveList = [flagParamSaveList;'flagParamSaveList';'fly'];


for fly = 1:length(chosenFlies)
    %Below borrowed from reduce_files
    clearvars('-except', flagParamSaveList{:});
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);
    if showDebugInfo
        thisFlyBlocks
    end

    if ~isempty(chosenBlocks) && ~isempty(chosenBlocks{fly})
        thisFlyBlocks = thisFlyBlocks(ismember(thisFlyBlocks.Block,chosenBlocks{fly}),:);
    end

    nBlocks = height(thisFlyBlocks);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));

    for b = chosenBlocks{fly}
            currentBlock = thisFlyBlocks( find( thisFlyBlocks.Block == b ) ,:); %New
            if showDebugInfo
                currentBlock
            end 
            %QA
            if isempty( currentBlock )
                ['## No blocks found for fly #', num2str(chosenFlies(fly)),' block ', num2str(b)]
                crash = yes
            end
            flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
            currentRDMDirectory = fullfile(rdmDirectory,currentDate,flyID);
            disp(['Fly: ',flyID]);
            if showDebugInfo
                currentRDMDirectory
            end
            
            codeStartTime = posixtime(datetime('now'));

            %yofukashi

            %Check for existence of green and red channel data
            if exist([currentRDMDirectory,filesep,'avg_z_green_aligned.mat']) ~= 2 || ...
                    exist([currentRDMDirectory,filesep,'avg_z_red_aligned.mat']) ~= 2
                ['## Alert: Either green or red avg aligned data missing ##']
                crash = yes
            end

            %Load data
            temp = cell(2,1);
            for chan = 1:2
                if chan == 1
                    temp{chan} = load([currentRDMDirectory,filesep,'avg_z_green_aligned.mat']);
                elseif chan == 2
                    temp{chan} = load([currentRDMDirectory,filesep,'avg_z_red_aligned.mat']);
                end
            end

            %DEBUG ONLY IMAGE TRIMMING
            input('Ack confirm image trimming')
            temp{1}.avg_z_green_aligned = temp{1}.avg_z_green_aligned(8:end-8,8:end-8,:);
            temp{2}.avg_z_red_aligned = temp{2}.avg_z_red_aligned(8:end-8,8:end-8,:);


            data.green = temp{1}.avg_z_green_aligned;
            data.red = temp{2}.avg_z_red_aligned;
            clear temp

            %Testatory plot
            figure
            subplot(2,1,1)
            plot(squeeze(nanmean(data.green,[1,2])),'Color','g')
            title('Green raw trace')
            subplot(2,1,2)
            plot(squeeze(nanmean(data.red,[1,2])),'Color','r')
            title('Red raw trace')


            %Subtract
            data.merged = data.green - data.red;
            data.mergedFlat = squeeze(nanmean(data.green,[1,2])) - smooth( squeeze(nanmean(data.red,[1,2])) , 5);

            %Another figure
            figure
            %plot(squeeze(nanmean(data.merged,[1,2])),'Color','k')
            hold on
            plot(squeeze(nanmean(data.green,[1,2])-nanmean(data.green,'all')),'Color','g')
            plot(squeeze(nanmean(data.red,[1,2])-nanmean(data.red,'all')),'Color','r')
            plot(squeeze(nanmean(data.merged,[1,2])),'Color','k')
            hold off
            title('Merged raw trace')

            %And a video
            %(Scaled)
            figure
            for i = 1:size(data.merged,3)
                subplot(2,2,2)
                imagesc( data.green(:,:,i) )
                title('Green')
                subplot(2,2,4)
                imagesc( data.red(:,:,i) )
                title('Red')
                subplot(2,2,[1,3])
                imagesc( data.merged(:,:,i) )
                title(['Merged (Frame ',num2str(i),')'])
                drawnow
            end

            %And a video
            %(Unscaled)
            temp = data.merged;
            %temp = temp * 25;
            temp = temp + 128;
            figure
            for i = 1:size(data.merged,3)
                subplot(2,2,2)
                image( data.green(:,:,i) )
                title('Green')
                subplot(2,2,4)
                image( data.red(:,:,i) )
                title('Red')
                subplot(2,2,[1,3])
                image( temp(:,:,i) )
                title(['Merged w/ manual scaling (Frame ',num2str(i),')'])
                drawnow
            end





    end


end