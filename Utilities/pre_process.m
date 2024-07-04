close all;

if ~exist('IJM','var')
    clear;
    addpath('D:\group_vanswinderen\Dinis\Fiji.app\scripts');
    ImageJ;
end

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable("D:\group_vanswinderen\Dinis\2P Record\2P_record");

%get rid of excluded flies
% blocks = blocks(~logical(blocks.Exclude),:);

%%
chosenFlies = [200:205,188];

chosenBlocks = {[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2,3],[1,2]};

for fly = 1:length(chosenFlies)
    
    BLOCKS = struct;
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);
    
    if ~isempty(thisFlyBlocks)
    
        currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));

        nBlocks = height(thisFlyBlocks);

        for b = chosenBlocks{fly}

            currentBlock = thisFlyBlocks(b,:);
            flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
            currentDirectory = fullfile(mainDirectory,currentDate,flyID);

            nSlices = currentBlock.Steps + currentBlock.FlybackFrames;
            totalFrames = currentBlock.realFrames;
            trials = totalFrames/nSlices;

%             disp(currentDirectory);
    
           if ~exist([currentDirectory '/green_channel.raw'],'file')

               disp(currentDirectory);

               imageFile = dir([currentDirectory '/Image*']);

               ij.IJ.run('Raw...',['open=' currentDirectory '/' imageFile.name '  width=512 height=512 little-endian number=200000']);

               if currentBlock.nChannels == 2
                   
                   %de-interleave channels
                   ij.IJ.run("Deinterleave", "how=2");
                   
                   ij.IJ.selectWindow([imageFile.name ' #1']);
                   
                   %save green channel
                   ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel.raw']);
                   
                   if ~exist([currentDirectory '/brain.jpg'],'file')             
                       ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all');
                       ij.IJ.run('Jpeg...',['path=' currentDirectory '/brain.jpg']);  
                       ij.IJ.run('Close');
                   end

                   ij.IJ.run('Close');
                   
                   ij.IJ.selectWindow([imageFile.name ' #2']);
                   
                   %save red channel
                   ij.IJ.run('Raw Data...',['path=' currentDirectory '/red_channel.raw']);
                   ij.IJ.run('Close');
                   
               else
               
                   %green channel
                   ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel.raw']);

    %                ij.IJ.run('Close');

    %                ij.IJ.run("Stack to Hyperstack...",['"channels=1 slices=' num2str(nSlices) ' frames=' num2str(trials)]);
    %                ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all');
    %                ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel_AVG.raw']);

    %                ij.IJ.selectWindow('green_channel.raw');
    %                ij.IJ.run('Close');

                    if ~exist([currentDirectory '/brain.jpg'],'file')             
                       ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all');
                       ij.IJ.run('Jpeg...',['path=' currentDirectory '/brain.jpg']);  
                       ij.IJ.run('Close');
                    end

                   ij.IJ.run('Close');
                   
               end

           else
               disp(['-# Preprocessed data already existing for ',flyID,' #-'])
           end
        end
    else
        disp(['-# No blocks found for fly ',num2str(chosenFlies(fly)),' #-'])
    end
end