close all;

if ~exist('IJM','var')
    clear;
    addpath('D:\group_swinderen\Dinis\Fiji.app\scripts');
    ImageJ;
end

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

blocks = readtable('../../2P Record/2P_record');

%get rid of excluded flies
blocks = blocks(~logical(blocks.Exclude),:);

%%

chosenFlies = 20:41;

for fly = 1:length(chosenFlies)
    
    BLOCKS = struct;
    
    % the blocks corresponding to this fly
    thisFlyBlocks = blocks(blocks.Fly == chosenFlies(fly),:);
    
    currentDate = char(datetime(thisFlyBlocks.Date(1),'Format','dMMMyy'));
    
    nBlocks = height(thisFlyBlocks);
    
    for b = 1:nBlocks
        
        currentBlock = thisFlyBlocks(b,:);
        flyID = ['fly' num2str(currentBlock.FlyOnDay) '_exp' num2str(currentBlock.Block) '_' currentDate];
        currentDirectory = fullfile(mainDirectory,currentDate,flyID);
        
        nSlices = currentBlock.Steps + currentBlock.FlybackFrames;
        totalFrames = currentBlock.realFrames;
        trials = totalFrames/nSlices;
        
%         disp(currentDirectory);
       
       if exist([currentDirectory '/green_channel.raw'],'file') && ~exist([currentDirectory '/brain.jpg'],'file')
           
           disp(currentDirectory);
       
%            imageFile = dir([currentDirectory '/Image*']);

           ij.IJ.run('Raw...',['open=' currentDirectory '/green_channel.raw  width=512 height=512 big-endian number=' num2str(totalFrames)]);

           %de-interleave channels
%            ij.IJ.run("Deinterleave", "how=2");

           %convert staks to hyperstacks
%            ij.IJ.selectWindow([imageFile.name ' #1']);

           %green channel
%            ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel.raw']);
           ij.IJ.run("Stack to Hyperstack...",['"channels=1 slices=' num2str(nSlices) ' frames=' num2str(trials)]);
           ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all');
%            ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel_AVG.raw']);
           ij.IJ.selectWindow('green_channel.raw');
           ij.IJ.run('Close');
           ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all');
           ij.IJ.run('Jpeg...',['path=' currentDirectory '/brain.jpg']);  
           
           ij.IJ.run('Close');
           ij.IJ.run('Close');

            %red channel
    %        ij.IJ.selectWindow('Image_0001_0001.raw #2');

    %        %save red channel
%            ij.IJ.run('Raw Data...',['path=' currentDirectory '/red_channel.raw']);
%            ij.IJ.run('Close');

    %        ij.IJ.run("Stack to Hyperstack...",'"channels=1 slices=4 frames=50');
    % %        ij.IJ.setTool("Rectangle");
    %        img = ij.IJ.getImage();
    %        img.setRoi(0,0,512,512);
    %        ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all')
    %        ij.IJ.run("Align slices in stack...",'"Matching method"="Normalized correlation coefficient"');

    %        ij.IJ.selectWindow('Image_0001_0001.raw #2');
    %        ij.IJ.run("Stack to Hyperstack...",'channels=1 slices=4 frames=50');
       
       end   
   end 
end