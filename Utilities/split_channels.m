close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

dirs = dir(mainDirectory);

addpath('D:\group_swinderen\Dinis\Fiji.app\scripts');

ImageJ;

%%

% this level is the list of dates
for date = [5:7 9]%setdiff(3:length(dirs),18)
   blocks = dir([mainDirectory dirs(date).name]);
   for block = 3:length(blocks)
       currentDirectory = [mainDirectory dirs(date).name '/' blocks(block).name];
       
       if ~exist([currentDirectory '/green_channel.raw'],'file')
           
           disp(currentDirectory);
       
           imageFile = dir([currentDirectory '/Image*']);

           ij.IJ.run('Raw...',['open=' currentDirectory '/' imageFile.name '  width=512 height=512 little-endian number=200000']);

           %de-interleave channels
           ij.IJ.run("Deinterleave", "how=2");

           %convert staks to hyperstacks
           ij.IJ.selectWindow([imageFile.name ' #1']);

           %green channel
           ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel.raw']);
%            ij.IJ.run("Stack to Hyperstack...",'"channels=1 slices=5 frames=11000');
%            ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all');
%            ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel_AVG.raw']);    
%            ij.IJ.run('Close');
           ij.IJ.run('Close');

            %red channel
    %        ij.IJ.selectWindow('Image_0001_0001.raw #2');

    %        %save red channel
%            ij.IJ.run('Raw Data...',['path=' currentDirectory '/red_channel.raw']);
           ij.IJ.run('Close');

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