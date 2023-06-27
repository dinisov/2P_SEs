close all; clear;

mainDirectory = '\\uq.edu.au\uq-inst-gateway1\RFDG2021-Q4413\2P_Data\Gcamp7s_CC\';

% dirs = dir(mainDirectory);

flies = readtable('../2P Record/2P_record');

% addpath('D:\group_swinderen\Dinis\Fiji.app\scripts');
% 
% ImageJ;

%%

% this level is the list of dates
for fly = 55:58%height(flies)
    
    currentFly = flies(fly,:);
    
    currentDate = char(datetime(currentFly.Date,'Format','dMMMyy'));
    
    currentFlyDirectory = ['fly' num2str(currentFly.FlyOnDay) '_exp' num2str(currentFly.Block) '_' currentDate];
       
    currentDirectory = fullfile(mainDirectory,currentDate,currentFlyDirectory);
       
    % open green channel
    ij.IJ.run('Raw...',['open=' currentDirectory '/green_channel.raw  width=512 height=512 big-endian number=' num2str(currentFly.realFrames)]);
%     ij.IJ.run('Raw...',['open=' currentDirectory '/green_channel.raw  width=512 height=512 big-endian number=110']);

    %calculate hyperstack parameters (for starters including flyback)
    slices = currentFly.Steps + currentFly.FlybackFrames;
    frames = currentFly.realFrames/slices;
%     slices = 5;
%     frames = 22;
    
    % build hyperstack
    ij.IJ.run("Stack to Hyperstack...",['"channels=1 slices=' num2str(slices) ' frames=' num2str(frames)]);
    
    %Z project
    ij.IJ.run('Z Project...','"Projection type"="Average Intensity" all');
    
    img = ij.IJ.getImage();
    img.setRoi(0,0,currentFly.pixelX,currentFly.pixelY);
    ij.IJ.run("Align slices in stack...",'"Matching method"="Normalized correlation coefficient"');
    
    
    ij.IJ.run('Raw Data...',['path=' currentDirectory '/green_channel_AVG_aligned.raw']);
    ij.IJ.run('Close All');
       
end