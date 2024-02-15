function plotStuff3D(R, outputDirectory, visibility)

brainImage = R.brainImage;

%% collapsed across volumes (time)

subDirectory = fullfile(outputDirectory,'Collapsed');
if ~exist(subDirectory,'dir')
   mkdir(subDirectory); 
end

% plot t-test results for AAAA vs AAAR and RRRR vs RRRA
if isfield(R,'AAAAvsAAAR')
    plotSignificance(R.AAAAvsAAAR,R.RRRRvsRRRA,brainImage,subDirectory,visibility);
end

% plot r_squared
if isfield(R,'r2')
    
    fieldNames = fieldnames(R.r2);
    
    for i = 1:length(fieldNames)
        plot3D(R.r2.(fieldNames{i}), visibility);
        title(upper(fieldNames{i})); 
        saveas(gcf,fullfile(subDirectory,[fieldNames{i} '.png']));
        close all;
    end
end

% plot r (collapsed across t)
if isfield(R,'r')
    
    fieldNames = fieldnames(R.r.r);
    
    for i = 1:length(fieldNames)
        plot3D(R.r.r.(fieldNames{i}), visibility);
        title(upper(fieldNames{i})); caxis([-1 1]); 
        saveas(gcf,fullfile(subDirectory,[fieldNames{i} '.png']));
        close all;
    end

        % plot r overlayed in fly's brain (thresholded by significance)
    %     signLevel = 0.05;
    %     plotBrain(R.r, brainImage, signLevel, subDirectory, visibility);

end

%% per volume

nVol = size(R.meanDataSeq,1);

subDirectory = fullfile(outputDirectory,'Per volume');
if ~exist(subDirectory,'dir')
   mkdir(subDirectory); 
end

for vol = 1:nVol
    thisVolDirectory = fullfile(subDirectory,num2str(vol));
    if ~exist(thisVolDirectory,'dir')
       mkdir(thisVolDirectory); 
    end
%     thisVolData = permute(squeeze(BLOCK.meanDataSeq(vol,:,:,:)),[2,3,1]);

    % plot t-test results for AAAA vs AAAR and RRRR vs RRRA
    if isfield(R,'AAAAvsAAARVol')
        plotSignificance(R.AAAAvsAAARVol(vol),R.RRRRvsRRRAVol(vol),brainImage,thisVolDirectory, visibility);
    end

    % plot r_squared
    if isfield(R,'r2Vol')
        fieldNames = fieldnames(R.r2Vol(vol));
    
        for i = 1:length(fieldNames)
            plot3D(R.r2Vol(vol).(fieldNames{i}), visibility);
            title(upper(fieldNames{i}));
            saveas(gcf,fullfile(thisVolDirectory,[fieldNames{i} '.png']));
            close all;
        end
    end

    % plot r
    if isfield(R,'rVol')
        fieldNames = fieldnames(R.rVol(vol).r);
    
        for i = 1:length(fieldNames)
            plot3D(R.rVol(vol).r.(fieldNames{i}), visibility);
            title(upper(fieldNames{i})); caxis([-1 1]); 
            saveas(gcf,fullfile(thisVolDirectory,[fieldNames{i} '.png']));
            close all;
        end

        % plot r overlayed in fly's brain (thresholded by significance)
%         signLevel = 0.05;
%         plotBrain(R.rVol(vol), brainImage, signLevel, thisVolDirectory, visibility);
    end
    
end