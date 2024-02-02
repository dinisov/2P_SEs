function plotStuff(R, outputDirectory, visibility)

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
    plotR2(R.r2, subDirectory, visibility);
end

if isfield(R,'r')
    % plot r
    plotR(R.r, subDirectory, visibility);
    
    % plot r overlayed in fly's brain (thresholded by significance)
    signLevel = 0.05;
    plotBrain(R.r, brainImage, signLevel, subDirectory, visibility);
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
        plotR2(R.r2Vol(vol), thisVolDirectory, visibility);
    end

    if isfield(R,'rVol')
        % plot r
        plotR(R.rVol(vol), thisVolDirectory, visibility);

        % plot r overlayed in fly's brain (thresholded by significance)
        signLevel = 0.05;
        plotBrain(R.rVol(vol), brainImage, signLevel, thisVolDirectory, visibility);
    end
    
end

 
% %% r plots thresholded by significance level (not overlayed on fly's brain)
% 
% % figure; h = imagesc(r_slrp); colorbar; set(h,'alphadata',p_slrp); set(gca,'color','k'); caxis([-1 1]);
% % figure; h = imagesc(r_lrpr); colorbar; set(h,'alphadata',p_lrpr); set(gca,'color','k'); caxis([-1 1]);
% % figure; h = imagesc(r_weird); colorbar; set(h,'alphadata',p_weird); set(gca,'color','k'); caxis([-1 1]);
% % figure; h = imagesc(r_ephys); colorbar; set(h,'alphadata',p_ephys); set(gca,'color','k'); caxis([-1 1]);