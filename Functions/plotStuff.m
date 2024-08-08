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

%Behav separated
if isfield( R , 'dataSeqBehav' ) %Note difference from usual architecture, because 'R' here is actually 'R(fly).BLOCK(b)', inherited 
    for statInd = 1:size( R.dataSeqBehav,2 )
        outputStateDirectory = [outputDirectory,filesep,'State_',num2str(R.dataSeqBehav(statInd).state)];
        disp(['Plotting stuff to ',outputStateDirectory])
        subDirectory = fullfile(outputStateDirectory,'Collapsed');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        
        %(Skip T-tests, because not implemented currently)
        
        if isfield(R.dataSeqBehav(statInd),'r2')
            % plot r_squared
            plotR2(R.dataSeqBehav(statInd).r2, subDirectory, visibility);
        end
        
        if isfield(R.dataSeqBehav(statInd),'r')
            % plot r
            plotR(R.dataSeqBehav(statInd).r, subDirectory, visibility);
            % plot r overlayed in fly's brain (thresholded by significance)
            signLevel = 0.05;
            plotBrain(R.dataSeqBehav(statInd).r, brainImage, signLevel, subDirectory, visibility);
        end
    end
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

%Behav separated
if isfield( R , 'dataSeqBehav' ) %Note difference from usual architecture, because 'R' here is actually 'R(fly).BLOCK(b)', inherited 
    for statInd = 1:size( R.dataSeqBehav,2 )
        nVol = size(R.meanDataSeq,1);

        outputStateDirectory = [outputDirectory,filesep,'State_',num2str(R.dataSeqBehav(statInd).state)];
        subDirectory = fullfile(outputStateDirectory,'Per volume');
        if ~exist(subDirectory,'dir')
           mkdir(subDirectory); 
        end
        
        for vol = 1:nVol
            thisVolDirectory = fullfile(subDirectory,num2str(vol));
            if ~exist(thisVolDirectory,'dir')
               mkdir(thisVolDirectory); 
            end

            % plot r_squared
            if isfield(R.dataSeqBehav(statInd),'r2Vol')
                plotR2(R.dataSeqBehav(statInd).r2Vol(vol), thisVolDirectory, visibility);
            end

            if isfield(R.dataSeqBehav(statInd),'rVol')
                % plot r
                plotR(R.dataSeqBehav(statInd).rVol(vol), thisVolDirectory, visibility);
                % plot r overlayed in fly's brain (thresholded by significance)
                signLevel = 0.05;
                plotBrain(R.dataSeqBehav(statInd).rVol(vol), brainImage, signLevel, thisVolDirectory, visibility);
            end
        
        end
        
    end
end

 
% %% r plots thresholded by significance level (not overlayed on fly's brain)
% 
% % figure; h = imagesc(r_slrp); colorbar; set(h,'alphadata',p_slrp); set(gca,'color','k'); caxis([-1 1]);
% % figure; h = imagesc(r_lrpr); colorbar; set(h,'alphadata',p_lrpr); set(gca,'color','k'); caxis([-1 1]);
% % figure; h = imagesc(r_weird); colorbar; set(h,'alphadata',p_weird); set(gca,'color','k'); caxis([-1 1]);
% % figure; h = imagesc(r_ephys); colorbar; set(h,'alphadata',p_ephys); set(gca,'color','k'); caxis([-1 1]);