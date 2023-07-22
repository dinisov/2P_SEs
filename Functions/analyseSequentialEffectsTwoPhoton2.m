function R = analyseSequentialEffectsTwoPhoton2(imageStack,randomSequence,nVol)

    %% sort data according to previous sequence

    n_back = 5; 
    
    % yields 5D matrix with (vol,seq,pixelX,pixelY,trial)
    dataSeq = sortSEs2P2(imageStack,randomSequence,nVol);
    
    % mean across fifth (trial) dimension (much faster than nan mean)
    meanDataSeq = sum(dataSeq,5)./sum(dataSeq~=0,5);
    
    % this fixes the possibility that all images in a stack have zero
    % values for a particular pixel (usually from alignment)
    meanDataSeq(isnan(meanDataSeq)) = 0;
    
    %number of images for each type of sequence
    nData = sum(dataSeq(1,:,16,16,:)~=0, 5);
    
    % reorder according to the literature
    meanDataSeq = meanDataSeq(:,seq_eff_order(n_back),:,:);
%     semERPs = semERPs(:,seq_eff_order(n_back));
    
    % we can make a figure showing the average image as a function of the
    % sequence, for each time point
%     figure;
%     plot(meanERPs);       

    % put all results into a neat structure
    R = struct;
     
    R.dataSeq = dataSeq(:,seq_eff_order(n_back),:,:,:);
    R.meanDataSeq = meanDataSeq;
    R.nData = nData;
       
%     %make a tick every x steps
%     y_ticks = 1:2:length(f);
%     y_tick_labels = f(y_ticks);
%     
%     % create 
%     figure; imagesc(squeeze(magnitudeSEs(:,16,:))-squeeze(magnitudeSEs(:,8,:)),'xdata',time_bounds); colorbar;
%     set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
%     xlabel('time (ms)'); ylabel('Frequency (Hz)');
%     title('AAAA_minus_AAAR');
% 
%     figure; imagesc(squeeze(magnitudeSEs(:,9,:))-squeeze(magnitudeSEs(:,1,:)),'xdata',time_bounds); colorbar;
%     set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
%     xlabel('time (ms)'); ylabel('Frequency (Hz)');
%     title('RRRR_minus_RRRA');
    
end