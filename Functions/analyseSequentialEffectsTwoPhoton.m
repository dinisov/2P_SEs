function R = analyseSequentialEffectsTwoPhoton(blocks, aux_plots)

    %% sort data according to previous sequence

    n_back = 5; n_seq = 2^n_back; %TODO: make n_back variable
    
    n_blocks = length(blocks);
    
    blocks = sortSEs2P(blocks, n_back);

    %% concatenate blocks
 
    total_length = 0;
 
    for b = 1:n_blocks
        total_length = total_length + size(blocks(b).rDataSeq,5);
    end 
     
    % concatenate images from different blocks
    auxSize = size(blocks(1).rDataSeq);
    allDataSeq = zeros([auxSize(1:4) total_length]);
     
    start_index = 1;
 
    % these are already separated by sequence so in order to group by block
    % it is only necessary to stack along FIFTH dimension
    for b = 1:n_blocks

        allDataSeq(:,:,:,:,start_index + 1:start_index + size(blocks(b).rDataSeq,5)) = blocks(b).rDataSeq;

        start_index = start_index + size(blocks(b).rDataSeq,5);

    end
     
    %%
     
    % join images corresponding to the same pattern (01001 and 10110 and so on)
    allDataSeq = allDataSeq + flip(allDataSeq,2);
    allDataSeq = allDataSeq(:,1:16,:,:,:);
    
    % in preparation for calculating the nan mean
    allDataSeq(allDataSeq == 0) = nan;
    
    % mean across fifth (trial) dimension
    meanDataSeq = mean(allDataSeq, 5, 'omitnan');
    
    %number of images for each type of sequence
    nData = sum(~isnan(allDataSeq(1,:,1,1,:)), 5);
    
    % reorder according to the literature
    meanDataSeq = meanDataSeq(:,seq_eff_order(n_back),:,:);
%     semERPs = semERPs(:,seq_eff_order(n_back));
    
    % we can make a figure showing the average image as a function of the
    % sequence, for each time point
%     figure;
%     plot(meanERPs);       

    % put all results into a neat structure
    R = struct;
     
    R.meanDataSeq = meanDataSeq;
       
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