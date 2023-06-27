function plotAR(meanDataSeq,threshold,I,outputDirectory)
%plotAR plot violations of repetitions and repetitions
%   Detailed explanation goes here

    %get the sequences
    AAAA = squeeze(meanDataSeq(:,16,:,:));
    AAAR = squeeze(meanDataSeq(:,8,:,:));
    RRRR = squeeze(meanDataSeq(:,1,:,:));
    RRRA = squeeze(meanDataSeq(:,9,:,:));

    %average across volumes
    AAAA = squeeze(mean(AAAA,1));
    AAAR = squeeze(mean(AAAR,1));
    RRRR = squeeze(mean(RRRR,1));
    RRRA = squeeze(mean(RRRA,1));

    %AAAA-AAAR
    figure; imagesc(AAAA-AAAR); title('AAAA-AAAR'); colorbar;
    saveas(gcf,[outputDirectory 'AAAA-AAAR.png']);
    figure; imagesc(RRRR-RRRA); title('RRRR-RRRA'); colorbar;
    saveas(gcf,[outputDirectory 'RRRR-RRRA.png']);

    figure; imshow(cat(3,I,I,I));%if image is rgb colormap does not apply
    hold on;
    h = image(xlim,ylim,ind2rgb(uint8(255 * (100+(AAAA-AAAR))/200), jet(256))); colormap(jet); caxis([-100 100]);
    c = colorbar;
    c.Label.String = '\bf{Difference}';
    set(h,'alphadata',abs(AAAA-AAAR) > threshold);
    title('AAAA-AAAR');
    saveas(gcf,[outputDirectory 'AAAA-AAAR' num2str(threshold) '.png'])

    figure; imshow(cat(3,I,I,I));%if image is rgb colormap does not apply
    hold on;
    h = image(xlim,ylim,ind2rgb(uint8(255 * (100+(RRRR-RRRA))/200), jet(256))); colormap(jet); caxis([-100 100]);
    c = colorbar;
    c.Label.String = '\bf{Difference}';
    set(h,'alphadata',abs(RRRR-RRRA) > threshold);
    title('RRRR-RRRA');
    saveas(gcf,[outputDirectory 'RRRR-RRRA' num2str(threshold) '.png']);

end

