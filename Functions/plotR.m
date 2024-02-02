function plotR(R, outputDirectory, visibility)

    fieldNames = fieldnames(R.r);

    for i = 1:length(fieldNames)
        figure('visible',visibility);
        imagesc(R.r.(fieldNames{i})); colorbar; colormap('jet'); caxis([-1 1]); title(upper(fieldNames{i}));
        xlabel('X'); ylabel('Y');
        saveas(gcf,fullfile(outputDirectory,[fieldNames{i} '.png']));
        close all;
    end

end