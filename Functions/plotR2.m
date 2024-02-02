function plotR2(R2, outputDirectory,visibility)

    fieldNames = fieldnames(R2);

    for i = 1:length(fieldNames)
        figure('visible',visibility);
        imagesc(R2.(fieldNames{i})); colorbar; colormap('jet'); title(upper(fieldNames{i}));
        xlabel('X'); ylabel('Y');
        saveas(gcf,fullfile(outputDirectory,[fieldNames{i} '.png']));
        close all;
    end

end