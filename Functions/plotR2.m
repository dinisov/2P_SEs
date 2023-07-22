function plotR2(R2, outputDirectory)

    fieldNames = fieldnames(R2);

    for i = 1:length(fieldNames)
        figure; 
        imagesc(R2.(fieldNames{i})); colorbar; colormap('jet'); title(upper(fieldNames{i}));
        xlabel('X'); ylabel('Y');
        saveas(gcf,fullfile(outputDirectory,[fieldNames{i} '.png']));
        close all;
    end

end