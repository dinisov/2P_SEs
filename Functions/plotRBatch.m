function plotRBatch(R, iter, outputDirectory)

    fieldNames = fieldnames(R.r);

    for i = 1:length(fieldNames)
        figure;
        imagesc(R.r.(fieldNames{i})); colorbar; colormap('jet'); caxis([-1 1]); title(upper(fieldNames{i}));
        xlabel('X'); ylabel('Y');
        saveas(gcf,fullfile(outputDirectory,[fieldNames{i} '_' num2str(iter) '.png']));
        close all;
    end

end