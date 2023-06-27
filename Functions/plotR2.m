function plotR2(R2, outputDirectory)

    fieldNames = fieldnames(R2);

    for i = 1:length(fieldNames)

        name = extractAfter(extractAfter(fieldNames{i},'_'),'_');

        figure; 
        imagesc(R2.(fieldNames{i})); colorbar; colormap('jet'); title(upper(name));
        xlabel('X'); ylabel('Y');
        saveas(gcf,[outputDirectory name '_r2.png']);

    end

end