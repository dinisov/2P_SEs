function plotR(R, outputDirectory)

    fieldNames = fieldnames(R.r);

    for i = 1:length(fieldNames)

        name = extractAfter(fieldNames{i},'_');

        figure;
        imagesc(R.r.(fieldNames{i})); colorbar; colormap('jet'); caxis([-1 1]); title(upper(name));
        xlabel('X'); ylabel('Y');
        saveas(gcf,[outputDirectory name '_r.png']);

    end

end