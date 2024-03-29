function  plotBrain(R, I, signLevel, outputDirectory, visibility)
%plotBrain Summary of this function goes here
%   Detailed explanation goes here

    fieldNames = fieldnames(R.r);
    
    for i = 1:length(fieldNames)
        
        name = extractAfter(fieldNames{i},'_');
    
        figure('visible',visibility); g = imshow(cat(3,I,I,I));%if image is rgb colormap does not apply
        hold on;
        h = image(xlim,ylim,ind2rgb(uint8(255 * (1+R.r.(fieldNames{i}))/2), jet(256))); colormap(jet); caxis([-1 1]);
        c = colorbar;
        c.Label.String = '\bf{Correlation coefficient}';
        set(h,'alphadata',R.p.(['p_' name]) < signLevel);
        saveas(gcf,fullfile(outputDirectory,[fieldNames{i} '_brain.png']));
        close all;
    end
    
end