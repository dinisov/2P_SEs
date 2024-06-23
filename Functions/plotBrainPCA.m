function  plotBrainPCA(R, I,visibility)
%plotBrain Summary of this function goes here
%   Detailed explanation goes here

    figure('visible',visibility); g = imshow(cat(3,I,I,I));%if image is rgb colormap does not apply
    hold on;
    
    % R mapped to [0 1]
    R_plot = R - min(R,[],'all');
    R_plot = R_plot./max(R_plot,[],'all');
    
    h = image(xlim,ylim,ind2rgb(uint8(255 * R_plot), jet(256))); colormap(jet); caxis([min(R,[],'all') max(R,[],'all')]);
    colorbar;
%     c.Label.String = '\bf{Correlation coefficient}';
    set(h,'alphadata',abs(R./max(abs(R),[],'all')));
    set(g.Parent,'position',[0 0 .9 .9]);
    set(h.Parent,'position',[0 0 .9 .9]);
%         saveas(gcf,fullfile(outputDirectory,[fieldNames{i} '_brain.png']));
%     close all;
end