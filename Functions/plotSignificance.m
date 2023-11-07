function  plotSignificance(AAAAvsAAAR, RRRRvsRRRA, I, outputDirectory, visibility)
%plotSignificance Summary of this function goes here
%   Detailed explanation goes here

    ploth(AAAAvsAAAR,'AAAAvsAAAR');
    ploth(RRRRvsRRRA,'RRRRvsRRRA');
    
    function ploth(data, name)
        data.h(isnan(data.h)) = 0;
        figure('visible',visibility); imshow(cat(3,I,I,I));%if image is rgb colormap does not apply
        hold on;
        h = image(xlim,ylim,ind2rgb(uint8(255 * data.h), [0.9290 0.6940 0.1250])); %caxis([-1 1]); %colormap(jet); 
        set(h,'alphadata',logical(data.h));
        saveas(gcf,fullfile(outputDirectory,[name '.png']));
        close all;        
    end
end