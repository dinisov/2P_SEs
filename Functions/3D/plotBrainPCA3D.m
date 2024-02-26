function  plotBrainPCA3D(data, brain,visibility)
%plotBrain Summary of this function goes here
%   Detailed explanation goes here

    brain = brain - min(brain,[],'all');
    brain = brain./max(brain,[],'all');

    figure('visible',visibility);
    h1 = slice(brain,[],[],1:size(brain,3),'nearest');
    colorbar('Visible','off');
    set(h1,'edgecolor','none');
    ax1 = gca;
    axes;
    h2 = slice(data,[],[],1:size(data,3),'nearest');
    colorbar;
    set(h2,'edgecolor','none');
    ax2 = gca;
    
    % flip axes
    set(ax1,'ydir','reverse','zdir','reverse','xdir','reverse');
    set(ax2,'ydir','reverse','zdir','reverse','xdir','reverse');
    
    % get rid of axes
    set(ax1,'Visible','off'); set(ax2,'Visible','off');
    
    % set camera position and data aspect ratio
    setCamera(ax1); setCamera(ax2);
    
    % gray for brain and jet for data
    set(ax1,'Colormap',gray); 
    set(ax2,'Colormap',jet);
    
    set(h2,'FaceColor','texturemap');
    set(h2,'FaceAlpha','texture');

%     data(data < 0.7) = 0;

    alphaChannel = abs(data./max(abs(data),[],'all'));
    
    for i = 1:size(data,3)
        h2(i).AlphaData = alphaChannel(:,:,i);
    end
    
    function setCamera(ax)
        set(ax,'DataAspectRatio',[11.6667 10.6667 .5]);
%         set(ax,'cameraposition',[781.4816 1.4335e+03 -5.8608]); 
%         set(ax,'cameraposition',[112.2586 -225.5818 10.0540]);
%         set(ax,'cameraposition',[220.2984 -547.1745 9.7920]);
        set(ax,'cameraposition',[-137.8676 -548.0298    9.7501]);
    end

end