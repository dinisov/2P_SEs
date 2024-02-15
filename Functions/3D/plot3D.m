function plot3D(data, visibility)

        figure('visible',visibility);
        h = slice(data,[],[],1:size(data,3),'nearest');
        set(h,'edgecolor','none');
        colorbar; colormap('jet');
        xlabel('X'); ylabel('Y'); zlabel('Z');
        
        % I just thought these looked good (I = Dinis)
        set(gca,'DataAspectRatio',[11.6667 10.6667 1]);
        set(gca,'cameraposition',[-87.2082 255.7671 9.8862]);

end