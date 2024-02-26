function plot3D(data, visibility)

        figure('visible',visibility);
        h = slice(data,[],[],1:size(data,3),'nearest');
        set(h,'edgecolor','none');
        colorbar; colormap('jet');
        xlabel('X'); ylabel('Y'); zlabel('Z');
        
        set(gca,'ydir','reverse');
        set(gca,'zdir','reverse');
        set(gca,'xdir','reverse');
        
        set(gca,'Visible','off')
        
        % I just thought these looked good (I = Dinis)
        set(gca,'DataAspectRatio',[11.6667 10.6667 1]);
%         set(gca,'cameraposition',[781.4816 1.4335e+03 -5.8608]);
%         set(gca,'cameraposition',[112.2586 -225.5818 10.0540]);
        set(gca,'cameraposition',[-63.6363 -208.1706   10.4675]);
end