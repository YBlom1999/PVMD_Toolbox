function h = flatplot3(V,F,C,h)
%2D plot of icohemishpere. Facet color C.  3D cartesian coordinates are
%converted to 2D cylincer coordinates first.

if isnumeric(h)                     %if handle is just a number
    Vcyl = cart2cyl(V);             %convert 3D to 2D coordinates
    figure(h);
    clf
    h = patch('Vertices',Vcyl,'Faces',F,...
        'FaceVertexCData',C,'FaceColor','flat');    %plot in 2D coordinates
    %---make it pretty---
    axis equal off
    shading flat
    colormap(parula(512))
    caxis([0,1000])
    hc = colorbar;
    set(get(hc,'Title'),'string','Irradiance [W/(m^2\cdot SR)]','FontSize',12)
    text(  0, 95,'North','HorizontalAlignment','center','FontSize',15)
    text( 95,  0,'East' ,'HorizontalAlignment','center','Rotation',-90,'FontSize',15)
    text(  0,-95,'South','HorizontalAlignment','center','FontSize',15)
    text(-95,  0,'West' ,'HorizontalAlignment','center','Rotation',90,'FontSize',15)
    ax =gca;
    ax.FontSize = 15;
else                                %if h is a real handle of existing figure
    set(h,'FaceVertexCData',C);     %just update the color
end
drawnow
%..........................................................................
    function Vcyl = cart2cyl(Vcart)
        %convert 3D cartesian to 2D cylinder coordinates
        %used to convert icohemisphere vertex coordinates
        %(not the light source, which is at the center of each vertex)
        zeni = atand(sqrt(Vcart(:,1).^2 + Vcart(:,2).^2)./Vcart(:,3));
        azi = atan2d(-Vcart(:,1),-Vcart(:,2));
        
        Vcyl_x = zeni .* -sind(azi);
        Vcyl_y = zeni .* -cosd(azi);
        Vcyl_z = zeros(size(Vcyl_x));               %z-coordinate is set to 0
        Vcyl = [Vcyl_x,Vcyl_y,Vcyl_z];              %combine
    end

end