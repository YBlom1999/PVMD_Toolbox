function [V,F] = rectabox(x,y,z,u,v,w)
%create a rectangular box with corner coordinates xyz and length
%width height uvw. It is assumed that uvw are positive. Facets:
%bot, north, east, south, west, top have normal pointing outward.
%If u==0 or v==0 or w==0 a 2D rectangle is created.

if u == 0 %make rectangle perpendicular to x-axis
    V = [x y z; x y+v z; x y+v z+w; x y z+w];
    F = [1 2 3 4]; %normal pointing in positive x-direction
elseif v == 0 %make rectangle perpendicular to y-axis
    V = [x y z; x y z+w; x+u v z+w; x+u y z];
    F = [1 2 3 4]; %normal pointing in positive y-direction
elseif w == 0 %make rectangle perpendicular to z-axis
    V = [x y z; x+u y z; x+u y+v z; x y+v z];
    F = [1 2 3 4]; %normal pointing in positive z-direction
else    %make real 3D box
    %define vertices (start at uvw)
    V = [x y z; x y+v z; x+u y+v z; x+u y z; x y z+w; x y+v z+w; x+u y+v z+w; x+u y z+w];
    %facets (bot, north, east, south, west, top) normal pointing outward
    F = [1 2 3 4; 6 7 3 2; 7 8 4 3; 8 5 1 4; 5 6 2 1; 8 7 6 5];
end
end
