function [V,F,T] = unitcell(Vm,DS,DF,MM,GA,GH,GD)
%create unit cell around module with periodic boundary condtions
%and ceiling as light source.

xmid = (max(Vm(:,1)) + min(Vm(:,1)))/2;     %mid-point x-coord
xmin = xmid - DS/2;                         %min x

ymid = (max(Vm(:,2)) + min(Vm(:,2)))/2;     %mid-point y-coord
ymin = ymid - DF/2;                         %min y

zmin = 0;                                   %ground is always at zero
zmax = max(Vm(:,3))+MM;                     %max z

[V,F] = rectabox(xmin,ymin,zmin,DS,DF,zmax);
%TODO: check whether box is larger than the module. User can make
%box arbitraly small.

%---floor---
T(1).Facet = 1;
T(1).RT = [GA 0];         %reflected albedo
T(1).Scat = [GH,GD];      %haze and diffuse exponent
%---front-back walls---
T(2).Facet = [2 4];
T(2).RT = [0 1];          %transparent
T(2).Teleport = -DF;      %with periodic boundary conditions
T(2).Invis = 1;           %make invisible in plot
%---left-right walls---
T(3).Facet = [3 5];
T(3).RT = [0 1];          %transparent
T(3).Teleport = -DS;      %with periodic boundary conditions
T(3).Invis = 1;           %make invisible in plot
%---ceiling---
T(4).Facet = 6;
T(4).RT = [0 0];          %absorbing
T(4).Emit = [0 0 0];      %light source: theta, phi and nr rays.
T(4).Invis = 1;           %make invisible in plot
end