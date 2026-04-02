function angle_albedo = angle_Albedo_finder(height,length,width,Mod_tilt,Sun_Azi)
%angle_Albedo_finder Finds the minimum angle that can still reach the back
%side of the module via albedo reflection.
%
% This function calculates the minimum angle that can still reach the
% module via albedo reflection
%
% Parameters
% ----------
% height : double
%   height of the module
% length : double
%   length of the module
% width : double
%   width of the module
% Mod_tilt : double
%   tilt of the module
% Sun_Azi : double
%   azimuth of the sun
%
% Returns
% -------
% angle_albedo : double
%   The minimum angle for the sun that can still reach the back of the
%   module via albedo reflection
%
% Developed by Y. Blom
length_actual = length*cosd(Mod_tilt);
Angle_side = atand(width/length_actual);
if abs(Sun_Azi) <= Angle_side
    heigth_angle = 2*height-sind(Mod_tilt)*0.5*length;
    length_angle = 0.5*length_actual/cosd(Sun_Azi);
elseif Sun_Azi > Angle_side && Sun_Azi <= 180-Angle_side
    heigth_angle = 2*height - cotd(Sun_Azi)*0.5*width*sind(Mod_tilt);
    length_angle = 0.5*width/sind(Sun_Azi);
elseif Sun_Azi > -180+Angle_side && Sun_Azi <= - Angle_side
    heigth_angle = 2*height - cotd(-Sun_Azi)*0.5*width*sind(Mod_tilt);
    length_angle = 0.5*width/sind(-Sun_Azi);
else
    heigth_angle = 2*height+sind(Mod_tilt)*0.5*length;
    length_angle = 0.5*length_actual/cosd(Sun_Azi-180);
end
angle_albedo 