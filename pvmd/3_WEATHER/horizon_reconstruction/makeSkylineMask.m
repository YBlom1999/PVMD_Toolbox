function skylineMask = makeSkylineMask(TOOLBOX_input,skydome)
% makeSkylineMask makes the skylineMask to account for horizon effects
%
% This function makes a mask that can be used to adjust the Sensitivity map
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% skydome: struct
%   Structure containing information on the skydome
%
% Returns
% -------
% skylineMask: double
%   The mask of the skyline that indicates whether the sky element is
%   blocked.
%
% Developed by unknown (E. Garcia?) Adjusted by Y. Blom. Commented by A. Alcaniz


load(TOOLBOX_input.irradiation.skyline_file,'skyline')
skyline_az = skyline{1}-180;
skyline_alt = skyline{2};
Azimuth = skydome.AZA(:,1);
Altitude = 90- skydome.AZA(:,2);

skylineMask = zeros(length(skydome.AZA),1);
for i = 1:length(skylineMask)
    Angle_skyline = interp1(skyline_az,skyline_alt,Azimuth(i));
    if  Angle_skyline < Altitude(i)
        skylineMask(i) = 1;
    end
end
end