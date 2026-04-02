clear all;
close all;

%% Input data
DNI = 0;
DHI = 100;

SunAz = 0;
SunAlt = 45;

Year = 2023;
Month = 1;
Day = 22;

ForceSkyType = 0; %If you want to simulate a specific sky type, you can set it to 1.
skytype_force = 1; %1 is clear sky, 2 is partly clouded, 3 is fully clouded.

spectra_choice = 2; %1 is SMARTS, 2 is SBDARTS

%% Load data
filename_MODULE_simulation = 'MonoPERC.mat';
load('Hemisphere_data.mat','Hemisphere_data','Vs','Fs'); %load data of the hemisphere
load(filename_MODULE_simulation,'MODULE_output');
SM_compact = squeeze(mean(MODULE_output.SM_f(:,:,end,:),2));

%% Make skymap
day_year = pvl_date2doy(Year, Month, Day);
day_angle = 2*pi*(day_year - 1)/365.25;
E0 = 1.00011 + 0.034221*cos(day_angle) + 0.00128*sin(day_angle) + ...
    0.000719*cos(2*day_angle) + 0.000077*sin(2*day_angle);
load('constants/weather_params.mat', 'AM0')
extra_sol_power = E0*AM0;


[skymap,ind_sun,skytype] = perez_model(SunAz, SunAlt, DNI, DHI, ...
    Hemisphere_data(:,1), Hemisphere_data(:,2), Hemisphere_data(:,3), extra_sol_power);

if ForceSkyType
    skytype = skytype_force;
end


flatplot3(Vs,Fs,skymap,1);

%% Calculate Spectral irradiance
air_mass = 1./sind(SunAlt); %The air mass is defined as 1/sin(SunAltitude)
wav_SMmap = (0.3:0.02:1.2); %The initial wavelength range

%The shape of the spectral irradiance is calculated. Both of these spectral
%irradiances are normalized, so they can be multiplied by the irradiance
%from the vertex.
if spectra_choice == 1
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSMARTS(wav_SMmap');
    rsd_i_dir = interp1(AM,RSD_i_dir',air_mass);
    rsd_i_dif = interp1(AM,RSD_i_dif',air_mass);
elseif spectra_choice == 2
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSBDarts(wav_SMmap');
    rsd_i_dir = interp1(AM,squeeze(RSD_i_dir(:,:,skytype))',air_mass);
    rsd_i_dif = interp1(AM,squeeze(RSD_i_dif(:,:,skytype))',air_mass);
end

%Calculate the irradiance from each vertex. The vertex of the sun will
%receive the direct spectral distribution, the other ones receive the
%diffuse spectral distribution.
Bi = skymap.*Hemisphere_data(:,3)*rsd_i_dif;
Bi(ind_sun,:) = skymap(ind_sun).*Hemisphere_data(ind_sun,3)*rsd_i_dir;

%Extend the SM to the desired wavelenghts
N_rep_before = find(wav_WEATHER == wav_SMmap(1))-1;
N_rep_after = length(wav_WEATHER) - find(wav_WEATHER == wav_SMmap(end));

SM_compact = [repelem(SM_compact(:,1),1,N_rep_before),SM_compact,repelem(SM_compact(:,end),1,N_rep_after)];


Spec_Irr = sum(Bi.*SM_compact)';

figure;
plot(wav_WEATHER,Spec_Irr);
xlabel('Wavelength [um]')
ylabel('Spectral Irradiance [W/m^2/um]')


%% Functions for plotting

function fig_handle = flatplot3(vertices,faces,facet_color,fig_handle)
%FLATPLOT3 Plot the icohemisphere in 2D
%
% The icohemisphere is divided intro triangles, whose vertices and faces
% are given as an input. The triangles are filled depending on the value of
% the facet color, creating sensitivity maps.
%
% Parameters
% ----------
% vertices : double
%   Vertices of the triangles
% faces : double
%   Position of the face of each triangle
% facet_color : double
%   Value to be filled by the triangles (e.g. irradiance)
% fig_handle : double/matlab.graphics.primitive.Patch
%   Indicates if the plot has already been generated or stores the patch
%   graphic
%
% Returns
% -------
% fig_handle : matlab.graphics.primitive.Patch
%   Stores the generated patch graphic
%
% Developed by R. Santbergen (2017). Commented by A. Alcaniz

% If handle is a number create the figure. Otherwise, just update the color
if isnumeric(fig_handle)
    % Convert 3D to 2D coordinates
    Vcyl = cart2cyl(vertices);
    figure(fig_handle);
    clf
    fig_handle = patch('Vertices',Vcyl,'Faces',faces,...
        'FaceVertexCData',facet_color,'FaceColor','flat');
    axis equal off
    shading flat
    colormap(parula(512))
    caxis([0,1000])
    hc = colorbar;
    set(get(hc,'Title'),'string','Sensitivity [-]')
    text(  0, 95,'North','HorizontalAlignment','center')
    text( 95,  0,'East' ,'HorizontalAlignment','center','Rotation',-90)
    text(  0,-95,'South','HorizontalAlignment','center')
    text(-95,  0,'West' ,'HorizontalAlignment','center','Rotation',90)
else
    set(fig_handle,'FaceVertexCData',facet_color);
end
drawnow

end


function Vcyl = cart2cyl(Vcart)
%CART2CYL Convert 3D cartesian to 2D cylinder coordinates
%
% Convert icohemisphere vertex coordinates (not the light source, which
% is at the center of each vertex)
%
% Parameters
% ----------
% Vcart : double
%   Vector with cartesian coordinates
%
% Returns
% -------
% Vcyl : double
%   Vector with cylinder coordinates

zenith = atand(sqrt(Vcart(:,1).^2 + Vcart(:,2).^2)./Vcart(:,3));
azimuth = atan2d(-Vcart(:,1),-Vcart(:,2));

Vcyl_x = zenith .* -sind(azimuth);
Vcyl_y = zenith .* -cosd(azimuth);
Vcyl_z = zeros(size(Vcyl_x));
Vcyl = [Vcyl_x,Vcyl_y,Vcyl_z];
end
