function [WEATHER_output, TOOLBOX_input] = WEATHER_main(TOOLBOX_input, MODULE_output)
%WEATHER_MAIN Main file for the Weather module in the PVMD toolbox
%
% This function calculates the total absorption and implied photocurrent of
% the system from the Module module considering the weather parameters from
% the location selected and the time period selected by the user
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% MODULE_output : struct
%   Output of the module block
%
% Returns
% -------
% WEATHER_output : struct
%   Output of this weather block
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by unknown (E. Garcia?) Adjusted by Y. Blom. Commented by A. Alcaniz

%---- Modify the toolbox input
if ~TOOLBOX_input.script
    TOOLBOX_input = set_config_weather(TOOLBOX_input);
end
absorptance_all = TOOLBOX_input.deviceOptic.exportAbsorptanceAll;

%---- Integrate the spectral distributions over wavelength
spectra_choice = TOOLBOX_input.irradiation.spectra_choice;
if spectra_choice == 1
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSMARTS(MODULE_output.wav);
elseif spectra_choice == 2
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSBDarts(MODULE_output.wav);
end

%---- Read weather data
[weather_data,number_hours] = load_meteonorm_data(TOOLBOX_input);

%---- Modify DNI and DHI with SVF and shading according to the horizon
if TOOLBOX_input.irradiation.include_horizon_reconstruction
    weather_data(:,7:8) = modify_irradiance(TOOLBOX_input, weather_data(:,5:8));
end

%---- Identify how much the wavelength range should be extended
wav_MODULE = MODULE_output.wav;
N_rep_before = find(wav_WEATHER == wav_MODULE(1))-1;
N_rep_after = length(wav_WEATHER) - find(wav_WEATHER == wav_MODULE(end));

%---- Take sensitivity (SM) information, adding rear side if module is bifacial
% Assume first/final value also represents data before/after wavelength range
if TOOLBOX_input.runPeriodic
    N_submodules = length(MODULE_output.SM_f);
    A = cell(1,N_submodules);
    J = cell(1,N_submodules);
    Irr = cell(1,N_submodules);
    UV = cell(1,N_submodules);
    for Submod_i = 1:N_submodules
        SM = MODULE_output.SM_f{Submod_i};
        if isfield(MODULE_output,'SM_r'), SM = SM + MODULE_output.SM_r{Submod_i}; end

        %---- Extend the SM to the wavelength range of MODULE to the wavelength
        %range of weather
        SM = cat(4,repelem(SM(:,:,:,1),1,1,1,N_rep_before),SM,repelem(SM(:,:,:,end),1,1,1,N_rep_after));

        %---- Compute the absorption and photocurrent at cell level
        [A{Submod_i}, J{Submod_i}, Irr{Submod_i}, UV{Submod_i},specIrr] = compute_absorption_and_photocurrent(...
            MODULE_output.skydome.Vs,MODULE_output.skydome.Fs,MODULE_output.skydome.AZA,...
            SM,number_hours,weather_data,RSD_i_dir, RSD_f_dir, RSD_i_dif, RSD_f_dif,...
            AM,TOOLBOX_input.irradiation.plot_weather,wav_WEATHER,spectra_choice,absorptance_all);
    end
else
    N_panels = TOOLBOX_input.Scene.N_panels;
    A = cell(N_panels, 1);
    J = cell(N_panels, 1);
    Irr = cell(N_panels, 1);
    UV = cell(N_panels, 1);
    for Panel_i = 1:N_panels
        SM = MODULE_output.SM_f{Panel_i};
        if isfield(MODULE_output,'SM_r'), SM = SM + MODULE_output.SM_r{Panel_i}; end
        
        %---- Extend the SM to the wavelength range of MODULE to the wavelength
        %range of weather
        SM = cat(4,repelem(SM(:,:,:,1),1,1,1,N_rep_before),SM,repelem(SM(:,:,:,end),1,1,1,N_rep_after));
        
        %---- Compute the absorption and photocurrent at cell level
        [A{Panel_i}, J{Panel_i}, Irr{Panel_i}, UV{Panel_i},specIrr] = compute_absorption_and_photocurrent(...
            MODULE_output.skydome.Vs,MODULE_output.skydome.Fs,MODULE_output.skydome.AZA,...
            SM,number_hours,weather_data,RSD_i_dir, RSD_f_dir, RSD_i_dif, RSD_f_dif,...
            AM,TOOLBOX_input.irradiation.plot_weather,wav_WEATHER,spectra_choice,absorptance_all);
    end
    
end
% Store the parameters that will be used later on
% The irradiance that is saved (GHI) is not modified according to the
% surroundings. However that irradiance is only employed for plotting, so
% it is not an issue
WEATHER_output.A = A;
WEATHER_output.J = J;
WEATHER_output.Irr = Irr;
WEATHER_output.specIrr = specIrr;
WEATHER_output.wav = wav_WEATHER;
WEATHER_output.Period = [...
    TOOLBOX_input.irradiation.init_day,...
    TOOLBOX_input.irradiation.init_month,...
    TOOLBOX_input.irradiation.end_day,...
    TOOLBOX_input.irradiation.end_month];
WEATHER_output.ambient_temperature = weather_data(:,9);
WEATHER_output.wind_speed = weather_data(:,10);
WEATHER_output.month = weather_data(:,2);
WEATHER_output.day = weather_data(:,3);
WEATHER_output.UV = UV;
if size(weather_data,2) > 11 %only for those files with this data
    WEATHER_output.RH = weather_data(:,12);
    WEATHER_output.UVa = weather_data(:,13);
    WEATHER_output.UVb = weather_data(:,14);
end

if ~TOOLBOX_input.script
    disp('Irradiance calculation finished.')
end
end


function [A, J, Irr, UV,specIrr] = compute_absorption_and_photocurrent(Vs, Fs, skydome,...
    SM, number_hours, weather_data, RSD_i_dir, RSD_f_dir, RSD_i_dif, RSD_f_dif, AM, plot_skymap,wav,spectra_choice,absorptance_all)
%COMPUTE_ABSORPTION_AND_PHOTOCURRENT Calculate the absorption and
%photocurrent at cell level
%
% Parameters
% ----------
% Vs : double
%   Vertices of the triangles in the skydome
% Fs : double
%   Faces of the triangles in the skydome
% skydome : double
%   Azimuth [deg], zenith [deg] and area [sr] of every triangle in the skydome
% SM : double
%   Sensitivity of each cell
% number_hours : double
%   Length of the weather_data variable
% weather_data : double
%   Weather data for the selected period
% DNI_mod : double
%   DNI corrected by the shading factor
% DHI_mod : double
%   DHI corrected by the SVF
% RSD_i_dir : double
%   Relative spectral distribution in terms of irradiance for direct irradiance[W/m2]
% RSD_f_dir : double
%   Relative spectral distribution in terms of photon flux for direct irraiance [W/m2]
% RSD_i_dif : double
%   Relative spectral distribution in terms of irradiance for diffuse irradiance[W/m2]
% RSD_f_dif : double
%   Relative spectral distribution in terms of photon flux for diffuse irraiance [W/m2]
% AM : double
%   Air mass numbers for which RSD is calculated
% plot_skymap: double
%   Indicator if the skymap should be plotted
% wav: double
%   The wavelength spectrum of the spectral distributions
% spectra_choice: double
%   The choice of the spectra that should be used (1 is SMARTS, 2 is
%   SBDarts).
%
% Returns
% -------
% A : double
%   Total absorption in every module-cell [W/m2]
% J : double
%   Implied photocurrent in every every absorber layer [mA/cm2]
% Irr : double
%   received irradiance in every cell [W/m2]

%---- Variables relocation
number_of_cells = size(SM,2);
number_of_absorber_layers = size(SM,3)-2;
if absorptance_all
    A = zeros(number_hours,number_of_cells,number_of_absorber_layers);
else
    A = zeros(number_hours,number_of_cells);
end
J = zeros(number_hours,number_of_cells,number_of_absorber_layers);
Irr = zeros(number_hours,number_of_cells);
UV = zeros(number_hours,number_of_cells);
figure_handle = 9;

% ---- Compute the air mass. Limit the maximum so that extrapolation is not
% needed
air_mass = 1./sind(weather_data(:,6));
air_mass(air_mass>AM(end)) = AM(end);

% ---- Calculate the extraterrestrial solar radiation
% Compute the day of the year for the Perez model
day_year = pvl_date2doy(weather_data(:,1), weather_data(:,2), ...
    weather_data(:,3));
% Expression used for the calculation is eq. 1.4.1b in "Solar engineering
% of thermal processes" by J. Duffie et al.
day_angle = 2*pi*(day_year - 1)/365.25;
E0 = 1.00011 + 0.034221*cos(day_angle) + 0.00128*sin(day_angle) + ...
    0.000719*cos(2*day_angle) + 0.000077*sin(2*day_angle);
load('constants/weather_params.mat', 'AM0')
extra_sol_power = E0*AM0;

%The wavelengths untill 400 nm are considered to be part of UV light.
wav_UV_ind = wav < 0.400;

specIrr = zeros(1,length(wav));

%---- Time-dependent calculation of absorption and photocurrent
for t = 1:number_hours
    if weather_data(t,8)>0
        % Create skymap [W/m2 SR]
        [skymap,ind_sun,skytype] = perez_model(weather_data(t,5), weather_data(t,6),...
            weather_data(t,7), weather_data(t,8),...
            skydome(:,1), skydome(:,2), skydome(:,3), extra_sol_power(t));
        % Get the spectral distribution over the interval. Interpolate the
        % spectral distribution depending on the air mass value
        if spectra_choice == 1
            rsd_i_dir = interp1(AM,RSD_i_dir',air_mass(t));
            rsd_f_dir = interp1(AM,RSD_f_dir',air_mass(t));
            
            rsd_i_dif = interp1(AM,RSD_i_dif',air_mass(t));
            rsd_f_dif = interp1(AM,RSD_f_dif',air_mass(t));
            
        elseif spectra_choice == 2
            rsd_i_dir = interp1(AM,squeeze(RSD_i_dir(:,:,skytype))',air_mass(t));
            rsd_f_dir = interp1(AM,squeeze(RSD_f_dir(:,:,skytype))',air_mass(t));
            
            rsd_i_dif = interp1(AM,squeeze(RSD_i_dif(:,:,skytype))',air_mass(t));
            rsd_f_dif = interp1(AM,squeeze(RSD_f_dif(:,:,skytype))',air_mass(t));
            
        end
        % Plot skymap if desired
        if plot_skymap
            figure_handle = flatplot3(Vs,Fs,skymap,figure_handle);
            title(mod(t,24))
        end
        
        % Sky spectral brightness matrix
        Bi = skymap.*skydome(:,3)*rsd_i_dif;
        Bf = skymap.*skydome(:,3)*rsd_f_dif;
        Bi(ind_sun,:) = skymap(ind_sun).*skydome(ind_sun,3)*rsd_i_dir;
        Bf(ind_sun,:) = skymap(ind_sun).*skydome(ind_sun,3)*rsd_f_dir;
        
        % Combine sensitivity and sky map. Calculate the total absorbed
        % power (for thermal model) in every cell and the absorbed photon
        % flux in every absorber layer
        for cell = 1:number_of_cells
            if absorptance_all
                for layer = 1:number_of_absorber_layers
                    A(t,cell,layer) = trapz(wav,sum(Bi.*squeeze(SM(:,cell,layer+1,:)))');
                end
            else
                A(t,cell) = trapz(wav,sum(Bi.*squeeze(SM(:,cell,1,:)))');
            end
            for layer = 1:number_of_absorber_layers
                J(t,cell,layer) = trapz(wav,sum(Bf.*squeeze(SM(:,cell,layer+1,:)))');
            end
            Irr(t,cell) = trapz(wav,sum(Bi.*squeeze(SM(:,cell,end,:)))');
            UV(t,cell) = trapz(wav(wav_UV_ind),sum(Bi(:,wav_UV_ind).*squeeze(SM(:,cell,end,wav_UV_ind)))');
        end
        specIrr = specIrr + sum(Bi.*squeeze(mean(SM(:,:,end,:),2)));
    end
end

end


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


