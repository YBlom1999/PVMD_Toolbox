function [A, J, Irr, UV] = calculateAbsorptionPhotocurrent(TOOLBOX_input,skydome,SM, number_hours, weather_data, RSD_i_dir, RSD_f_dir, RSD_i_dif, RSD_f_dif, AM,wav,spectra_choice,absorptance_all)
%COMPUTE_ABSORPTION_AND_PHOTOCURRENT Calculate the absorption and
%photocurrent at cell level
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% skydome: struct
%   The information (Vertices, faces, AZA) of the skydome
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
%   The choice of the spectra that should be used (1 is SMARTS, 2 is SBDarts).
% settingsFilePath: string
%   The path of the file where the settings are stored
%
% Returns
% -------
% A : double
%   Total absorption in every module-cell [W/m2]
% J : double
%   Implied photocurrent in every every absorber layer [mA/cm2]
% Irr : double
%   received irradiance in every cell [W/m2]

%---- Obtain skydome parameters
Vs = skydome.Vs;
Fs = skydome.Fs;
AZA = skydome.AZA;

%---- Read inputs from TOOLBOX_input
plot_skymap = TOOLBOX_input.irradiation.plot_weather;
settingsFilePath = TOOLBOX_input.settings;
AM0 = TOOLBOX_input.constants.AM0;

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
day_year = pvl_date2doy(weather_data(:,1), weather_data(:,2),weather_data(:,3));
% Expression used for the calculation is eq. 1.4.1b in "Solar engineering
% of thermal processes" by J. Duffie et al.
day_angle = 2*pi*(day_year - 1)/365.25;
E0 = 1.00011 + 0.034221*cos(day_angle) + 0.00128*sin(day_angle) + ...
    0.000719*cos(2*day_angle) + 0.000077*sin(2*day_angle);
extra_sol_power = E0*AM0;

%The wavelengths untill 400 nm are considered to be part of UV light.
wav_UV_ind = wav < 400e-9;

%---- Time-dependent calculation of absorption and photocurrent
for t = 1:number_hours
    if weather_data(t,8)>0
        % Create skymap [W/m2 SR]
        [skymap,ind_sun,skytype] = perez_model(weather_data(t,5), weather_data(t,6),...
            weather_data(t,7), weather_data(t,8),AZA(:,1), AZA(:,2), AZA(:,3), extra_sol_power(t),settingsFilePath);
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
            figure_handle = flatplot3(Vs,Fs,skymap,[0,1000],parula(512),'Irradiance [W/m^2]',figure_handle);
            title(mod(t,24))
        end
        
        % Sky spectral brightness matrix
        Bi = skymap.*AZA(:,3)*rsd_i_dif;
        Bf = skymap.*AZA(:,3)*rsd_f_dif;
        Bi(ind_sun,:) = skymap(ind_sun).*AZA(ind_sun,3)*rsd_i_dir;
        Bf(ind_sun,:) = skymap(ind_sun).*AZA(ind_sun,3)*rsd_f_dir;
        
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
    end
end

end