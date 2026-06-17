function [WEATHER_output] = WEATHER_main(TOOLBOX_input, MODULE_output)
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
%
% Developed by unknown (E. Garcia?) Adjusted by Y. Blom. Commented by A. Alcaniz

absorptance_all = TOOLBOX_input.deviceOptic.exportAbsorptanceAll;

%---- Integrate the spectral distributions over wavelength
spectra_choice = TOOLBOX_input.irradiation.spectra_choice;
if spectra_choice == 1
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSMARTS(MODULE_output.wav,TOOLBOX_input.constants);
elseif spectra_choice == 2
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSBDart(MODULE_output.wav,TOOLBOX_input.constants);
end

%---- Read weather data
[weather_data,number_hours] = load_meteonorm_data(TOOLBOX_input);

%---- Modify DNI and DHI with SVF and shading according to the horizon
if TOOLBOX_input.irradiation.include_horizon_reconstruction
    skylineMask = makeSkylineMask(TOOLBOX_input,MODULE_output.skydome);
else
    skylineMask = ones(length(MODULE_output.skydome.AZA),1);
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
        SM = SM.*skylineMask;

        %---- Compute the absorption and photocurrent at cell level
        [A{Submod_i}, J{Submod_i}, Irr{Submod_i}, UV{Submod_i}] = calculateAbsorptionPhotocurrent(TOOLBOX_input,...
            MODULE_output.skydome,SM,number_hours,weather_data,RSD_i_dir, RSD_f_dir, RSD_i_dif, RSD_f_dif,...
            AM,wav_WEATHER,spectra_choice,absorptance_all);
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
        SM = SM.*skylineMask;
        
        %---- Compute the absorption and photocurrent at cell level
        [A{Panel_i}, J{Panel_i}, Irr{Panel_i}, UV{Panel_i}] = calculateAbsorptionPhotocurrent(TOOLBOX_input,...
            MODULE_output.skydome,SM,number_hours,weather_data,RSD_i_dir, RSD_f_dir, RSD_i_dif, RSD_f_dif,...
            AM,wav_WEATHER,spectra_choice,absorptance_all);
    end
    
end
% Store the parameters that will be used later on
% The irradiance that is saved (GHI) is not modified according to the
% surroundings. However that irradiance is only employed for plotting, so
% it is not an issue
WEATHER_output.A = A;
WEATHER_output.J = J;
WEATHER_output.Irr = Irr;
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