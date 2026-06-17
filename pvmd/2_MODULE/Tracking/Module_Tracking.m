function [WEATHER_output, MODULE_output] = Module_Tracking(TOOLBOX_input, CELL_output)
%WEATHER_periodic Calculates the sensitivity map for periodic simulations
%
% This function calculates the sensitivity map of the PV system
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Output of the CELL block
%
% Returns
% -------
% MODULE_output : struct
%   Output of this module block
%
% Developed by Orestis Chatzilampos, integrated by Youri Blom, Tibo will
% improve

[Vs,Fs,A,~, zenith_vertex,azim_vertex] = icohemisphere(2,1);
%icosphere is the function that creates vertices andcleatc all of their properties
%rl=2 because vertices=160 I care only about these outputs not all MODULE_output
%Vs=the points cordinates, Fs the triplet of points forming the triangle

%%
%Keep the same format as Module_output as before and pass that on to weather.
% I create module output using only the function I need, not the ray traycing.
% Call only the functions needed to create Module_output_structure
%this the function that builds the module frame, properties and scheme

if strcmp(TOOLBOX_input.Scene.trackingType,'HSATS')
    [Vf,Ff,~,Bf,Ncells,Acell,Amod,ML,MW] = moduleGeometryTracking_hsats(TOOLBOX_input,CELL_output,1,1);
elseif strcmp(TOOLBOX_input.Scene.trackingType,'TSATS')
    [Vf,Ff,~,Bf,Ncells,Acell,Amod,ML,MW] = moduleGeometryTracking_tsats(TOOLBOX_input,CELL_output,1,1);
elseif strcmp(TOOLBOX_input.Scene.trackingType,'DATS')
    [Vf,Ff,~,Bf,Ncells,Acell,Amod,ML,MW] = moduleGeometryTracking_dats(TOOLBOX_input,CELL_output,1,1);
end

%create Module_output similar to the original
MODULE_output.skydome.AZA=[azim_vertex zenith_vertex A];
MODULE_output.skydome.Vs = Vs;                      %skydome vertices (for plotting SKYmap) each row represents x,y,z coord
MODULE_output.skydome.Fs = Fs;                      %triangular shapes each row represent 3 vertices
MODULE_output.wav = CELL_output.CELL_FRONT.wav;     %pass on wavelength information 
MODULE_output.N=Ncells;
MODULE_output.A=Acell;
MODULE_output.Amod=Amod;
MODULE_output.ML=ML*1e-2;
MODULE_output.MW=MW*1e-2;

%% Add iCellContainingPlane to Toolbox_input as this should only be calculated once
% the faces will not change only the vertices will
N_faces = size(Ff,1);%number of triangles in the scene
cellCorners = reshape(Vf(1:Ncells*4,:)',3,4,Ncells);
coor_faces = single(reshape(Vf(Ff',:)',3,3,N_faces));
cell_inplane = zeros(Ncells,N_faces);

for iCell=1:Ncells
    iCellContainingPlane = nan(1,N_faces,'single');
    for iTriangle = 1:N_faces
        for iTriCorner = 1:3
            if sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,1,iCell)).^2)<1 || ...
                    sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,2,iCell)).^2)<1 || ...
                    sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,3,iCell)).^2)<1 || ...
                    sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,4,iCell)).^2)<1
                iCellContainingPlane(iTriangle) = iTriangle;
            end
        end
    end
    cell_inplane(iCell,:) = iCellContainingPlane;
end

TOOLBOX_input.Scene.module_mounting.cell_inplane_front=cell_inplane;

if Bf == 1
    cell_inplane = zeros(Ncells,N_faces);
    cellCorners = reshape(Vf(Ncells*4+1:2*Ncells*4,:)',3,4,Ncells);
    for iCell=1:Ncells
        iCellContainingPlane = nan(1,N_faces,'single');
        for iTriangle = 1:N_faces
            for iTriCorner = 1:3
                if sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,1,iCell)).^2)<1 || ...
                        sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,2,iCell)).^2)<1 || ...
                        sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,3,iCell)).^2)<1 || ...
                        sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,4,iCell)).^2)<1
                    iCellContainingPlane(iTriangle) = iTriangle;
                end
            end
        end
        cell_inplane(iCell,:) = iCellContainingPlane;
    end

    TOOLBOX_input.Scene.module_mounting.cell_inplane_rear=cell_inplane;
end


%% Initialize matrix where the sensitivity maps will be stored

az = 0:359;
tilt = 0:90;

SM_stored = cell(numel(az), numel(tilt));
        
%% Adjust the period tested & Location

%---- Integrate the spectral distributions over wavelength
spectra_choice = TOOLBOX_input.irradiation.spectra_choice;
plot_skymap = TOOLBOX_input.irradiation.plot_weather;
figure_handle = 9;
if spectra_choice == 1
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSMARTS(MODULE_output.wav,TOOLBOX_input.constants);         % relative spectral distribution
elseif spectra_choice == 2
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSBDart(MODULE_output.wav,TOOLBOX_input.constants);
end

%---- Read weather data
[weather_data,number_hours] = load_meteonorm_data(TOOLBOX_input);

%---- Modify DNI and DHI with SVF and shading according to the horizon % in
%my case this is defined false from toolbox input loaded before
if TOOLBOX_input.irradiation.include_horizon_reconstruction
    weather_data(:,7:8) = modify_irradiance(TOOLBOX_input, weather_data(:,5:8));
end

%% Compute_absorption_and_photocurrent edited function.
% No need to run the weather_main. This is an editted version

%---- Variables reassigned
skydome=MODULE_output.skydome.AZA;

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

layers=length(CELL_output.CELL_FRONT.lay)-2;  % SOS define them before or make parametric

%% Algorith selction choice
% The code for the two methods and the switch criterion. First we find the
% optimal orientation and then then calculate Jph= Jsc for each layer!

all_opt_azimuths = TOOLBOX_input.Scene.module_mounting.ModAzimuth*ones(number_hours,1); %initial case
all_opt_tilts = TOOLBOX_input.Scene.module_mounting.ModTilt*ones(number_hours,1); %intitial case

%individual cell sensitivty
A = zeros(number_hours, Ncells); % Absorption in W/m2
J=zeros(number_hours, Ncells,layers); % only the absorber layers, Jph in mA/cm2
Irr = zeros(number_hours,Ncells);% Incoming irradiance in W/m2

%trackers data
P_tracker= TOOLBOX_input.Scene.module_mounting.P_tracker; % in W
eta_tracker=TOOLBOX_input.Scene.module_mounting.eta_tracker;
omega_tracker=TOOLBOX_input.Scene.module_mounting.omega_tracker;%in deg/s
E_tracker=zeros(number_hours,3); % col 1: azim consumption, col 2: tilt consumption, col 3:both axis consumption

delta_azim=zeros(number_hours,1); % no need to store them but to have an idea
delta_tilt=zeros(number_hours,1);

for t=1:1:number_hours

    if weather_data(t,8)>0  % get in this loop only for daylight
        
        [skymap,ind_sun,skytype] = perez_model(weather_data(t,5), weather_data(t,6),...
            weather_data(t,7), weather_data(t,8),skydome(:,1), skydome(:,2), skydome(:,3), extra_sol_power(t),TOOLBOX_input.settings);

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

        % Sky spectral brightness matrix
        Bi = skymap.*skydome(:,3)*rsd_i_dif;
        Bf = skymap.*skydome(:,3)*rsd_f_dif;
        Bi(ind_sun,:) = skymap(ind_sun).*skydome(ind_sun,3)*rsd_i_dir;          %index sun, triangle closest to position sun
        Bf(ind_sun,:) = skymap(ind_sun).*skydome(ind_sun,3)*rsd_f_dir;

        if plot_skymap
            figure_handle = flatplot3(Vs,Fs,skymap,figure_handle);
            title(mod(t,24))
        end

        %  SOS THIS IS THE SWITCH CRITERION FOR ALGORITHMS CAN BE CHANGED
        % CROSS METHOD


        % A lot of inputs when calling are needed for the SM_new calculation, happening within the function of the cross method
        %use the initial guess of last_hour for cross_method for somecases might cause faster convergence
        % for some other irrelevant

        %in cross call method below t-1= 0 will never happen as always it starts with 01.00 hour
        % that is always night across the year for Delft. and also in
        % stockholm which is further north, so it want hit for t=1. I
        % think this is same for all locations I choose


        [all_opt_azimuths(t),all_opt_tilts(t),A(t,:),J(t,:,:),Irr(t,:),SM_stored]= CrossMethod(CELL_output, wav_WEATHER,Bi,Bf,all_opt_azimuths(t-1),all_opt_tilts(t-1), SM_stored,TOOLBOX_input,Ncells);

    end
    
    % THE TRACKER CONSUMPTION WORKS!! CHECKED, NO NEED TO SET AZIM,TILT=0 IN
    % THE BEGINNING OF THE DAY IT IS PREALLOCATED. ONLY CHANGE TO AZIM=180
    % PREALLOCATION FOR SOUTH HEM LOCATION !!!!!!!
  
    if t>1 % it gets in for t=2, so I can calculate the difference

        azim_diff = abs(all_opt_azimuths(t) - all_opt_azimuths(t-1)); % always absolute values
        delta_azim(t) =abs(mod(azim_diff + 180, 360) - 180);  % to ensure that delta_azim is always between 0-180
        % circular nature of azim=0-360, delta azim always positive % no need to
        % store it

        delta_tilt(t)=abs(all_opt_tilts(t)-all_opt_tilts(t-1)); % no need for smth tilt is always between 0-90, no circular nature.

        E_tracker(t,1)=P_tracker*eta_tracker*delta_azim(t)/(omega_tracker*3600); % in Wh  consumption in azim axis
        E_tracker(t,2)=P_tracker*eta_tracker*delta_tilt(t)/(omega_tracker*3600); % in Wh  consumption in tilt axix

    end
end

E_tracker(:,3)=E_tracker(:,1)+E_tracker(:,2);   % both axis tracker consump per hour Wh


% delta_tilt=diff(all_opt_tilts);  % no need for these keep in mind the diff function
% delta_azim=diff(all_opt_azimuths);

% no need for alba's corrections they are fixed in the perez_model_script
% the rest of Weather_main refer to the skymap plots and are not needed !!
%% Redefine the Weather_output to match the orginal
% store mod.tilt_optimums

MODULE_output.ModTilt = all_opt_tilts; % to be used in thermal
MODULE_output.ModAz  = all_opt_azimuths;
MODULE_output.SM_f = SM_stored(TOOLBOX_input.Scene.module_mounting.ModAzimuth+1,TOOLBOX_input.Scene.module_mounting.ModTilt+1);
MODULE_output.E_tracker = E_tracker;
% 

WEATHER_output.J={J};

WEATHER_output.A = {A};
WEATHER_output.Irr = {Irr};
WEATHER_output.Period = [...
    TOOLBOX_input.irradiation.init_day,...
    TOOLBOX_input.irradiation.init_month,...
    TOOLBOX_input.irradiation.end_day,...
    TOOLBOX_input.irradiation.end_month];
WEATHER_output.ambient_temperature = weather_data(:,9);
WEATHER_output.wind_speed = weather_data(:,10);
WEATHER_output.month = weather_data(:,2);
WEATHER_output.day = weather_data(:,3);

if size(weather_data,2) > 11 %only for those files with these data
    WEATHER_output.RH = weather_data(:,12);
    WEATHER_output.UVa = weather_data(:,13);
    WEATHER_output.UVb = weather_data(:,14);
end

end