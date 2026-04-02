function [WEATHER_output, MODULE_output, TOOLBOX_input] = TwoAxisTracking(TOOLBOX_input, CELL_output)
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
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by Orestis Chatzilampos, integrated by Youri Blom

[Vs,Fs,A,~, zenith_vertex,azim_vertex] = icohemisphere(2,1);
%icosphere is the function that creates vertices andcleatc all of their properties
%rl=2 because vertices=160 I care only about these outputs not all MODULE_output
%Vs=the points cordinates, Fs the triplet of points forming the triangle


%Keep the same format as Module_output as before and pass that on to weather.
% I create module output using only the function I need, not the ray traycing.
% Call only the functions needed to create Module_output_structure
%this the function that builds the module frame, properties and scheme
[~,~,~,~,Ncells,Acell,Amod,~,~,TOOLBOX_input] = moduleGeometry(TOOLBOX_input,CELL_output,0);

%create Module_output similar to the original
MODULE_output.skydome.AZA=[azim_vertex zenith_vertex A];
MODULE_output.skydome.Vs = Vs;                      %skydome vertices (for plotting SKYmap)
MODULE_output.skydome.Fs = Fs;
MODULE_output.wav = CELL_output.CELL_FRONT.wav;     %pass on wavelength information
MODULE_output.N=Ncells;
MODULE_output.A=Acell;
MODULE_output.Amod=Amod;
%MODULE_output.SM_f=  %missing called later in optimization
%MODULE_output.ModTilt = TOOLBOX_input.Scene.module_mounting.ModTilt; %
%changes all the time, add it after optimization to locate it!

%% Adjust the period tested & Location
% The rest are already defined (overwrite toolbox input)
if ~TOOLBOX_input.script
    TOOLBOX_input = set_config_weather(TOOLBOX_input);
end

%---- Integrate the spectral distributions over wavelength
spectra_choice = TOOLBOX_input.irradiation.spectra_choice;
plot_skymap = TOOLBOX_input.irradiation.plot_weather;
figure_handle = 9;
if spectra_choice == 1
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSMARTS(MODULE_output.wav);
elseif spectra_choice == 2
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,AM,wav_WEATHER] = spectral_distrSBDarts(MODULE_output.wav);
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
wav=wav_WEATHER;

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

%% Algorith selction choice
% The code for the two methods and the switch criterion. First we find the
% optimal orientation and then then calculate Jph= Jsc for each layer!

%Both algorithms input
if TOOLBOX_input.Scene.module_mounting.Albedo_eff
    albedo=TOOLBOX_input.Scene.module_mounting.Albedo;
else
    [~,~,data_folder] = get_folder_structure;
    Ground_material = TOOLBOX_input.Scene.module_mounting.Ground_material;
    lib_path = fullfile(data_folder, 'Material Library','SpectralReflectivityLibrary',Ground_material);
    load(lib_path,'lambda','specRefl');
    wav_al = CELL_output.CELL_FRONT.wav;
    albedo = interp1(lambda,specRefl,wav_al*1e3);
end
layers=length(CELL_output.CELL_FRONT.lay)-2;  % SOS define them before or make parametric

all_opt_azimuths=zeros(number_hours,1); % SOS CHANGE IT PER LOCATION  default south=0 deg, for Northern hem, to be ready for next day
all_opt_tilts=zeros(number_hours,1); % default tilt=zero so that is ready for the next day

A = zeros(number_hours, 1); % Absorption in W/m2
J=zeros(number_hours, layers); % only the absorber layers, Jph in mA/cm2
Irr = zeros(number_hours,1);% Incoming irradiance in W/m2
counter_cross=0; % to check how many times cross method is used
counter_surr=0; % similar

%find all overcast hours %ΜΑΥΒΕ REMOVE THOSE WERE DNO=0, DHI=1 PUT DHI>1 I
%Want during the day not overcast hour on the sunset
% overcast_hours=find(weather_data(:,7)==0 & weather_data(:,8)>6); % not the # hour of the year but that of the selected period
% consecutive_overcast_hours_first_hour=find(diff(overcast_hours)==1); % That is the index number of the consecutive overcast hours in the period given % first of the two number of the pair
% consecutive_overcast_hours_second_hour=consecutive_overcast_hours_first_hour+1 ;% returns the index of the second number (last hour) of the pair
% overcast_hours_second_hour_non_index=overcast_hours(consecutive_overcast_hours_second_hour); % SOS retuns the hours not the indices of the second in the pair

% consecutive_pairs=[consecutive_overcast_hours_first_hour; consecutive_overcast_hours_first_hour+1];
%col 1: first_number, col 2:second_number_ of the pair % the left number has always a consecutive the col 2

%trackers data
P_tracker=25; % in W
eta_tracker=0.97;
omega_tracker=0.3;%in deg/s
E_tracker=zeros(number_hours,3); % col 1: azim consumption, col 2: tilt consumption, col 3:both axis consumption

delta_azim=zeros(number_hours,1); % no need to store them but to have an idea
delta_tilt=zeros(number_hours,1);

for t=1:1:number_hours

    if weather_data(t,8)>0  % get in this loop only for daylight

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

        % returns skymap same for single hour 160x1 the calculate skymap function has the if statement when DHI>0 to return
        % only irradiance values when there is sun.

        %         [skymap,period_tested]=calculate_skymap_V2(TOOLBOX_input.irradiation.init_day,TOOLBOX_input.irradiation.end_day, ...
        %             TOOLBOX_input.irradiation.init_month,TOOLBOX_input.irradiation.end_month,TOOLBOX_input.irradiation.year_choice, ...
        %             weather_data,t);
        [skymap,ind_sun,skytype] = perez_model(weather_data(t,5), weather_data(t,6),...
            weather_data(t,7), weather_data(t,8),...
            skydome(:,1), skydome(:,2), skydome(:,3), extra_sol_power(t));

        % Sky spectral brightness matrix
        Bi = skymap.*skydome(:,3)*rsd_i_dif;
        Bf = skymap.*skydome(:,3)*rsd_f_dif;
        Bi(ind_sun,:) = skymap(ind_sun).*skydome(ind_sun,3)*rsd_i_dir;
        Bf(ind_sun,:) = skymap(ind_sun).*skydome(ind_sun,3)*rsd_f_dir;

        if plot_skymap
            figure_handle = flatplot3(Vs,Fs,skymap,figure_handle);
            title(mod(t,24))
        end

        %  SOS THIS IS THE SWITCH CRITERION FOR ALGORITHMS CAN BE CHANGED
        if      weather_data(t,7)<2*weather_data(t,7)
            %weather_data(t,7)/weather_data(t,8) <= 0.7 && (weather_data(t,7) <= 150 || weather_data(t,8) <= 150) % OLDER SWITCH CRITERION


            % SURROGATE OPTIMIZATION

            counter_surr=counter_surr+1;
            % define the objective
            objective = @(x) -surrogate_fun_per_hour(x,CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
                albedo,Bi);

            % Define the bounds for the optimizer
            lb = [0 0];
            ub = [359 90];

            iterations_limit=150;

            options = optimoptions('surrogateopt', 'MaxFunctionEvaluations', iterations_limit,'PlotFcn', [], 'Display', 'off'); % no plots are displayed and nothing in the command window

            intcon = [1, 2];  % indicates that both x(1)=azim and x(2)=tilt are integer variables

            % Optimize the objective funtion with surrogate optimization
            [x_opt, fval, ~, ~] = surrogateopt(objective, lb, ub, intcon,options);

            % calculate the maximum absorbed irradiance
            max_abs_irradiance = -fval;

            % calculate Jph absorbers in optimum position, because it is not stored

            SM_J_surr=calculate_SM_new_single_orientation(CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
                albedo,x_opt(1),x_opt(2));

            A(t)=trapz(wav,sum(Bi.*squeeze(SM_J_surr(:,1,:)))');
            for lay_i = 1:layers
                J(t,lay_i)=trapz(wav,sum(Bf.*squeeze(SM_J_surr(:,lay_i+1,:)))');
            end
            all_opt_azimuths(t)=x_opt(1);
            all_opt_tilts(t)=x_opt(2);
            Irr(t)=max_abs_irradiance;

        else
            % CROSS METHOD

            counter_cross=counter_cross+1;

            % A lot of inputs when calling are needed for the SM_new calculation, happening within the function of the cross method
            %use the initial guess of last_hour for cross_method for somecases might cause faster convergence
            % for some other irrelevant

            %in cross call method below t-1= 0 will never happen as always it starts with 01.00 hour
            % that is always night across the year for Delft. and also in
            % stockholm which is further north, so it want hit for t=1. I
            % think this is same for all locations I choose


            [all_opt_azimuths(t),all_opt_tilts(t),A(t),J(t),Irr(t)]= cross_method_fun_per_hour(CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
                albedo,Bi,Bf,all_opt_azimuths(t-1),all_opt_tilts(t-1));

        end

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

MODULE_output.ModTilt=all_opt_tilts; % to be used in thermal
MODULE_output.E_tracker = E_tracker;

% add cell dimension to Absorption and Jsc
Irr=repmat(Irr,1,Ncells);
A=repmat(A,1,Ncells);

WEATHER_output.J=zeros(number_hours,Ncells,layers);

for t=1:Ncells
    WEATHER_output.J(:,t,:) = J;
end

WEATHER_output.A = A;
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

if size(weather_data,2) > 11 %only for those files with these data
    WEATHER_output.RH = weather_data(:,12);
    WEATHER_output.UVa = weather_data(:,13);
    WEATHER_output.UVb = weather_data(:,14);
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