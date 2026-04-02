function [Flux_angles,Irradiance_angles] = incoming_irradiance_BIF_TAN(Incoming_Irr_input,TOOLBOX_input,CELL_output,SpecData,weather_data,CONSTANTS)
%incoming_irradiance_BIF_TAN Calculates the incoming irradiance for different angles for bifacial tandems.
%
% This function calculates the incoming irradiance a certain hour in the
% year for bifacial tandems
%
% Parameters
% ----------
% Incoming_Irr_input : struc
%   The needed inputs for calculating the incoming irradiance.
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Simulation results of the CELL module
% SpecData : double
%   The spectral irradiance for different air masses
% weather_data : double
%   Weather data from meteonorm
% CONSTANTS : struct
%   Physical constants
%
% Returns
% -------
% Flux_angles : double
%   The flux density for each angle per wavelength [#/(m^2 nm)]
% Irradiance_angles : double
%   The irradiance for each angle per wavelength [W/(m^2 nm)]
% factor : double
%   The correction factor that is needed to match the incoming irradiance
%   with the simulated photo-generated current.
%
% Developed by Y. Blom

wav = Incoming_Irr_input.wav;
Mod_alti = Incoming_Irr_input.Mod_alti;
Mod_azi = Incoming_Irr_input.Mod_azi;
J_abs1 = Incoming_Irr_input.J_abs(1,1,1);
J_abs2 = Incoming_Irr_input.J_abs(1,1,2);
A1 = Incoming_Irr_input.A1;
A2 = Incoming_Irr_input.A2;
AZA = Incoming_Irr_input.AZA;
index = Incoming_Irr_input.index;
extra_sol_power = Incoming_Irr_input.extra_sol_power;
SM_f = Incoming_Irr_input.SM_f;
SM_r = Incoming_Irr_input.SM_r;

sun_alti = weather_data(index,6);
wav_GENPRO = CELL_output.CELL_FRONT.wav;
angles_GENPRO = CELL_output.CELL_FRONT.aoi;


%The spectral distribution is determined by looking at the air mass and
%interpolating

AM = SpecData.AM;
air_mass = 1./sind(sun_alti);
air_mass = min(air_mass,AM(end));
RSD_f_dir = SpecData.RSD_f_dir*1e6;
RSD_f_dif = SpecData.RSD_f_dif*1e6;
RSD_i_dir = SpecData.RSD_i_dir*1e6;
RSD_i_dif = SpecData.RSD_i_dif*1e6;
spectra_choice = SpecData.spectra_choice;

% The perez model is used to calcualte the incoming irradiance
[skymap,ind_sun,skytype] = perez_model(weather_data(index,5), weather_data(index,6),...
    weather_data(index,7), weather_data(index,8),...
    AZA(:,1), AZA(:,2), AZA(:,3),extra_sol_power(index));


if spectra_choice == 1
    rsd_i_dir = interp1(AM,RSD_i_dir',air_mass);
    rsd_f_dir = interp1(AM,RSD_f_dir',air_mass);
    
    rsd_i_dif = interp1(AM,RSD_i_dif',air_mass);
    rsd_f_dif = interp1(AM,RSD_f_dif',air_mass);
    
elseif spectra_choice == 2
    rsd_i_dir = interp1(AM,squeeze(RSD_i_dir(:,:,skytype))',air_mass);
    rsd_f_dir = interp1(AM,squeeze(RSD_f_dir(:,:,skytype))',air_mass);
    
    rsd_i_dif = interp1(AM,squeeze(RSD_i_dif(:,:,skytype))',air_mass);
    rsd_f_dif = interp1(AM,squeeze(RSD_f_dif(:,:,skytype))',air_mass);
    
end

% Sky spectral brightness matrix
Bi = skymap.*AZA(:,3)*rsd_i_dif;
Bf = skymap.*AZA(:,3)*rsd_f_dif;
Bi(ind_sun,:) = skymap(ind_sun).*AZA(ind_sun,3)*rsd_i_dir;
Bf(ind_sun,:) = skymap(ind_sun).*AZA(ind_sun,3)*rsd_f_dir;

% The irradiance and fluxes for each angles are initialised
Irradiance_angles = zeros(length(angles_GENPRO),length(wav),2);
Flux_angles = zeros(length(angles_GENPRO),length(wav),2);

Mod_alti_rear = -Mod_alti;
Mod_azi_rear = Mod_azi+180;

% module settings
d = TOOLBOX_input.Scene.module_mounting.ModMountHeight/100;
Mod_tilt = TOOLBOX_input.Scene.module_mounting.ModTilt;
CR = TOOLBOX_input.Scene.module_mounting.CellRows;          %Number of cell rows
CC = TOOLBOX_input.Scene.module_mounting.CellColumns;       %Number of cell columns
CS = TOOLBOX_input.Scene.module_mounting.CellSpacing/100;   %Cell spacing [m]
ES = TOOLBOX_input.Scene.module_mounting.EdgeSpacing/100;   %Edge spacing [m]
CL = TOOLBOX_input.Scene.module_mounting.CellLength/100;    %Cell length [m]
CW = TOOLBOX_input.Scene.module_mounting.CellWidth/100;     %Cell width [m]
L_mod = CR*CL + (CR-1)*CS + 2*ES;                           %Module length [m] cell + intercel spacing + edge spacing
W_mod = CC*CW + (CC-1)*CS + 2*ES;                           %Module width [m]


%% For each vertex, the contribtion to the incoming irradiance is calculated
%The angle closest to the AOI is selected. The irradiance/flux for this
% vertix is added for this angle.
for i = 1:160
    Sens = squeeze(mean(SM_f(i,:,end,:)));
    vertex_azi = AZA(i,1);
    vertex_alti = 90-AZA(i,2);
    cosAOI = cosd(Mod_alti)*cosd(vertex_alti)*cosd(Mod_azi-vertex_azi)+sind(Mod_alti)*sind(vertex_alti);
    if cosAOI > 0
        AOI = acosd(cosAOI);
        [~,angle_index] = min(abs(AOI-angles_GENPRO));
        Irradiance_angles(angle_index,:,1) = Irradiance_angles(angle_index,:,1)+Bi(i,:).*Sens';
        Flux_angles(angle_index,:,1) = Flux_angles(angle_index,:,1)+Bf(i,:).*Sens';
    end
    
    
    Sens = squeeze(mean(SM_r(i,:,end,:)));
    cosAOI_rear = cosd(Mod_alti_rear)*cosd(vertex_alti)*cosd(Mod_azi_rear-vertex_azi)+sind(Mod_alti_rear)*sind(vertex_alti);
    cosAOI_rear2 = cosd(Mod_alti_rear)*cosd(-vertex_alti)*cosd(Mod_azi_rear-vertex_azi)+sind(Mod_alti_rear)*sind(-vertex_alti);
    angle_Albedo = angle_Albedo_finder(d,L_mod,W_mod,Mod_tilt,vertex_azi);
    %Determine the weight to both angles
    [weigth1, weigth2] = Define_Weights(cosAOI_rear,cosAOI_rear2,Sens,vertex_alti,angle_Albedo);
    if cosAOI_rear > 0
        AOI_rear = acosd(cosAOI_rear);
        %input power from direct back side
        [~,angle_index] = min(abs(AOI_rear-angles_GENPRO));
        Irradiance_angles(angle_index,:,2) = Irradiance_angles(angle_index,:,2)+Bi(i,:).*Sens'*weigth1;
        Flux_angles(angle_index,:,2) = Flux_angles(angle_index,:,2)+Bf(i,:).*Sens'*weigth1;
    end
    %Input power from albedo back side
    if vertex_alti < angle_Albedo
        if cosAOI_rear2 > 0
            AOI_rear_Albedo2 = acosd(cosAOI_rear2);
            [~,angle_index] = min(abs(AOI_rear_Albedo2-angles_GENPRO));
            Irradiance_angles(angle_index,:,2) = Irradiance_angles(angle_index,:,2)+Bi(i,:).*Sens'*weigth2;
            Flux_angles(angle_index,:,2) = Flux_angles(angle_index,:,2)+Bf(i,:).*Sens'*weigth2;
        end
    end
end

%% The correction factor is found based on the absorbed current density
A1_fr = max(interp1(wav_GENPRO*1e-6,A1(:,:,1)',wav),0);
A1_rr = max(interp1(wav_GENPRO*1e-6,A1(:,:,2)',wav),0);
A2_fr = max(interp1(wav_GENPRO*1e-6,A2(:,:,1)',wav),0);
A2_rr = max(interp1(wav_GENPRO*1e-6,A2(:,:,2)',wav),0);
J_test1 = sum(trapz(wav,A1_fr.*Flux_angles(:,:,1)'+A1_rr.*Flux_angles(:,:,2)'));
factor1 = J_test1/J_abs1;
Irradiance_angles = Irradiance_angles/factor1;
Flux_angles = Flux_angles/factor1;
J_test2 = sum(trapz(wav,A2_fr.*Flux_angles(:,:,1)'+A2_rr.*Flux_angles(:,:,2)'));
J_fr2 = sum(trapz(wav,A2_fr.*Flux_angles(:,:,1)'));
factor2 = (J_test2-J_fr2)/(J_abs2-J_fr2);
error = J_test2./J_abs2;
Irradiance_angles(:,:,2) = Irradiance_angles(:,:,2)/factor2;
Flux_angles(:,:,2) = Flux_angles(:,:,2)/factor2;

end

function [weigth1, weigth2] = Define_Weights(cosAOI_rear,cosAOI_rear2,Sens,vertex_alti,angle_Albedo)
    %Define_Weights calculates the weigth for the two different angles.
%
% This function calculates what weight should be given to the two angles
% when the incoming irradiance is calculated.
%
% Parameters
% ----------
% cosAOI_rear : double
%   The angle for the first path (directly).
% cosAOI_rear2 : double
%   The angle for the second path (via ground).
% Sens : double
%   The sensitivity of the vertext as calculated by the MODULE part.
% vertex_alti : double
%   The altitude for the specific vertex.
% angle_Albedo : double
%   The angle for which albedo reflection is still possible.
%
% Returns
% -------
% weigth1 : double
%   The weight for the first path (directly).
% weigth2 : double
%   The weight for the second path (via ground).
%
% Developed by Y. Blom
    if cosAOI_rear > 0 && (cosAOI_rear2 < 0 || vertex_alti > angle_Albedo)
        weigth1 = 1;
        weigth2 = 0;
    elseif cosAOI_rear < 0 && cosAOI_rear2 > 0
        weigth1 = 0;
        weigth2 = 1;
