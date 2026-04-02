function Losses_Operating = Analysis_Operating(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_output)
%Analysis_Operating Calculates the losses at operating conditions
%
% This function calculates the losses that are present in the PV system for
% operating conditions
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Simulation results of the CELL module
% MODULE_output : struct
%   Simulation results of the MODULE module
% WEATHER_output : struct
%   Simulation results of the WEATHER module
% THERMAL_output : struct
%   Simulation results of the THERMAL module
% ELECTRIC_output : struct
%   Simulation results of the ELECTRIC module
% CONVERSION_output : struct
%   Simulation results of the CONVERSION module
%
% Returns
% -------
% Losses_operating : struct
%   Simulation results of the loss analysis at operating conditions
%
% Developed by Y. Blom



%% constants
CONSTANTS.h = 6.62607004e-34;
CONSTANTS.q = 1.60217662e-19;
CONSTANTS.c = 299792458;
CONSTANTS.k = 1.380649e-23;
CONSTANTS.T_S = 5778;


%% Loaded input from the structures
DCP = ELECTRIC_output.DCP;
type = CELL_output.TYPE;
if (isfield(TOOLBOX_input.Scene.module_mounting, 'avgSensitivity')==1)
    avg_sen = TOOLBOX_input.Scene.module_mounting.avgSensitivity;
else
    avg_sen = 1;
end

%% Calculate the Reflectance, Parasitic absorption, and Absorption
RAT_front = CELL_output.CELL_FRONT.RAT;
if strcmp(type,'SHJ')|| strcmp(type,'T-F') %For a monofacial cell
    R = RAT_front(:,:,1)'+RAT_front(:,:,end)';
    A_diff = 1 - sum(RAT_front(:,:,:),3)';
    A = RAT_front(:,:,2)';
    Bifacial = 0;
    Angle_emit = 2*pi;
elseif strcmp(type,'BIF')  %For a bifiacial cell
    RAT_rear = CELL_output.CELL_REAR.RAT;
    R = cat(3,RAT_front(:,:,1)'+RAT_front(:,:,end)',RAT_rear(:,:,1)'+RAT_rear(:,:,end)');
    A_diff = cat(3,1 - sum(RAT_front(:,:,:),3)',1 - sum(RAT_rear(:,:,:),3)');
    A = cat(3,RAT_front(:,:,2)',RAT_rear(:,:,2)');
    Bifacial = 1;
    Angle_emit = 4*pi;
elseif strcmp(type,'Tan') % For a tandem cell
    R = RAT_front(:,:,1)'+RAT_front(:,:,end)';
    A_diff = 1 - sum(RAT_front(:,:,:),3)';
    A1 = RAT_front(:,:,2)';
    A = A1;
    A2 = RAT_front(:,:,3)';
    Bifacial = 0;
    Angle_emit = 2*pi;
elseif strcmp(type,'BIF-Tan')
    RAT_rear = CELL_output.CELL_REAR.RAT;
    R = cat(3,RAT_front(:,:,1)'+RAT_front(:,:,end)',RAT_rear(:,:,1)'+RAT_rear(:,:,end)');
    A_diff = cat(3,1 - sum(RAT_front(:,:,:),3)',1 - sum(RAT_rear(:,:,:),3)');
    A1 = cat(3,RAT_front(:,:,2)',RAT_rear(:,:,2)');
    A2 = cat(3,RAT_front(:,:,3)',RAT_rear(:,:,3)');
    Bifacial = 1;
    Angle_emit = 4*pi;
end

WEATHER_J = WEATHER_output.J;
THERMAL_T = THERMAL_output.T;
Parameters_1 = ELECTRIC_output.Parameters_1;
Parameters_2 = ELECTRIC_output.Parameters_2;
ModTilt = MODULE_output.ModTilt;
Ncells = MODULE_output.N;
Amod = MODULE_output.Amod;
AZA = MODULE_output.skydome.AZA;
SM_f = MODULE_output.SM_f;
if Bifacial; SM_r = MODULE_output.SM_r; end
V_mpp = ELECTRIC_output.Vmpp;
I_mpp = ELECTRIC_output.Impp;
if (isfield(TOOLBOX_input, 'runACConversionPart')==1 && TOOLBOX_input.runACConversionPart == 1)
    CONVERSION_Pac = CONVERSION_output.Pac;
else
    CONVERSION_Pac = [];
end
angles_GENPRO = CELL_output.CELL_FRONT.aoi;
if (isfield(TOOLBOX_input.electric, 'Terminals')==1)
    Terminals = TOOLBOX_input.electric.Terminals;
else
    Terminals = 2;
end

if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F') %If it is a single junction
    days = find(sum(WEATHER_J,2));
    E_g = TOOLBOX_input.LossAnalysis.E_g;
    Parameters = zeros(Ncells,length(DCP),5);
    Parameters(:,days,:) = reshape(Parameters_1,Ncells,length(days),5);
elseif strcmp(type,'Tan')  || strcmp(type,'BIF-Tan')%For a tandem cell
    days = find(sum(WEATHER_J(:,:,1),2));
    E_g1 = TOOLBOX_input.LossAnalysis.E_g1;
    E_g2 = TOOLBOX_input.LossAnalysis.E_g2;
    Parameters1 = zeros(Ncells,length(DCP),5);
    Parameters1(:,days,:) = reshape(Parameters_1,Ncells,length(days),5);
    Parameters2 = zeros(Ncells,length(DCP),5);
    Parameters2(:,days,:) = reshape(Parameters_2,Ncells,length(days),5);
end


%% Prepare weather data
weather_data = load_meteonorm_data(TOOLBOX_input);
wav_orig = MODULE_output.wav*1e3;
if TOOLBOX_input.irradiation.spectra_choice == 1
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,air_mass,wav] = spectral_distrSMARTS(wav_orig*1e-3);
    wav = wav*1e-6;
elseif TOOLBOX_input.irradiation.spectra_choice == 2
    [RSD_i_dir,RSD_f_dir,RSD_i_dif,RSD_f_dif,air_mass,wav] = spectral_distrSBDarts(wav_orig*1e-3);
    wav = wav*1e-6;
end

N_rep_before = find(wav*1e9 == wav_orig(1))-1;
N_rep_after = length(wav) - find(wav*1e9 == wav_orig(end));
SM_f = cat(4,repelem(SM_f(:,:,:,1),1,1,1,N_rep_before),SM_f,repelem(SM_f(:,:,:,end),1,1,1,N_rep_after));
if Bifacial; SM_r = cat(4,repelem(SM_r(:,:,:,1),1,1,1,N_rep_before),SM_r,repelem(SM_r(:,:,:,end),1,1,1,N_rep_after)); end

SpecData.AM = air_mass;
SpecData.WAV = wav;
SpecData.RSD_f_dir = RSD_f_dir;
SpecData.RSD_f_dif = RSD_f_dif;
SpecData.RSD_i_dir = RSD_i_dir;
SpecData.RSD_i_dif = RSD_i_dif;
SpecData.spectra_choice = TOOLBOX_input.irradiation.spectra_choice;

%%  Module parameters
Mod_tilt = ModTilt;
Mod_alti = 90-Mod_tilt;
Mod_azi = TOOLBOX_input.Scene.module_mounting.ModAzimuth;
T_cell = mean(THERMAL_T,2)+273.15;
Angle_abs = 67.7e-6; %https://en.wikipedia.org/wiki/Solid_angle



%% Incoming irradiance
% ---- Calculate the extraterrestrial solar radiation
% Compute the day of the year for the Perez model
day_year = pvl_date2doy(weather_data(:,1), weather_data(:,2), ...
    weather_data(:,3));
% Expression used for the calculation is eq. 1.4.1b in "Solar engineering
% of thermal processes" by J. Duffie et al.
day_angle = 2*pi*(day_year - 1)/365.25;
E0 = 1.00011 + 0.034221*cos(day_angle) + 0.00128*sin(day_angle) + ...
    0.000719*cos(2*day_angle) + 0.000077*sin(2*day_angle);
[src_folder,~,~] = get_folder_structure;
load(fullfile(src_folder,'3_WEATHER','constants','weather_params.mat'), 'AM0')
extra_sol_power = E0*AM0;
if ~Bifacial
    Irradiance_Year = zeros(length(DCP),15,length(wav));
    Photon_flux_Year = zeros(length(DCP),15,length(wav));
    for hour_index = 1:length(DCP)
        if DCP(hour_index) == 0
            continue
        end
        Incoming_Irr_input = [];
        Incoming_Irr_input.wav = wav;
        Incoming_Irr_input.Mod_alti = Mod_alti;
        Incoming_Irr_input.Mod_azi = Mod_azi;
        Incoming_Irr_input.Bifacial = Bifacial;
        Incoming_Irr_input.J_abs = mean(WEATHER_J(hour_index,:,1));
        Incoming_Irr_input.AZA = AZA;
        Incoming_Irr_input.index = hour_index;
        Incoming_Irr_input.A = A;
        Incoming_Irr_input.extra_sol_power = extra_sol_power;
        Incoming_Irr_input.SM_f = SM_f;
        Incoming_Irr_input.Ncells = Ncells;
        % Calculate the incoming irradiance
        [Flux_angles,Irradiance_angles] = incoming_irradiance(Incoming_Irr_input,TOOLBOX_input,CELL_output,SpecData,weather_data,CONSTANTS);
        Irradiance_Year(hour_index,:,:) = Irradiance_angles;
        Photon_flux_Year(hour_index,:,:) = Flux_angles;
    end
else
    Irradiance_Year = zeros(length(DCP),15,length(wav),2);
    Photon_flux_Year = zeros(length(DCP),15,length(wav),2);
    if strcmp(type,'BIF')
        for hour_index = 1:length(DCP)
            if DCP(hour_index) == 0
                continue
            end
            Incoming_Irr_input = [];
            Incoming_Irr_input.wav = wav;
            Incoming_Irr_input.Mod_alti = Mod_alti;
            Incoming_Irr_input.Mod_azi = Mod_azi;
            Incoming_Irr_input.Bifacial = Bifacial;
            Incoming_Irr_input.J_abs = mean(WEATHER_J(hour_index,:,1));
            Incoming_Irr_input.AZA = AZA;
            Incoming_Irr_input.index = hour_index;
            Incoming_Irr_input.A = A;
            Incoming_Irr_input.extra_sol_power = extra_sol_power;
            Incoming_Irr_input.SM_f = SM_f;
            Incoming_Irr_input.SM_r = SM_r;
            % Calculate the incoming irradiance
            [Flux_angles,Irradiance_angles] = incoming_irradiance(Incoming_Irr_input,TOOLBOX_input,CELL_output,SpecData,weather_data,CONSTANTS);
            Irradiance_Year(hour_index,:,:,:) = Irradiance_angles;
            Photon_flux_Year(hour_index,:,:,:) = Flux_angles;
        end
    elseif strcmp(type,'BIF-Tan')
        for hour_index = 1:length(DCP)
            if DCP(hour_index) == 0
                continue
            end
            Incoming_Irr_input = [];
            Incoming_Irr_input.wav = wav;
            Incoming_Irr_input.Mod_alti = Mod_alti;
            Incoming_Irr_input.Mod_azi = Mod_azi;
            Incoming_Irr_input.J_abs = mean(WEATHER_J(hour_index,:,:));
            Incoming_Irr_input.AZA = AZA;
            Incoming_Irr_input.index = hour_index;
            Incoming_Irr_input.A1 = A1;
            Incoming_Irr_input.A2 = A2;
            Incoming_Irr_input.extra_sol_power = extra_sol_power;
            Incoming_Irr_input.SM_f = SM_f;
            Incoming_Irr_input.SM_r = SM_r;
            % Calculate the incoming irradiance
            [Flux_angles,Irradiance_angles] = incoming_irradiance_BIF_TAN(Incoming_Irr_input,TOOLBOX_input,CELL_output,SpecData,weather_data,CONSTANTS);
            Irradiance_Year(hour_index,:,:,:) = Irradiance_angles;
            Photon_flux_Year(hour_index,:,:,:) = Flux_angles;
        end
    end
end
disp('The irradiance has been calculated')


if strcmp(type,'SHJ') || strcmp(type,'BIF')
    I_abs = zeros(length(DCP),length(angles_GENPRO));
    % The optimal voltage is calculated
    V_opt1 = Vopt_calculator(wav,squeeze(sum(sum(Photon_flux_Year,4),2)),E_g,T_cell,Angle_emit);
elseif strcmp(type,'Tan')
    I_abs1 = zeros(length(DCP),length(angles_GENPRO));
    I_abs2 = zeros(length(DCP),length(angles_GENPRO));
    % The optimal voltages are calculated
    V_opt1 = Vopt_calculator(wav,squeeze(sum(sum(Photon_flux_Year,4),2)),E_g1,T_cell,Angle_emit);
    V_opt2 = Vopt_calculator(wav,squeeze(sum(sum(Photon_flux_Year,4),2)),E_g2,T_cell,Angle_emit);
elseif strcmp(type,'BIF-Tan')
    I_abs1 = zeros(length(DCP),length(angles_GENPRO));
    I_abs2 = zeros(length(DCP),length(angles_GENPRO));
    % The optimal voltages are calculated
    V_opt1 = Vopt_calculator(wav,squeeze(sum(sum(Photon_flux_Year,4),2)),E_g1,T_cell,Angle_emit);
    V_opt2 = Vopt_calculator(wav,squeeze(sum(sum(Photon_flux_Year,4),2)),E_g2,T_cell,Angle_emit);

    %Find the voltage over cell 2
    Vth = CONSTANTS.k*T_cell/CONSTANTS.q;
    Iph1 = mean(Parameters1(:,:,1))';
    Rs1 = mean(Parameters1(:,:,2))';
    Rsh1 = mean(Parameters1(:,:,3))';
    n1 = mean(Parameters1(:,:,4))';
    I01 = mean(Parameters1(:,:,5))';
    V_cell1 = Iph1.*Rsh1+I01.*Rsh1-I_mpp.*(Rs1+Rsh1)-n1.*Vth.*lambertw(Rsh1.*I01./(n1.*Vth).*exp(Rsh1./(n1.*Vth).*(Iph1+I01-I_mpp)));
    V_cell2 = V_mpp-MODULE_output.N*V_cell1;
    V_cell2(isnan(V_cell2)) = 0;
end

% The fundamental and optical losses are calculated per angle. First,
% all values are initialised
P_in = zeros(length(DCP),length(angles_GENPRO));
P_term = zeros(length(DCP),length(angles_GENPRO));
P_below = zeros(length(DCP),length(angles_GENPRO));
P_carnot = zeros(length(DCP),length(angles_GENPRO));
P_emission = zeros(length(DCP),length(angles_GENPRO));
P_angle = zeros(length(DCP),length(angles_GENPRO));
P_gain = zeros(length(DCP),length(angles_GENPRO));
P_fund = zeros(length(DCP),length(angles_GENPRO));
P_cell = zeros(length(DCP),length(angles_GENPRO));
P_metal = zeros(length(DCP),length(angles_GENPRO));
P_ref = zeros(length(DCP),length(angles_GENPRO));
P_diffA = zeros(length(DCP),length(angles_GENPRO));


%For each angle, the optical and fundamental losses are calculated
for i = 1:length(angles_GENPRO)
    AOI = angles_GENPRO(i);
    %input power
    if Bifacial == 0
        P_in(:,i)=trapz(wav,squeeze(Irradiance_Year(:,i,:))')*Amod;
    elseif Bifacial ==1
        P_in_fr = trapz(wav,squeeze(Irradiance_Year(:,i,:,1))')*Amod;
        P_in_rr = trapz(wav,squeeze(Irradiance_Year(:,i,:,2))')*Amod;
        P_in(:,i)= P_in_fr+P_in_rr;
    end
    % The optical parameters for a certain angle are calculated
    if strcmp(type,'SHJ') || strcmp(type,'BIF')
        A_angle = interp1(angles_GENPRO,A,AOI);
        R_angle = interp1(angles_GENPRO,R,AOI);
        A_diff_angle = interp1(angles_GENPRO,A_diff,AOI);
    elseif strcmp(type,'Tan')|| strcmp(type,'BIF-Tan')
        A1_angle = interp1(angles_GENPRO,A1,AOI);
        A2_angle = interp1(angles_GENPRO,A2,AOI);
        R_angle = interp1(angles_GENPRO,R,AOI);
        A_diff_angle = interp1(angles_GENPRO,A_diff,AOI);
    end
    %Fundamental losses
    if strcmp(type,'SHJ')
        Fund_Losses_input.wav = wav;
        Fund_Losses_input.Irr_spec = squeeze(Irradiance_Year(:,i,:))';
        Fund_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:))';
        Fund_Losses_input.Angle_abs = Angle_abs;
        Fund_Losses_input.Angle_emit_Vopt = Angle_emit;
        Fund_Losses_input.Angle_emit_emission = Angle_emit/length(angles_GENPRO);
        Fund_Losses_input.A = A_angle;
        Fund_Losses_input.E_g = E_g;
        Fund_Losses_input.I_mpp = I_mpp';
        Fund_Losses_input.V_mpp = V_mpp';
        Fund_Losses_input.T_cell = T_cell';
        Fund_Losses_input.V_opt1 = V_opt1;
        [P_term(:,i),P_below(:,i),P_emission(:,i),P_carnot(:,i),P_angle(:,i),P_gain(:,i)] = FundamentalLossesSingle(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);
    elseif strcmp(type,'BIF')
        Fund_Losses_input.wav = wav;
        Fund_Losses_input.Irr_spec = squeeze(Irradiance_Year(:,i,:,1))';
        Fund_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,1))';
        Fund_Losses_input.Angle_abs = Angle_abs;
        Fund_Losses_input.Angle_emit_Vopt = Angle_emit;
        Fund_Losses_input.Angle_emit_emission = 2*pi/length(angles_GENPRO);
        Fund_Losses_input.A = A_angle(:,:,1);
        Fund_Losses_input.E_g = E_g;
        Fund_Losses_input.I_mpp = I_mpp';
        Fund_Losses_input.V_mpp = V_mpp';
        Fund_Losses_input.T_cell = T_cell';
        Fund_Losses_input.V_opt1 = V_opt1;
        [P_term_fr,P_below_fr,P_emission_fr,P_carnot_fr,P_angle_fr,P_gain_fr] = FundamentalLossesSingle(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);

        Fund_Losses_input.Irr_spec = squeeze(Irradiance_Year(:,i,:,2))';
        Fund_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,2))';
        Fund_Losses_input.A = A_angle(:,:,2);
        [P_term_rr,P_below_rr,P_emission_rr,P_carnot_rr,P_angle_rr,P_gain_rr] = FundamentalLossesSingle(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);

        P_term(:,i) = P_term_fr+P_term_rr;
        P_below(:,i) = P_below_fr+P_below_rr;
        P_emission(:,i) = P_emission_fr+P_emission_rr;
        P_carnot(:,i) = P_carnot_fr + P_carnot_rr;
        P_angle(:,i) = P_angle_fr+P_angle_rr;
        P_gain(:,i) = P_gain_fr + P_gain_rr;
        P_fund_fr = P_term_fr+P_below_fr+P_emission_fr+P_carnot_fr+P_angle_fr-P_gain_fr;
        P_fund_rr = P_term_rr+P_below_rr+P_emission_rr+P_carnot_rr+P_angle_rr-P_gain_rr;
    elseif strcmp(type,'Tan')
        Fund_Losses_input.wav = wav;
        Fund_Losses_input.Irr_spec = squeeze(Irradiance_Year(:,i,:))';
        Fund_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:))';
        Fund_Losses_input.Angle_abs = Angle_abs;
        Fund_Losses_input.Angle_emit_Vopt = Angle_emit;
        Fund_Losses_input.Angle_emit_emission = Angle_emit/length(angles_GENPRO);
        Fund_Losses_input.A1 = A1_angle;
        Fund_Losses_input.A2 = A2_angle;
        Fund_Losses_input.E_g1 = E_g1;
        Fund_Losses_input.E_g2 = E_g2;
        Fund_Losses_input.Parameters1 = squeeze(mean(Parameters1));
        Fund_Losses_input.I_mpp = I_mpp';
        Fund_Losses_input.V_mpp = V_mpp';
        Fund_Losses_input.T_cell = T_cell';
        Fund_Losses_input.V_opt1 = V_opt1;
        Fund_Losses_input.V_opt2 = V_opt2;
        [P_term(:,i),P_below(:,i),P_emission(:,i),P_carnot(:,i),P_angle(:,i),P_gain(:,i)] = FundamentalLossesTandem(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);
    elseif strcmp(type,'BIF-Tan')
        Fund_Losses_input.wav = wav;
        Fund_Losses_input.Irr_spec = squeeze(Irradiance_Year(:,i,:,1))';
        Fund_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,1))';
        Fund_Losses_input.Angle_abs = Angle_abs;
        Fund_Losses_input.Angle_emit_Vopt = Angle_emit;
        Fund_Losses_input.Angle_emit_emission = 2*pi/length(angles_GENPRO);
        Fund_Losses_input.A1 = A1_angle(:,:,1);
        Fund_Losses_input.A2 = A2_angle(:,:,1);
        Fund_Losses_input.E_g1 = E_g1;
        Fund_Losses_input.E_g2 = E_g2;
        Fund_Losses_input.Parameters1 = squeeze(mean(Parameters1));
        Fund_Losses_input.I_mpp = I_mpp';
        Fund_Losses_input.V_mpp = V_mpp';
        Fund_Losses_input.T_cell = T_cell';
        Fund_Losses_input.V_opt1 = V_opt1;
        Fund_Losses_input.V_opt2 = V_opt2;
        [P_term_fr,P_below_fr,P_emission_fr,P_carnot_fr,P_angle_fr,P_gain_fr] = FundamentalLossesTandem(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);

        Fund_Losses_input.Irr_spec = squeeze(Irradiance_Year(:,i,:,2))';
        Fund_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,2))';
        Fund_Losses_input.A = A2_angle(:,:,2);
        Fund_Losses_input.V_mpp = V_cell2';
        Fund_Losses_input.V_opt1 = V_opt2;
        Fund_Losses_input.E_g = E_g2;
        [P_term_rr,P_below_rr,P_emission_rr,P_carnot_rr,P_angle_rr,P_gain_rr] = FundamentalLossesSingle(Fund_Losses_input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS);

        P_term(:,i) = P_term_fr+P_term_rr;
        P_below(:,i) = P_below_fr+P_below_rr;
        P_emission(:,i) = P_emission_fr+P_emission_rr;
        P_carnot(:,i) = P_carnot_fr + P_carnot_rr;
        P_angle(:,i) = P_angle_fr+P_angle_rr;
        P_gain(:,i) = P_gain_fr + P_gain_rr;
        P_fund_fr = P_term_fr+P_below_fr+P_emission_fr+P_carnot_fr+P_angle_fr-P_gain_fr;
        P_fund_rr = P_term_rr+P_below_rr+P_emission_rr+P_carnot_rr+P_angle_rr-P_gain_rr;
    end
    P_fund(:,i) = P_term(:,i)+P_below(:,i)+P_emission(:,i)+P_carnot(:,i)+P_angle(:,i)-P_gain(:,i);
    %Optical losses
    if strcmp(type,'SHJ')
        Opt_Losses_input.P_in = P_in(:,i);
        Opt_Losses_input.P_fund = P_fund(:,i);
        Opt_Losses_input.wav = wav;
        Opt_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:))';
        Opt_Losses_input.R = R_angle;
        Opt_Losses_input.A = A_angle;
        Opt_Losses_input.A_diff = A_diff_angle;
        Opt_Losses_input.E_g = E_g;
        Opt_Losses_input.V_opt1 = V_opt1;
        Opt_Losses_input.Angle_emit = Angle_emit/length(angles_GENPRO);
        Opt_Losses_input.V_mpp = V_mpp';
        Opt_Losses_input.I_mpp = I_mpp';
        Opt_Losses_input.T_cell = T_cell';
        [P_cell(:,i),P_metal(:,i),P_ref(:,i),P_diffA(:,i),I_abs(:,i)] = OpticalLossesSingle(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);

    elseif strcmp(type,'BIF')
        Opt_Losses_input.P_in = P_in_fr;
        Opt_Losses_input.P_fund = P_fund_fr;
        Opt_Losses_input.wav = wav;
        Opt_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,1))';
        Opt_Losses_input.R = R_angle(:,:,1);
        Opt_Losses_input.A = A_angle(:,:,1);
        Opt_Losses_input.A_diff = A_diff_angle(:,:,1);
        Opt_Losses_input.E_g = E_g;
        Opt_Losses_input.V_opt1 = V_opt1;
        Opt_Losses_input.Angle_emit = 2*pi/length(angles_GENPRO);
        Opt_Losses_input.V_mpp = V_mpp';
        Opt_Losses_input.I_mpp = I_mpp';
        Opt_Losses_input.T_cell = T_cell';
        [P_cell_fr,P_metal_fr,P_ref_fr,P_diffA_fr,I_abs_fr] = OpticalLossesSingle(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);

        Opt_Losses_input.P_in = P_in_rr;
        Opt_Losses_input.P_fund = P_fund_rr;
        Opt_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,2))';
        Opt_Losses_input.R = R_angle(:,:,2);
        Opt_Losses_input.A = A_angle(:,:,2);
        Opt_Losses_input.A_diff = A_diff_angle(:,:,2);
        [P_cell_rr,P_metal_rr,P_ref_rr,P_diffA_rr,I_abs_rr] = OpticalLossesSingle(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);

        P_cell(:,i) = P_cell_fr + P_cell_rr;
        P_metal(:,i) = P_metal_fr + P_metal_rr;
        P_ref(:,i) = P_ref_fr + P_ref_rr;
        P_diffA(:,i) = P_diffA_fr + P_diffA_rr;
        I_abs(:,i) = I_abs_fr + I_abs_rr;
    elseif strcmp(type,'Tan')
        Opt_Losses_input.P_in = P_in(:,i);
        Opt_Losses_input.P_fund = P_fund(:,i);
        Opt_Losses_input.wav = wav;
        Opt_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:))';
        Opt_Losses_input.R = R_angle;
        Opt_Losses_input.A1 = A1_angle;
        Opt_Losses_input.A2 = A2_angle;
        Opt_Losses_input.A_diff = A_diff_angle;
        Opt_Losses_input.E_g1 = E_g1;
        Opt_Losses_input.E_g2 = E_g2;
        Opt_Losses_input.V_opt1 = V_opt1;
        Opt_Losses_input.V_opt2 = V_opt2;
        Opt_Losses_input.Angle_emit = Angle_emit/length(angles_GENPRO);
        Opt_Losses_input.Parameters1 = squeeze(mean(Parameters1));
        Opt_Losses_input.V_mpp = V_mpp';
        Opt_Losses_input.I_mpp = I_mpp';
        Opt_Losses_input.T_cell = T_cell';
        [P_cell(:,i),P_metal(:,i),P_ref(:,i),P_diffA(:,i),I_abs1(:,i),I_abs2(:,i)] = OpticalLossesTandem(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);
    elseif strcmp(type,'BIF-Tan')
        Opt_Losses_input.P_in = P_in_fr;
        Opt_Losses_input.P_fund = P_fund_fr;
        Opt_Losses_input.wav = wav;
        Opt_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,1))';
        Opt_Losses_input.R = R_angle(:,:,1);
        Opt_Losses_input.A1 = A1_angle(:,:,1);
        Opt_Losses_input.A2 = A2_angle(:,:,1);
        Opt_Losses_input.A_diff = A_diff_angle(:,:,1);
        Opt_Losses_input.E_g1 = E_g1;
        Opt_Losses_input.E_g2 = E_g2;
        Opt_Losses_input.V_opt1 = V_opt1;
        Opt_Losses_input.V_opt2 = V_opt2;
        Opt_Losses_input.Angle_emit = 2*pi/length(angles_GENPRO);
        Opt_Losses_input.Parameters1 = squeeze(mean(Parameters1));
        Opt_Losses_input.V_mpp = V_mpp';
        Opt_Losses_input.I_mpp = I_mpp';
        Opt_Losses_input.T_cell = T_cell';
        [P_cell_fr,P_metal_fr,P_ref_fr,P_diffA_fr,I_abs1_fr,I_abs2_fr] = OpticalLossesTandem(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);

        Opt_Losses_input.P_in = P_in_rr;
        Opt_Losses_input.P_fund = P_fund_rr;
        Opt_Losses_input.photon_spec = squeeze(Photon_flux_Year(:,i,:,2))';
        Opt_Losses_input.E_g = E_g2;
        Opt_Losses_input.V_opt1 = V_opt2;
        Opt_Losses_input.V_mpp = V_cell2';
        Opt_Losses_input.R = R_angle(:,:,2);
        Opt_Losses_input.A = A2_angle(:,:,2);
        Opt_Losses_input.A_diff = A_diff_angle(:,:,2);
        [P_cell_rr,P_metal_rr,P_ref_rr,P_diffA_rr,I_abs2_rr] = OpticalLossesSingle(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS);

        P_cell(:,i) = P_cell_fr + P_cell_rr;
        P_metal(:,i) = P_metal_fr + P_metal_rr;
        P_ref(:,i) = P_ref_fr + P_ref_rr;
        P_diffA(:,i) = P_diffA_fr + P_diffA_rr;
        I_abs1(:,i) = I_abs1_fr;
        I_abs2(:,i) = I_abs2_fr+ I_abs2_rr;
    end
end

P_in = sum(P_in,2);
P_term = sum(P_term,2);
P_below = sum(P_below,2);
P_carnot = sum(P_carnot,2);
P_emission = sum(P_emission,2);
P_angle = sum(P_angle,2);
P_gain = sum(P_gain,2);
P_cell = sum(P_cell,2);
P_metal = sum(P_metal,2);
P_ref = sum(P_ref,2);
P_diffA = sum(P_diffA,2);
if strcmp(type,'SHJ') || strcmp(type,'BIF')
    I_abs = sum(I_abs,2);
    IV_curvename = TOOLBOX_input.electric.IVtype;
    [P_diffA, I_abs] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs,T_cell,V_opt1',IV_curvename);
elseif strcmp(type,'Tan')|| strcmp(type,'BIF-Tan')    
    I_abs1 = sum(I_abs1,2);
    I_abs2 = sum(I_abs2,2);
    IV_curvename_top = TOOLBOX_input.electric.IVtypeTop;
    IV_curvename_bot = TOOLBOX_input.electric.IVtypeBot;
    [P_diffA, I_abs1] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs1,T_cell,V_opt1',IV_curvename_top);
    [P_diffA, I_abs2] = Correct_PdiffA(TOOLBOX_input,MODULE_output,P_diffA,I_abs2,T_cell,V_opt2',IV_curvename_bot);
end


%% Electrical losses

if strcmp(type,'SHJ') || strcmp(type,'BIF')
    Elec_Losses_input.I_abs = I_abs;
    Elec_Losses_input.V_opt1 = V_opt1;
    Elec_Losses_input.Parameters = Parameters;
    Elec_Losses_input.T_cell = T_cell;
    [P_series,P_shunt,P_NRRI,P_NRRV] = ElectricLossesSingle_avg(Elec_Losses_input,MODULE_output);
    PowerRatio = zeros(size(DCP));
elseif strcmp(type,'Tan')|| strcmp(type,'BIF-Tan')
    Elec_Losses_input.I_abs1 = I_abs1;
    Elec_Losses_input.I_abs2 = I_abs2;
    Elec_Losses_input.V_opt1 = V_opt1;
    Elec_Losses_input.V_opt2 = V_opt2;
    Elec_Losses_input.Parameters1 = Parameters1;
    Elec_Losses_input.Parameters2 = Parameters2;
    Elec_Losses_input.T_cell = T_cell;
    [P_series,P_shunt,P_NRRI,P_NRRV,PowerRatio] = ElectricLossesTandem_avg(Elec_Losses_input,TOOLBOX_input,MODULE_output,CONSTANTS);

end

%% System losses

if strcmp(type,'SHJ') || strcmp(type,'BIF')
    Sys_Losses_input.I_mpp = I_mpp;
    Sys_Losses_input.V_mpp = V_mpp;
    Sys_Losses_input.I_abs = I_abs;
    Sys_Losses_input.Parameters = Parameters;
    Sys_Losses_input.T_cell = T_cell;
    Sys_Losses_input.Pdc = DCP;
    if (isfield(TOOLBOX_input, 'runACConversionPart')==1) %To check whether the AC simulation is performed
        if TOOLBOX_input.runACConversionPart == 1
            Sys_Losses_input.Pac = CONVERSION_Pac;
        end
    end
    [P_con,P_mismatch,P_cable,P_inv] = SystemLossesSingle(Sys_Losses_input,TOOLBOX_input,ELECTRIC_output,MODULE_output);
elseif strcmp(type,'Tan')|| strcmp(type,'BIF-Tan')
    if Terminals == 2 || Terminals ==3
        Sys_Losses_input.I_mpp = I_mpp;
        Sys_Losses_input.V_mpp = V_mpp;
        Sys_Losses_input.I_abs1 = I_abs1;
        Sys_Losses_input.I_abs2 = I_abs2;
        Sys_Losses_input.Parameters1 = Parameters1;
        Sys_Losses_input.Parameters2 = Parameters2;
        Sys_Losses_input.T_cell = T_cell;
        Sys_Losses_input.Pdc = DCP;
        if TOOLBOX_input.runACConversionPart == 1
            Sys_Losses_input.Pac = CONVERSION_Pac;
        end
        [P_con,P_mismatch,P_cable,P_inv] = SystemLossesTandem2T(Sys_Losses_input,TOOLBOX_input,ELECTRIC_output,MODULE_output);
    end
end
Index_zero = find(~DCP);
P_term(Index_zero) = 0;
P_below(Index_zero) = 0;
P_carnot(Index_zero) = 0;
P_emission(Index_zero) = 0;
P_angle(Index_zero) = 0;
P_gain(Index_zero) = 0;
P_cell(Index_zero) = 0;
P_metal(Index_zero) = 0;
P_ref(Index_zero) = 0;
P_diffA(Index_zero) = 0;
P_NRRI(Index_zero) = 0;
P_shunt(Index_zero) = 0;
P_series(Index_zero) = 0;
P_con(Index_zero) = 0;
P_NRRV(Index_zero) = 0;
P_mismatch(Index_zero) = 0;
P_cable(Index_zero) = 0;
P_inv(Index_zero) = 0;
P_in(Index_zero) = 0;
PowerRatio(Index_zero) = 0;
P_out_DC = DCP.*[PowerRatio, (1-PowerRatio)];
if (isfield(TOOLBOX_input, 'runACConversionPart')==1)&& TOOLBOX_input.runACConversionPart == 1 %To check whether the AC simulation is performed
    P_out_AC = max(CONVERSION_Pac.*[PowerRatio, (1-PowerRatio)],0);
end


%% Total power
disp('DC losses')
Power = [sum(P_term);sum(P_below);sum(P_emission);sum(P_carnot);sum(P_angle);-1*sum(P_gain);sum(P_cell);sum(P_metal);sum(P_ref);sum(P_diffA);sum(P_NRRI);sum(P_shunt);sum(P_series);sum(P_con);sum(P_NRRV);sum(P_mismatch);sum(sum(P_out_DC))];
P_total = sum(Power);
Power = [Power; P_total];
Components = ["Thermalization";"Below bandgap";"Emission";"Carnot losses";"Angle mismatch";"Gain";"Cell spacing";"Metal shading";"Reflection";"Parasitic absorption";"Recombination I";"Shunt resistance";"Series Resistance";"Connection loss";"Recombination V";"Mismatch losses";"DC power"; "Total"];
Percentage = 100*Power/sum(P_in);
disp(table(Components,Power,Percentage))


if (isfield(TOOLBOX_input, 'runACConversionPart')==1 && TOOLBOX_input.runACConversionPart == 1) %To check whether the AC simulation is performed
    disp('AC losses')
    panels = TOOLBOX_input.Conversion.Parallel_Modules*TOOLBOX_input.Conversion.Series_Modules;
    Power = panels*[sum(P_term);sum(P_below);sum(P_emission);sum(P_carnot);sum(P_angle);-1*sum(P_gain);sum(P_cell);sum(P_metal);sum(P_ref);sum(P_diffA);sum(P_NRRI);sum(P_shunt);sum(P_series);sum(P_con);sum(P_NRRV);sum(P_mismatch)];
    Power = [Power;sum(P_cable);sum(P_inv);sum(sum(P_out_AC))];
    P_total = sum(Power);
    Power = [Power; P_total];
    Components = ["Thermalization";"Below bandgap";"Emission";"Carnot losses";"Angle mismatch";"Gain";"Cell spacing";"Metal shading";"Reflection";"Parasitic absorption";"Recombination I";"Shunt resistance";"Series Resistance";"Connection loss";"Recombination V";"Mismatch losses";"Cable losses";"Inverter losses" ;"AC Power"; "Total"];
    Percentage = 100*Power/(sum(P_in)*panels);
    disp(table(Components,Power,Percentage))

end

Losses_Operating.P_term_full = P_term;
Losses_Operating.P_below_full = P_below;
Losses_Operating.P_carnot_full = P_carnot;
Losses_Operating.P_emission_full = P_emission;
Losses_Operating.P_angle_full = P_angle;
Losses_Operating.P_gain_full = P_gain;
Losses_Operating.P_cell_full = P_cell;
Losses_Operating.P_metal_full = P_metal;
Losses_Operating.P_ref_full = P_ref;
Losses_Operating.P_diffA_full = P_diffA;
Losses_Operating.P_NRRI_full = P_NRRI;
Losses_Operating.P_shunt_full = P_shunt;
Losses_Operating.P_series_full = P_series;
Losses_Operating.P_con_full = P_con;
Losses_Operating.P_NRRV_full = P_NRRV;
Losses_Operating.P_mismatch_full = P_mismatch;
Losses_Operating.P_cable_full = P_cable;
Losses_Operating.P_inv_full = P_inv;
Losses_Operating.P_in_full = P_in;
Losses_Operating.P_out_DC_full = P_out_DC;
if (isfield(TOOLBOX_input, 'runACConversionPart')==1 && TOOLBOX_input.runACConversionPart == 1)
    Losses_Operating.P_out_AC_full = P_out_AC;
end

%% Plot results
if TOOLBOX_input.LossAnalysis.plotFigures == 1
    if (isfield(TOOLBOX_input, 'runACConversionPart')==1 && TOOLBOX_input.runACConversionPart == 1)
        if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F') %AC simulation for single junction
            Plot_Components = ["AC power";"Inverter losses";"Cable losses";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Power = 100*[sum(P_out_AC(:,2));sum(P_inv);sum(P_cable);sum(P_mismatch)*panels;sum(P_con)*panels;sum(P_NRRI)*panels;sum(P_NRRV)*panels;sum(P_shunt)*panels;sum(P_series)*panels;sum(P_diffA)*panels;sum(P_ref)*panels;sum(P_metal)*panels;sum(P_cell)*panels;-sum(P_gain)*panels;sum(P_angle)*panels;sum(P_carnot)*panels;sum(P_emission)*panels;sum(P_below)*panels;sum(P_term)*panels]/(sum(P_in)*panels);
            Plot_Categories = ["Fundamental losses:", "Optical losses:", "Electrical losses:","System losses:","Output Power"];
            Plot_Power_Categories = 100*[sum(P_term)*panels+sum(P_below)*panels+sum(P_angle)*panels+sum(P_carnot)*panels+sum(P_emission)*panels-sum(P_gain)*panels;sum(P_cell)*panels+sum(P_metal)*panels+sum(P_ref)*panels+sum(P_diffA)*panels;sum(P_NRRV)*panels+sum(P_NRRI)*panels+sum(P_series)*panels+sum(P_shunt)*panels;sum(P_mismatch)*panels+sum(P_con)*panels+sum(P_cable)+sum(P_inv);sum(P_out_AC(:,2))]/(sum(P_in)*panels);
            Plot_type = 1;
        elseif strcmp(type,'Tan')  || strcmp(type,'BIF-Tan')%DC simulation for tandem cell
            Plot_Components = ["AC power bottom cell";"AC power top cell";"Inverter losses";"Cable losses";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Power = 100*[sum(P_out_AC(:,2));sum(P_out_AC(:,1));sum(P_inv);sum(P_cable);sum(P_mismatch)*panels;sum(P_con)*panels;sum(P_NRRI)*panels;sum(P_NRRV)*panels;sum(P_shunt)*panels;sum(P_series)*panels;sum(P_diffA)*panels;sum(P_ref)*panels;sum(P_metal)*panels;sum(P_cell)*panels;-sum(P_gain)*panels;sum(P_angle)*panels;sum(P_carnot)*panels;sum(P_emission)*panels;sum(P_below)*panels;sum(P_term)*panels]/(sum(P_in)*panels);
            Plot_Categories = ["Fundamental losses:", "Optical losses:", "Electrical losses:","System losses:","Output Power:"];
            Plot_Power_Categories = 100*[sum(P_term)*panels+sum(P_below)*panels+sum(P_angle)*panels+sum(P_carnot)*panels+sum(P_emission)*panels-sum(P_gain)*panels;sum(P_cell)*panels+sum(P_metal)*panels+sum(P_ref)*panels+sum(P_diffA)*panels;sum(P_NRRV)*panels+sum(P_NRRI)*panels+sum(P_series)*panels+sum(P_shunt)*panels;sum(P_mismatch)*panels+sum(P_con)*panels+sum(P_cable)+sum(P_inv);sum(sum(P_out_AC))]/(sum(P_in)*panels);
            Plot_type = 2;
        end
    else
        if strcmp(type,'SHJ')|| strcmp(type,'BIF')|| strcmp(type,'T-F') %DC simulation for single junction
            Plot_Components = ["DC power";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Power = 100*[sum(P_out_DC(:,2));sum(P_mismatch);sum(P_con);sum(P_NRRI);sum(P_NRRV);sum(P_shunt);sum(P_series);sum(P_diffA);sum(P_ref);sum(P_metal);sum(P_cell);-1*sum(P_gain);sum(P_angle);sum(P_carnot);sum(P_emission);sum(P_below);sum(P_term)]/sum(P_in);
            Plot_Categories = ["Fundamental losses:", "Optical losses:", "Electrical losses:","System losses:","Output Power"];
            Plot_Power_Categories = 100*[sum(P_term)+sum(P_below)+sum(P_angle)+sum(P_carnot)+sum(P_emission)-sum(P_gain);sum(P_cell)+sum(P_metal)+sum(P_ref)+sum(P_diffA);sum(P_NRRV)+sum(P_NRRI)+sum(P_series)+sum(P_shunt);sum(P_mismatch)+sum(P_con);sum(P_out_DC(:,2))]/sum(P_in);
            Plot_type = 3;
        elseif strcmp(type,'Tan') || strcmp(type,'BIF-Tan') %DC simulation for tandem cell
            Plot_Components = ["DC power bottom cell";"DC power top cell";"Mismatch losses";"Cell interconnection";'Recombination I';'Recombination V';"Shunt resistance";"Series resistance";"Parasitic absorption";"Reflection";"Metal shading";'Cell spacing';"Non-ideality effect";'Angle mismatch';'Carnot losses';'Emission losses';'Below bandgap';'Thermalization'];
            Plot_Power = 100*[sum(P_out_DC(:,2));sum(P_out_DC(:,1));sum(P_mismatch);sum(P_con);sum(P_NRRI);sum(P_NRRV);sum(P_shunt);sum(P_series);sum(P_diffA);sum(P_ref);sum(P_metal);sum(P_cell);-1*sum(P_gain