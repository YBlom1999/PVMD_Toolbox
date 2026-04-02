function [P_term,P_below,P_emission,P_carnot,P_angle,P_gain] = FundamentalLossesSingle(Fund_Losses_Input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS)
%FundamentalLossesSingle Calculates the fundamental losses of the system 
% with single junction modules.
%
% This function calculates the fundamental losses given the incoming irradiance and
% fundamental properties of the solar module for a single junction.
%
% Parameters
% ----------
% Fund_Losses_input : struc
%   The needed inputs for the fundamental losses calculation
% CELL_output : struct
%   Simulation results of the CELL module
% MODULE_output : struct
%   Simulation results of the MODULE module
% CONSTANTS : struct
%   Physical constants
%
% Returns
% -------
% P_term : double
%   The losses due to thermalization
% P_below : double
%   The losses due to below bandgap non absorption
% P_emission : double
%   The losses due to emission
% P_carnot : double
%   The losses due to the Carnot limit
% P_angle : double
%   The losses due to an angle mismatch
% P_gain : double
%   The gain in power due to non-idealities
% V_mpp_cell1 : double
%   The maximum power point voltage of cell 1
% V_mpp_cell2 : double
%   The maximum power point voltage of cell 2
%
% Developed by Y. Blom

wav = Fund_Losses_Input.wav;
Irr_spec = Fund_Losses_Input.Irr_spec;
photon_spec = Fund_Losses_Input.photon_spec;
Angle_abs = Fund_Losses_Input.Angle_abs;
Angle_emit_Vopt = Fund_Losses_Input.Angle_emit_Vopt;
Angle_emit_emission = Fund_Losses_Input.Angle_emit_emission;
A = Fund_Losses_Input.A;
E_g = Fund_Losses_Input.E_g;
I_mpp = Fund_Losses_Input.I_mpp;
V_mpp = Fund_Losses_Input.V_mpp;
T_cell = Fund_Losses_Input.T_cell;
V_opt1 = Fund_Losses_Input.V_opt1;


%% constants
h = CONSTANTS.h;
q = CONSTANTS.q;
c = CONSTANTS.c;
k = CONSTANTS.k;
T_S = CONSTANTS.T_S;

%Properties module
A_mod = MODULE_output.Amod;
N_cells = MODULE_output.N;
R_con = TOOLBOX_input.electric.resistance;
V_mpp_cell1 = V_mpp/N_cells+I_mpp*R_con;
V_mpp_cell2 = 0;

V_opt1_ideal = E_g*(1-T_cell/T_S)-k*T_cell/q*log(Angle_emit_Vopt/Angle_abs);
V_opt2 = 0;

factor_Vopt = (E_g-V_opt1)./(E_g-V_opt1_ideal);

%The wavelengths that are present at GENPRO are loaded
wav_GENPRO = CELL_output.CELL_FRONT.wav;
N = find(wav > wav_GENPRO(end)*1e-6,1)-1;  %N is the last wavelength that is in GENPRO
N_1 = find(wav > h*c/(q*E_g),1)-1; %N1 is the wavelength corresponding to the bandgap of silicon

%% Thermalization & Emission

%A numberical integration is done to calculate the thermalization losses
photon_emission = 2*Angle_emit_emission*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt1))*h*c./wav-q*V_opt1)./(k*T_cell))-1);
P_term_f = A_mod*photon_spec.*(h*c./wav-E_g*q);
P_emission_f = E_g*q*A_mod*photon_emission;
I_max_f = A_mod*(photon_spec-photon_emission)*q;
P_term = trapz(wav(1:N_1),P_term_f(1:N_1,:));
P_emission = trapz(wav(1:N_1),P_emission_f(1:N_1,:));
I_max = trapz(wav(1:N_1),I_max_f(1:N_1,:));

%% Below bandgap losses
%A numberical integration is done to calculate the below bandgap lossses.
P_below_f = A_mod*Irr_spec;
P_below = trapz(wav(N_1:end),P_below_f(N_1:end,:));


%% Carnot losses
V_carnot = E_g*T_cell/T_S.*factor_Vopt;
P_carnot = V_carnot.*I_max;

%% Angle mismatch losses
V_angle = k*T_cell/q*log(Angle_emit_Vopt/Angle_abs).*factor_Vopt;
P_angle = V_angle.*I_max;


%% Gains
A = interp1(wav_GENPRO,A,wav(1:N)*1e6);
A(isnan(A)) = 0;

photon_emission_actual = 2*Angle_emit_emission*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt1))*h*c./wav-q*V_mpp_cell1)./(k*T_cell))-1);
I_gain_f1 = A_mod*(photon_emission-photon_emission_actual)*q;
I_gain_f2 = A_mod*(photon_spec(1:N,:)-photon_emission_actual(1:N,:)).*A*q;
I_gain = trapz(wav(1