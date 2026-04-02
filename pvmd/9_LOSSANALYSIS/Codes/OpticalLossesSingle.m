function [P_cell,P_metal,P_ref,P_diffA,I_abs] = OpticalLossesSingle(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS)
%OpticalLossesSingle Calculates the optical losses of the system 
% with single junction modules.
%
% This function calculates the optical losses given the incoming irradiance and
% optical properties of the solar module for a single junction.
%
% Parameters
% ----------
% Opt_Losses_input : struc
%   The needed inputs for the fundamental losses calculation
% MODULE_output : struct
%   Simulation results of the MODULE module
% CELL_output : struct
%   Simulation results of the CELL module
% CONSTANTS : struct
%   Physical constants
%
% Returns
% -------
% P_cell : double
%   The losses due to cell spacing
% P_metal : double
%   The losses due to metal shading
% P_ref : double
%   The losses due to reflection
% P_diffA : double
%   The losses due to the Parasitic absorption
% I_abs : double
%   The absorbed current
%
% Developed by Y. Blom

P_in = Opt_Losses_input.P_in;
P_fund = Opt_Losses_input.P_fund;
wav = Opt_Losses_input.wav;
photon_spec = Opt_Losses_input.photon_spec;
R = Opt_Losses_input.R;
A = Opt_Losses_input.A;
A_diff = Opt_Losses_input.A_diff;
E_g = Opt_Losses_input.E_g;
V_opt1 = Opt_Losses_input.V_opt1;
Angle_emit = Opt_Losses_input.Angle_emit;
V_mpp = Opt_Losses_input.V_mpp;
I_mpp = Opt_Losses_input.I_mpp;
T_cell = Opt_Losses_input.T_cell;

% constants
h = CONSTANTS.h;
q = CONSTANTS.q;
c = CONSTANTS.c;
k = CONSTANTS.k;

%Properties of the module
N_cells = MODULE_output.N;
A_cell = MODULE_output.A;
A_mod = MODULE_output.Amod;
Sf = TOOLBOX_input.electric.shading/100;
A_eff = N_cells*A_cell*(1-Sf);

%Electrical properties
R_con = TOOLBOX_input.electric.resistance;
V_mpp_cell = V_mpp/N_cells+I_mpp*R_con;

%The wavelengths that are present at GENPRO are loaded
wav_GENPRO = CELL_output.CELL_FRONT.wav;
angles_GENPRO = CELL_output.CELL_FRONT.aoi;
N_1 = find(wav > h*c/(q*E_g),1)-1; %N1 is the wavelength corresponding to the bandgap of silicon
N = find(wav > wav_GENPRO(end)*1e-6,1)-1;  %N is the last wavelength that is in GENPRO

A = interp1(wav_GENPRO,A,wav(1:N)*1e6);
R = interp1(wav_GENPRO,R,wav(1:N)*1e6);
A_diff = interp1(wav_GENPRO,A_diff,wav(1:N)*1e6);
A(isnan(A)) = 0;
R(isnan(R)) = 0;
A_diff(isnan(A_diff)) = 0;
%% Cell spacing
P_cell = (P_in-P_fund)*(1-N_cells*A_cell/A_mod);

%% Metalisation
P_metal = (P_in-P_fund)*N_cells*A_cell/A_mod*Sf;

%% Reflecion/Transmission and Absorption
photon_emission_actual = 2*Angle_emit*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt1))*h*c./wav-q*V_mpp_cell)./(k*T_cell))-1);
P_ref_f = A_eff*(photon_spec(1:N,:)-photon_emission_actual(1:N,:)).*R.*V_opt1*q;
P_diffA_f = A_eff*(photon_spec(1:N,:)-photon_emission_actual(1:N,:)).*A_diff.*V_opt1*q;
I_abs_f = A_eff/N_cells*(photon_spec(1:N,:)-photon_emission_actual(1:N,:)).*A*q;
P_ref = trapz(wav(1:N_1,:),P_ref_f(1:N_1,:));
P_diffA = 