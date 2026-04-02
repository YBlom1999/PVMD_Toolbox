function [P_cell,P_metal,P_ref,P_diffA,I_abs1,I_abs2] = OpticalLossesTandem(Opt_Losses_input,TOOLBOX_input,MODULE_output,CELL_output,CONSTANTS)
%OpticalLossesTandem Calculates the optical losses of the system 
% with tandem modules.
%
% This function calculates the optical losses given the incoming irradiance and
% optical properties of the solar module for with a tandem cell.
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
% I_abs1 : double
%   The absorbed current in the top cell
% I_abs2 : double
%   The absorbed current in the bottom cell
%
% Developed by Y. Blom


P_in = Opt_Losses_input.P_in;
P_fund = Opt_Losses_input.P_fund;
wav = Opt_Losses_input.wav;
photon_spec = Opt_Losses_input.photon_spec;
R = Opt_Losses_input.R;
A1 = Opt_Losses_input.A1;
A2 = Opt_Losses_input.A2;
A_diff = Opt_Losses_input.A_diff;
E_g1 = Opt_Losses_input.E_g1;
E_g2 = Opt_Losses_input.E_g2;
V_opt1 = Opt_Losses_input.V_opt1;
V_opt2 = Opt_Losses_input.V_opt2;
Angle_emit = Opt_Losses_input.Angle_emit;
Parameters1 = Opt_Losses_input.Parameters1;
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

Vth = k*T_cell/q;
Iph1 = Parameters1(:,1);
Rs1 = Parameters1(:,2);
Rsh1 = Parameters1(:,3);
n1 = Parameters1(:,4);
I01 = Parameters1(:,5);
if TOOLBOX_input.electric.Terminals == 2
    V_cell1 = Iph1.*Rsh1+I01.*Rsh1-I_mpp'.*(Rs1+Rsh1)-n1.*Vth'.*lambertw(Rsh1.*I01./(n1.*Vth').*exp(Rsh1./(n1.*Vth').*(Iph1+I01-I_mpp')));
    V_cell2 = V_mpp_cell'-V_cell1;
    V_cell1(isnan(V_cell1)) = 0;
    V_cell2(isnan(V_cell2)) = 0;
elseif TOOLBOX_input.electric.Terminals == 3
    V_cell1 = 2*V_mpp_cell/3;
    V_cell2 = V_mpp_cell/3;
elseif TOOLBOX_input.electric.Terminals ==4
    V_cell1 = V_mpp_cell(1);
    V_cell2 = V_mpp_cell(2);
end

%The wavelengths that are present at GENPRO are loaded
wav_GENPRO = CELL_output.CELL_FRONT.wav;
N = find(wav > wav_GENPRO(end)*1e-6,1)-1;  %N is the last wavelength that is in GENPRO
N_1 = find(wav > h*c/(q*E_g1),1)-1; %N1 is the wavelength corresponding to the bandgap of perovskite
N_2 = find(wav > h*c/(q*E_g2),1)-1; %N1 is the wavelength corresponding to the bandgap of silicon

A1 = interp1(wav_GENPRO,A1,wav(1:N)*1e6);
A2 = interp1(wav_GENPRO,A2,wav(1:N)*1e6);
R = interp1(wav_GENPRO,R,wav(1:N)*1e6);
A_diff = interp1(wav_GENPRO,A_diff,wav(1:N)*1e6);
A1(isnan(A1)) = 0;
A2(isnan(A2)) = 0;
R(isnan(R)) = 0;
A_diff(isnan(A_diff)) = 0;

%% Cell spacing
P_cell = (P_in-P_fund)*(1-N_cells*A_cell/A_mod);

%% Metalisation
P_metal = (P_in-P_fund)*N_cells*A_cell/A_mod*Sf;

%% Reflecion/Transmission and Absorption
photon_emission_actual1 = 2*Angle_emit*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt1))*h*c./wav-q*V_cell1')./(k*T_cell))-1);
photon_emission_actual2 = 2*Angle_emit*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt2))*h*c./wav-q*V_cell2')./(k*T_cell))-1);
P_ref_f1 = A_eff*(photon_spec(1:N,:)-photon_emission_actual1(1:N,:)).*R.*V_opt1*q;
P_ref_f2 = A_eff*(photon_spec(1:N,:)-photon_emission_actual2(1:N,:)).*R.*V_opt2*q;
P_ref = trapz(wav(1:N_1,:),P_ref_f1(1:N_1,:))+trapz(wav(N_1:N_2,:),P_ref_f2(N_1:N_2,:));

P_diffA_f1 = A_eff*(photon_spec(1:N,:)-photon_emission_actual1(1:N,:)).*A_diff.*V_opt1*q;
P_diffA_f2 = A_eff*(photon_spec(1:N,:)-photon_emission_actual2(1:N,:)).*A_diff.*V_opt2*q;
P_diffA = trapz(wav(1:N_1,:),P_diffA_f1(1:N_1,:))+trapz(wav(N_1:N_2,:),P_diffA_f2(N_1:N_2,:));

I_abs1_f = A_eff/N_cells*(photon_spec(1:N,:)-photon_emission_actual1(1:N,:)).*A1.*q;
I_abs2_f = A_eff/N_cells*(photon_spec(1:N,:)-photon_em