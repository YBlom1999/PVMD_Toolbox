function [P_term,P_below,P_emission,P_carnot,P_angle,P_gain] = FundamentalLossesTandem(Fund_Losses_Input,TOOLBOX_input,CELL_output,MODULE_output,CONSTANTS)
%FundamentalLossesSingle Calculates the fundamental losses of the system 
% with tandem modules.
%
% This function calculates the fundamental losses given the incoming irradiance and
% fundamental properties of the solar module for a tandem cell.
%
% Parameters
% ----------
% Fund_Losses_input : struc
%   The needed inputs for the fundamental losses calculation
% TOOLBOX_input : struct
%   Simulation parameters
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
% V_opt1 : double
%   The optimal voltage of cell 1
% V_opt2 : double
%   The optimal voltage of cell 2
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
A1 = Fund_Losses_Input.A1;
A2 = Fund_Losses_Input.A2;
E_g1 = Fund_Losses_Input.E_g1;
E_g2 = Fund_Losses_Input.E_g2;
Parameters1 = Fund_Losses_Input.Parameters1;
I_mpp = Fund_Losses_Input.I_mpp;
V_mpp = Fund_Losses_Input.V_mpp;
T_cell = Fund_Losses_Input.T_cell;
V_opt1 = Fund_Losses_Input.V_opt1;
V_opt2 = Fund_Losses_Input.V_opt2;


%% constants
h = CONSTANTS.h;
q = CONSTANTS.q;
c = CONSTANTS.c;
k = CONSTANTS.k;
T_S = CONSTANTS.T_S;

%The wavelengths that are present at GENPRO are loaded
wav_GENPRO = CELL_output.CELL_FRONT.wav;
N = find(wav > wav_GENPRO(end)*1e-6,1)-1;  %N is the last wavelength that is in GENPRO
N_1 = find(wav > h*c/(q*E_g1),1)-1; %N1 is the wavelength corresponding to the bandgap of perovskite
N_2 = find(wav > h*c/(q*E_g2),1)-1; %N1 is the wavelength corresponding to the bandgap of silicon

%Properties of the module
N_cells = MODULE_output.N;
A_mod = MODULE_output.Amod;


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
    V_cell1(isinf(V_cell1)) = 0;
    
    V_cell2(isnan(V_cell2)) = 0;
    V_cell2(isinf(V_cell2)) = 0;
elseif TOOLBOX_input.electric.Terminals == 3
    V_cell1 = 2*V_mpp_cell/3;
    V_cell2 = V_mpp_cell/3;
elseif TOOLBOX_input.electric.Terminals ==4
    V_cell1 = V_mpp_cell(1);
    V_cell2 = V_mpp_cell(2);
end

V_opt1_ideal = E_g1*(1-T_cell/T_S)-k*T_cell/q*log(Angle_emit_Vopt/Angle_abs);
V_opt2_ideal = E_g2*(1-T_cell/T_S)-k*T_cell/q*log(Angle_emit_Vopt/Angle_abs);
factor_Vopt1 = (E_g1-V_opt1)./(E_g1-V_opt1_ideal);
factor_Vopt2 = (E_g2-V_opt2)./(E_g2-V_opt2_ideal);

%% Thermalization & Emission

%A numberical integration is done to calculate the thermalization losses
photon_emission1 = 2*Angle_emit_emission*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt1))*h*c./wav-q*V_opt1)./(k*T_cell))-1);
photon_emission2 = 2*Angle_emit_emission*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt2))*h*c./wav-q*V_opt2)./(k*T_cell))-1);

P_term_f1 = A_mod*photon_spec.*(h*c./wav-E_g1*q);
P_term_f2 = A_mod*photon_spec.*(h*c./wav-E_g2*q);
P_term = trapz(wav(1:N_1),P_term_f1(1:N_1,:))+trapz(wav(N_1:N_2),P_term_f2(N_1:1:N_2,:));

P_emission1_f = E_g1*q*A_mod*photon_emission1;
P_emission2_f = E_g2*q*A_mod*photon_emission2;
P_emission = trapz(wav(1:N_1),P_emission1_f(1:N_1,:))+trapz(wav(N_1:N_2),P_emission2_f(N_1:N_2,:));

I_max1_f = A_mod*(photon_spec-photon_emission1)*q;
I_max2_f = A_mod*(photon_spec-photon_emission2)*q;
I_max1 = trapz(wav(1:N_1),I_max1_f(1:N_1,:));
I_max2 = trapz(wav(N_1:N_2),I_max2_f(N_1:N_2,:));


%% Below bandgap losses
%A numberical integration is done to calculate the below bandgap lossses.
P_below_f = A_mod*Irr_spec;
P_below = trapz(wav(N_2:end),P_below_f(N_2:end,:));

%% Carnot losses
V_carnot1 = E_g1*T_cell/T_S.*factor_Vopt1;
V_carnot2 = E_g2*T_cell/T_S.*factor_Vopt2;
P_carnot = V_carnot1.*I_max1 + V_carnot2.*I_max2;

%% Angle mismatch losses
V_angle1 = k*T_cell/q*log(Angle_emit_Vopt/Angle_abs).*factor_Vopt1;
V_angle2 = k*T_cell/q*log(Angle_emit_Vopt/Angle_abs).*factor_Vopt2;
P_angle = V_angle1.*I_max1+V_angle2.*I_max2;


%% Gains
A1 = interp1(wav_GENPRO,A1,wav(1:N)*1e6);
A1(isnan(A1)) = 0;
A2 = interp1(wav_GENPRO,A2,wav(1:N)*1e6);
A2(isnan(A2)) = 0;

photon_emission_actual1 = 2*Angle_emit_emission*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt1))*h*c./wav-q*V_cell1')./(k*T_cell))-1);
photon_emission_actual2 = 2*Angle_emit_emission*c./(wav.^4)*1./(exp((ones(length(wav),length(V_opt1))*h*c./wav-q*V_cell2')./(k*T_cell))-1);
I_gain1_f1 = A_mod*(photon_emission1-photon_emission_actual1)*q;
I_gain1_f2 = A_mod*(photon_spec(1:N,:)-photon_emission_actual1(1:N,:)).*A1*q;
I_gain1_f3 = A_mod*(photon_spec(1:N,:)-photon_emission_actual1(1:N,:)).*A2*q;
I_gain1 = trapz(wav(1:N_1),I_gain1_f1(1:N_1,:))-trapz(wav(1:N_1),I_gain1_f3(1:N_1,:))+trapz(wav(N_1:N),I_gain1_f2(N_1:N,:));

I_gain2_f1 = A_mod*(photon_emission2-photon_emission_actual2)*q;
I_gain2_f2 = A_mod*(photon_spec(1:N,:)-photon_emission_actual2(1:N,:)).*A2*q;
I_gain2_f3 = A_mod*(photon_spec(1:N,:)-photon_emission_actual2(1:N,:)).*A1*q;
I_gain2 = trapz(wav(N_