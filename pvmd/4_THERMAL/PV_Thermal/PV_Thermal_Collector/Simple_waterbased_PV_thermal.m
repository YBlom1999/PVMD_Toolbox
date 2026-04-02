function [pvt_collector_output] = Simple_waterbased_PV_thermal(MODULE_output, WEATHER_output)
% Calculates the thermal and PV yield of simple_water-based PVT collector
% (Rd)
%
% Parameters
% ----------
% MODULE_output : struct
%   Output of this module block
% WEATHER_output : struct
%   Simulation results of the WEATHER module
%
% Returns
% -------
% pvt_collector_output : struct
%   Simulation results of the PVT module
% 
% Developed by ZUA.

I = mean(WEATHER_output.Irr,2)';
Tam = [WEATHER_output.ambient_temperature]';
L = MODULE_output.ML;
W = MODULE_output.MW;
% Vw = [WEATHER_output.wind_speed]';
% Tc = zeros(size(Tam));

idx = find(I~=0);                          % find the indices of non-zero elements in G
if length(idx) > length(I)                 % if there are non-zero elements in G
    idx = idx(1:length(I));                % select the non-zero elements
end
I1 = I(idx);                               % select the corresponding elements from G
T_win = Tam + 1;                           % Temperature of the fluid entering the collector (*C)

eta_E_ref = 0.15;                          % Electrical efficiency at standard test conditions (tau_eff=>92.5, eta_el 15.52(cel) to 12.97(col).)
tau_cel = 0.90;                            % transperrency of PV module
alpha_cel = 0.95;                          % Absorptivity of cell
delta_cel = 0.9;                           % packing factor
A_eff = tau_cel*alpha_cel*delta_cel;       % Effective absorption factor (total absorption factor(0.85-0.93)-eta_el)

w = 0.5*W;                                 % Receiver width
dm_w = 0.025;                              % Water flow rate
C_w = 4186;                                % Specific heat of water
beta = 0.002;                              % Temperature coefficient (0.002 & 0.0045 - depending on type....)

h_cs = 120;                                % Heat transfer coefficient cell to sheet -   W/(m^2.K)
h_st = 200;                                % Heat transfer coefficient sheet to tube -   W/(m^2.K) 
h_tw = 700;                                % Heat transfer coefficient tube to fluid -   W/(m^2.K) 
h_cw = 1/((1/h_cs) + (1/h_st) + (1/h_tw)); % Heat transfer coefficient cell to water -   Range: 40 to 500 W/(m^2.K) => cts:120,stt:200,ttf:700 => 70 turbulent, 40 laminar
h_ca = 25;                                 % Heat transfer coefficient cell to ambient - Range: 5 to 25 W/(m^2.K)

% % % h_v = 2.8 + 3.*Vw;
% % % % h_v= 5.7+3.8.*Vw;                  % Convective heat coefficient due to wind (Vw<5)
% % % % h_v(Vw>=5) = 6.47+Vw.^0.78;          
% % % F = 0.93;                            % Sky view factor
% % % em_c = 0.90;                         % far infrared emissivity of cell
% % % sigma = 5.67e-8;
% % % Tsk = 0.0552*(Tam+273.15).^1.5;      % Temperature of Sky
% % % T_s = Tsk - 273.15;
% % % h_r = F*em_c*sigma*((Tc+T_s)./2).^3;
% % % h_ca = h_v + h_r;

T_stag = Tam + ((A_eff*I*(W/w))/h_ca);
T_red = (T_win - Tam)./I;
l = (dm_w*C_w/w)*((h_ca + h_cw)./(h_ca*h_cw));
Two = T_stag + (T_win - T_stag).*exp(-L./l);
Tmean = (T_win + Two) ./ 2;
Tc = (Tam.*h_ca + A_eff.*I.*(W/w) + Two.*h_cw)./(h_ca + h_cw);
eta_E = zeros(size(I));                            % Preallocate eta_E with zeros
eta_E(idx) = eta_E_ref*(1 - beta*(Tc(idx) - 25));  % Calculate eta_E only for the valid indices

eta_tha = ((A_eff./h_ca) - T_red.*w./W).*((dm_w.*C_w)/(w.*L)).*(1-exp(-L./l));
eta_tha(h_cw>10.*h_ca) = A_eff - h_ca*(w*W).*T_red(h_cw>10.*h_ca);
eta_tha(eta_tha<0) = 0;

% Total solar irradiance received by the PVT system in kWh/m^2/year
total_I_sun = sum(mean(WEATHER_output.Irr,2))/1000; 
% Average thermal efficiency
etaSum = sum(eta_tha)/length(I1);
% Average thermal efficiency
etaE = sum(eta_E)/length(I1);

%% change if kWh or kWh/m2
% Yearly total thermal output in kWh/m2
total_thermal_output = etaSum*total_I_sun;
% Yearly total electrical output in kWh/m2
total_electrical_output = etaE*total_I_sun;

pvt_collector_output.eta_tha=eta_tha;
pvt_collector_output.etaSum=etaSum;
pvt_collector_output.eta_E=eta_E;
pvt_collector_output.etaE=etaE;
pvt_collector_output.Two=Two;
pvt_collector_output.T=Tc;
pvt_collector_output.Tmean=Tmean;
pvt_collector_output.total_thermal_output=total_thermal_output;
pvt_collector_output.total_electrical_output=total_electrical_output;
