function [k_pred,Rcon] = DegThermalCycling(ELECTRIC_output,T,A_c,Ea,n,b,Rcon_init,NumCells,CONSTANTS)
% DegThermalCycling Calculates the degradation due to thermal cycling (TC)
%
% This function calculates the degradation rate caused by TC. Due to
% differences in thermal expansion coefficients, the solder joints
% experience mechanical stress.
%
% Parameters
% ----------
% ELECTRIC_output: struct
%   Result of the electrical simulation
% T : double
%   Temperature of the module
% A_c : double
%   The pre-exponential constant of thermal cycling
% Ea : double
%   The activation energy for thermal cycling
% n : double
%   The exponent related to delta T
% b : double
%   The exponent related to the ammount of thermal cycles
% Rcon_init: double
%   The original resistance for the interconnection
% NumCells: double
%   The number of cells in the module
% CONSTANTS: struct
%   physical constants
%
% Returns
% -------
% k_pred: double
%   The predicted degradation rate due to LID
% Rcon: double
%   The updated resistance of the interconnection
%
% Developed by by Youri Blom
k = CONSTANTS.k;
q = CONSTANTS.q;

T_daily = reshape(T(1:8760),[24,8760/24]);
T_max = mean(max(T_daily));
Delta_T = mean(max(T_daily)-min(T_daily));

TC = T_max-Delta_T/2;

T_diff = (T-TC)>0;
crosses = abs(T_diff(2:end)-T_diff(1:end-1));
Rate = [0,cumsum(crosses)];

Damage = A_c*(Delta_T)^n*Rate.^b*exp(-q*Ea/k/T_max);
k_pred = [0,diff(Damage)];

P_STC = ELECTRIC_output.P_STC;
Impp_STC = ELECTRIC_output.Impp_STC;
P_loss = P_STC*sum(k_pred);

Rcon = Rcon_init+P_loss/Impp_STC.^2/NumCells;

end