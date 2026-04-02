function [k_pred,Rcon] = DegThermalCycling(ELECTRIC_output,T,C,Ea,n,b,Time,Rcon_init,NumCells,CONSTANTS)
%DegDiscoloration Calculates the degradation due to LID
%
% This function calculates the degradation rate caused by LID. A decrease
% in minority carrier lifetime is simulated, which has an effect on the
% equivalent circuit parameters
%
% Parameters
% ----------
% ELECTRIC_output: struct
%   Result of the electrical simulation
% Jph_abs : double
%   Absorbed current of the module
% factor_max : double
%   The maximum value the saturation current should be multiplied with.
% C : double
%   The constant that determines how fast the degradation goes.
% Time: double
%   The time for which the degradation needs to be simulated
% CONSTANTS : struct
%   Structure of physical constants
%
% Returns
% -------
% k_pred: double
%   The predicted degradation rate due to LID
%
% Developed by by Youri Blom
k_b = CONSTANTS.k_b;
q = CONSTANTS.q;

T_daily = reshape(T(1:8760),[24,8760/24]);
T_max = mean(max(T_daily));
Delta_T = mean(max(T_daily)-min(T_daily));

TC = T_max-Delta_T/2;

T_diff = (T-TC)>0;
crosses = abs(T_diff(2:end)-T_diff(1:end-1));
Rate = [0,cumsum(crosses)];

Damage = C*(Delta_T)^n*Rate.^b*exp(-q*Ea/k_b/T_max);
k_pred = [0,diff(Damage)];

P_STC = ELECTRIC_output.P_STC;
Impp_STC = ELECTRIC_output.Impp_STC;
P_loss = P_STC*sum(k_pred);

Rcon = Rcon_init+P_loss/Impp_STC.^2/NumCells;

end