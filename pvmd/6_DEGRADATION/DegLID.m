function [k_pred,factor] = DegLID(ELECTRIC_output,Jph_abs,factor_max,C,Time,CONSTANTS)
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
T_STC = CONSTANTS.T_STC;

Voltage = 0:0.01:0.9;
Vth=k_b*T_STC/q;

Iph = ELECTRIC_output.Parameters_STC_1(1,1);
Rs = ELECTRIC_output.Parameters_STC_1(1,2);
Rsh = ELECTRIC_output.Parameters_STC_1(1,3);
n = ELECTRIC_output.Parameters_STC_1(1,4);
I0_init = ELECTRIC_output.Parameters_STC_1(1,5);


Pmpp = zeros(1,length(Time));
k_pred = zeros(1,length(Time));

for i = 1:length(Time)
    factor_tau = factor_max+(1-factor_max)*exp(-sum(Jph_abs(1:Time(i)))/C);
    I0 = I0_init/factor_tau;

    z=(Rs*I0/(n*Vth*(1+Rs/Rsh)))*exp((Rs*(Iph+I0)+Voltage)./(n*Vth*(1+Rs/Rsh)));
    Current=(Iph+I0-Voltage/(Rsh))/(1+Rs/Rsh)-lambertw(z).*(n*Vth)/Rs;
    
    Pmpp(i) = max(Voltage.*Current);
    
    if i > 1
        k_pred(i) = (1-Pmpp(i)/Pmpp(i-1))/(Time(i)-Time(i-1));
    end
end

factor = log10(I0)/log10(I0_init);
end