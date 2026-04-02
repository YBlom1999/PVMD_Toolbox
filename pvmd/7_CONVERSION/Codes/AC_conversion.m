function [Pac,Pac_total,eff,eff_overall,Pac_STC] = AC_conversion(system,...
    inv_constants,inverter,run_periodic)
%AC_CONVERSION Conversion of inverter output from DC to AC using SNL model
%
% Parameters
% ----------
% system : struct
%   DC electrical parameters (P,V,I) arriving at the inverter
% inv_constants : double
%   SNL constants of the user-selected inverter
% inverter : struct
%   User inputs related to the inverter characteristics
% run_periodic : logical
%   Indicate if the simulations are periodic or not
%
% Returns
% -------
% Pac : double
%   AC power out of the inverter
% Pac_total : double
%   Total AC power out of the inverter
% eff : double
%   Efficiency of the inverter
% eff_overall : double
%   Overall efficiency of the inverter

if run_periodic || strcmp(inverter.type,'CEN')
    [Pac,eff,eff_overall] = SNL(inv_constants,system.power,system.voltage);
    [Pac_STC,~,~] = SNL(inv_constants,system.power_STC,system.voltage_STC);
    Pac_total = sum(Pac);
else
    [Pac,eff,eff_overall] = cellfun(@(x,y) SNL(inv_constants,x,y),...
        system.power,system.voltage,'UniformOutput',false);
    [Pac_STC,~,~] = cellfun(@(x,y) SNL(inv_constants,x,y),...
        system.power_STC,system.voltage_STC,'UniformOutput',false);
    Pac_total = sum(cell2mat(Pac));
    eff_overall = mean(cell2mat(eff_overall));
end
end

function [Pac,eff,eff_overall] = SNL(constants,Pdc,Vdc)
%SNL Apply SNL inverter model to obtain the AC power and efficiency
%
% Parameters
% ----------
% constants : double
%   SNL constants of the user-selected inverter
% Pdc : double
%   DC power arriving at the inverter
% Vdc : double
%   DC voltage arriving at the inverter
%
% Returns
% -------
% Pac : double
%   AC power out of the inverter
% eff : double
%   efficiency of the inverter
% eff_overall : double
%   overall efficiency of the inverter
%
% Author: Tim Stark?

Paco = constants(4);
Pdco = constants(5);
Vdco = constants(6);
Pso = constants(7); % DC power required to start the inversion process
Pnt = constants(12); % DC power loss during night
    
C0 = constants(8); 
C1 = constants(9); 
C2 = constants(10); 
C3 = constants(11);

A = Pdco*(1+C1*(Vdc-Vdco));
B = Pso*(1+C2*(Vdc-Vdco));
C = C0*(1+C3*(Vdc-Vdco));

Pac = ((Paco./(A-B))-C.*(A-B)).*(Pdc-B) + C.*(Pdc-B).^2;        

% Inverter clipping at maximum rated AC Power
Pac(Pac>Paco) = Paco;
% Inverter night tare losses
Pac(Pdc<Pso)= -1*abs(Pnt);

eff = Pac./Pdc*100;
eff(eff<0) = 0;

eff_overall = sum(Pac)/sum(Pdc)*100;
end