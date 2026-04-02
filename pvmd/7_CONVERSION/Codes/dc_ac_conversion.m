function CONVERSION_output = dc_ac_conversion(inverter, system, string, ...
    module, cable_losses, inv_constants, is_periodic)
%DC_AC_CONVERSION Convert DC power to AC including cable losses
%
% Parameters
% ----------
% inverter : struct
%   User inputs related to the inverter characteristics
% system : struct
%   Electrical parameters at system level
% string : struct
%   Electrical parameters at string level
% module : struct
%   Electrical parameters at module level
% cable_losses : struct
%   User inputs related to the cable losses
% inv_constants : double
%   SNL constants of the chosen inverter
% is_periodic : logical
%   Indicates whether the simulation is periodic or not
%
% Returns
% -------
% CONVERSION_output : struct
%   Main results from the conversion block

if strcmp(inverter.type,'CEN')
    system.power = include_cable_losses(system.power,string.current, ...
        system.current,cable_losses,'CEN');
    system.power_STC = include_cable_losses(system.power_STC,string.current_STC, ...
        system.current_STC,cable_losses,'CEN');

    system.voltage = get_voltage_from_power_current(system.power,system.current);
    system.voltage_STC = get_voltage_from_power_current(system.power_STC,system.current_STC);

    
    [Pac,Pac_tot,eff,eff_tot,Pac_STC] = AC_conversion(system,inv_constants,...
        inverter,is_periodic);
elseif strcmp(inverter.type,'STR')
    string.power = include_cable_losses(string.power,string.current, ...
        nan,cable_losses,'STR');
    string.voltage = get_voltage_from_power_current(string.power,string.current);
    
    [Pac,Pac_tot,eff,eff_tot,Pac_STC] = AC_conversion(string,inv_constants,...
        inverter,is_periodic);
    if is_periodic
        system.power = string.power*inverter.num_strings;
        Pac = Pac*inverter.num_strings;
        Pac_tot = Pac_tot*inverter.num_strings;
    else
        system.power = string.power;
    end
elseif strcmp(inverter.type,'MIC')
    module.power = include_cable_losses(module.power,module.current, ...
        nan,cable_losses,'MIC');
    module.voltage = get_voltage_from_power_current(module.power,module.current);
    
    [Pac,Pac_tot,eff,eff_tot,Pac_STC] = AC_conversion(module,inv_constants,...
        inverter,is_periodic);
    if is_periodic
        system.power = module.power*inverter.num_modules;
        Pac = Pac*inverter.num_modules;
        Pac_tot = Pac_tot*inverter.num_modules;
    else
        system.power = module.power;
    end
end

% Print the output results
fprintf('The predicted inverter output energy yield is %0.2f Wh [AC].\n',...
    sum(Pac_tot));
fprintf('The overall inverter efficiency is %0.2f %%\n',eff_tot)

% Build the output structure
CONVERSION_output.Pac_total = Pac_tot;
CONVERSION_output.Pac = Pac;
CONVERSION_output.Pdc = system.power;
CONVERSION_output.eff = eff;
CONVERSION_output.eff_overall = eff_tot;
CONVERSION_output.Pac_STC = Pac_STC;

end


function voltage = get_voltage_from_power_current(power, current)
%GET_VOLTAGE_FROM_POWER_CURRENT Apply P = V*I to obtain V
%
% Parameters
% ----------
% power : cell/array
%   electric power
% current: cell/array
%   electric current
%
% Returns
% -------
% voltage: cell/array
%   electric voltage

if iscell(power)
    voltage = cellfun(@(x,y) x./y,power,current,'UniformOutput',false);
    voltage = cellfun(@(x) fillmissing(x,'constant',0), voltage,...
        'UniformOutput', false);
else
    voltage = power./current;
    voltage(isnan(voltage)) = 0;
end
end