function Vmax = compute_max_system_voltage(inv_type,system,string,module)
%COMPUTE_MAX_SYSTEM_VOLTAGE calculate the maximum voltage of the system
% This is needed to filter the suitable inverters from the database
%
% Parameters
% ----------
% inv_type : char
%   Type of inverter
% system : struct
%   Electrical parameters at system level
% string : struct
%   Electrical parameters at string level
% module : struct
%   Electrical parameters at module level
%
% Returns
% -------
% Vmax : double
%   Maximum voltage measured in the system


if strcmp(inv_type,'CEN')
    Vmax = max(system.voltage);
elseif strcmp(inv_type,'STR')
    if iscell(string.voltage)
        Vmax = max([string.voltage{:}],[],'all');
    else
        Vmax = max(string.voltage);
    end
elseif strcmp(inv_type,'MIC')
    if iscell(module.voltage)
        Vmax = max([module.voltage{:}],[],'all');
    else
        Vmax = max(module.voltage);
    end
end
end