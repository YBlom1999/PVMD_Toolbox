function [CONVERSION_output,TOOLBOX_input] = CONVERSION_MAIN(TOOLBOX_input,ELECTRIC_output)
%CONVERSION_MAIN Block that converts the DC power to AC
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
% ELECTRIC_output : struct
%   Results from the previous electric block
%
% Returns
% -------
% conversion_out : struct
%   Main results from the conversion block

%---- Build useful structures for the code
[module, inverter, cable_losses, simulation, power_optimizer] = ...
    get_conversion_inputs(TOOLBOX_input,ELECTRIC_output);

%---- Modify the toolbox input
if ~TOOLBOX_input.script
    TOOLBOX_input = modify_toolbox_input(TOOLBOX_input, inverter, ...
        power_optimizer, cable_losses);
end

%---- Code for power optimizers
if strcmp(inverter.type,'POW-OPT')
    CONVERSION_output = power_optimizer_main(module, inverter, cable_losses, ...
        simulation, power_optimizer);
    return;
end

%---- Parameters precalculation for the main code
[system, string] = compute_pv_power(module,inverter,simulation.periodic);

Vmax = compute_max_system_voltage(inverter.type,system,string,module);

inv_constants = get_inverter_constants(Vmax,inverter);

%---- DC to AC conversion including cable losses
CONVERSION_output = dc_ac_conversion(inverter, system, string, module, ...
    cable_losses, inv_constants, simulation.periodic);

%---- Plot the results if desired
if inverter.plot_graphs
    plot_inverter_graphs(CONVERSION_output.Pac, ...
        CONVERSION_output.Pdc, CONVERSION_output.eff, ...
        inverter.type, simulation.periodic)
end

end

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

function inv_constants = get_inverter_constants(Vdc_max,inverter)
%GET_INVERTER_CONSTANTS Return SNL constants of the inverter selected
% Read the inverters characteristics on SNL database and select those that
% are compatible: if max array voltage does not exceed max inverter input
% voltage. Compare that selected database with the inverter chosen by the
% user
%
% Parameters
% ----------
% Vdc_max : double
%   Maximum DC voltage in the time period
% inverter : struct
%   User inputs related to the inverter characteristics
%
% Returns
% -------
% inv_constants : double
%   SNL constants of the chosen inverter if compatible

% Import data from SNL database
data_folder = get_data_path('conv');
file = fullfile(data_folder,'constants.xlsx');
[all_constants, all_models] = xlsread(file);
all_models = all_models(4:end,1);

%Remove incompatible inverter types
compatible = all_constants(:,13) > Vdc_max;
models = all_models(compatible);
constants = all_constants(compatible,:);

model_idx = find(strcmp(models,inverter.model),1);
if isempty(find(strcmp(all_models,inverter.model),1))
    error("conversion:InvalidInverterModel", "The selected inverter model %s is not in the list",...
        inverter.model);
elseif isempty(model_idx)
    error("conversion:ExceedMaximumVoltage", ['The maximum PV system input voltage of %0.1f V exceeds the ',...
        'maximum voltage of the selected inverter model %s.'],...
        Vdc_max, inverter.model);
else
    inv_constants = constants(model_idx,:);
end

end