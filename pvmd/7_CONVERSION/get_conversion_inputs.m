function [module, inverter, cable_losses, simulation, power_optimizer] = ...
    get_conversion_inputs(TOOLBOX_input,ELECTRIC_output)
%GET_CONVERSION_INPUT Prepare the inputs to the main conversion function
% Finds the inputs to the conversion block, either by asking the user
% through the GUI or by preparing the TOOLBOX_input variable. Should be
% temporal until refactoring of the toolbox is complete
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
% module : struct
%   Electrical parameters at module level
% inverter : struct
%   User inputs related to the inverter characteristics
% cable_losses : struct
%   User inputs related to the cable losses
% simulation : struct
%   Details of the simulation
% power_optimizer : struct
%   User inputs related to the power optimizer characteristics
%
% Author: A. Alcaniz

if ~TOOLBOX_input.runPeriodic
    nummod = TOOLBOX_input.Scene.N_panels;
else
    nummod = nan;
end


inv_type = determine_inverter_type(TOOLBOX_input.Conversion.ConversionType);
if strcmp(inv_type,'POW-OPT')
    pow_opt_model = TOOLBOX_input.Conversion.PowerOpt.Model;
    add_central_inv = TOOLBOX_input.Conversion.PowerOpt.AddCentralInverter;
    if add_central_inv
        inv_model = TOOLBOX_input.Conversion.Model;
        fixed_voltage = TOOLBOX_input.Conversion.PowerOpt.FixedVoltage;
    else
        fixed_voltage = nan;
    end
else
    inv_model = TOOLBOX_input.Conversion.Model;
    add_central_inv = true;
end

if add_central_inv
    parallel = TOOLBOX_input.Conversion.Parallel_Modules;
    series = TOOLBOX_input.Conversion.Series_Modules;

    detailed_calc = TOOLBOX_input.Conversion.CableLossesDetailedCalculation;
    if detailed_calc
        if contains(['CEN','POW-OPT'],inv_type)
            str_resistance = TOOLBOX_input.Conversion.StringCableResistivity*...
                TOOLBOX_input.Conversion.StringCableLength/...
                TOOLBOX_input.Conversion.StringCableCrossSection;
        end
        inv_resistance = TOOLBOX_input.Conversion.InverterCableResistivity*...
            TOOLBOX_input.Conversion.InverterCableLength/...
            TOOLBOX_input.Conversion.InverterCableCrossSection;
        fixed_percentage = nan;
    else
        str_resistance = nan; inv_resistance = nan;
        user_cable_ans = nan;
        fixed_percentage = TOOLBOX_input.Conversion.CableLoss;
    end
end

plot_conversion = TOOLBOX_input.Conversion.plot_conversion;


if TOOLBOX_input.runPeriodic && add_central_inv
    nummod = parallel*series;
end

if contains(['STR','CEN','POW-OPT'],inv_type) && add_central_inv
    numstring = parallel;
elseif strcmp(inv_type,'MIC')
    numstring = nummod;
end

%% Build the output structures
simulation.periodic = TOOLBOX_input.runPeriodic;

module.power = ELECTRIC_output.DCP;
module.voltage = ELECTRIC_output.Vmpp;
module.current = ELECTRIC_output.Impp;
module.power_STC = ELECTRIC_output.P_STC;
module.voltage_STC = ELECTRIC_output.Vmpp_STC;
module.current_STC = ELECTRIC_output.Impp_STC;

if strcmp(inv_type,'POW-OPT')
    power_optimizer.model = pow_opt_model;
    power_optimizer.central_inverter = add_central_inv;
    power_optimizer.fixed_voltage = fixed_voltage; %Only for periodic simulations
else
    power_optimizer = struct([]);
end

inverter.type = inv_type;
inverter.plot_graphs = plot_conversion;
if add_central_inv
    inverter.model = inv_model;
    inverter.parallel = parallel;
    inverter.series = series;
    inverter.num_modules = nummod;
    inverter.num_strings = numstring;

    cable_losses.detailed_calc = detailed_calc;
    cable_losses.fixed_percentage = fixed_percentage;
    if contains(['CEN','POW-OPT'],inv_type)
        cable_losses.str_resistance = str_resistance;
    end
    cable_losses.inv_resistance = inv_resistance;
    if ~TOOLBOX_input.script
        cable_losses.user_ans = user_cable_ans;
    end
else
    if ~TOOLBOX_input.runPeriodic
        inverter.num_modules = nummod;
    end
    cable_losses = nan;
end

end

function inv_type = determine_inverter_type(user_choice)
%Determine the type of inverter depending on the user choice
%
% Parameters
% ----------
% user choice : char
%   type of inverter as written in the script version
%
% Returns
% -------
% inv_type : char
%   type of inverter

inv_container = containers.Map;
inv_container('Central') = 'CEN'; inv_container('String') = 'STR';
inv_container('Micro') = 'MIC'; inv_container('Power opt.') = 'POW-OPT';
try
    inv_type = inv_container(user_choice);
catch
    error('conversion:InvalidInverterType', 'Wrong inverter type input')
end
end
