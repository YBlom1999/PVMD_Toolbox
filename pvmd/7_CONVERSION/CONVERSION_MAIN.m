function [CONVERSION_output] = CONVERSION_MAIN(TOOLBOX_input,ELECTRIC_output)
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

%---- Code for power optimizers
if strcmp(inverter.type,'POW-OPT')
    CONVERSION_output = power_optimizer_main(module, inverter, cable_losses, ...
        simulation, power_optimizer,TOOLBOX_input.electric.Mod_Con);
    return;
end

%---- Parameters precalculation for the main code
[system, string] = compute_pv_power(module,inverter,simulation.periodic,TOOLBOX_input.electric.Mod_Con);

Vmax = compute_max_system_voltage(inverter.type,system,string,module);

inv_constants = TOOLBOX_input.Conversion.inv_constants;
if Vmax > TOOLBOX_input.Conversion.inv_constants(13)
    error("conversion:ExceedMaximumVoltage", ['The maximum PV system input voltage of %0.1f V exceeds the ',...
        'maximum voltage of the selected inverter model %s.'],...
        Vmax, inverter.model);
end

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


