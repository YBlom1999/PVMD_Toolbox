function TOOLBOX_input = modify_toolbox_input(TOOLBOX_input, inverter, ...
    power_optimizer, cable_losses)
%Modify the struct toolbox input by adding the user simulation choices
% Used only if the simulations are run with the GUI
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
% inverter : struct
%   User inputs related to the inverter characteristics
% power_optimizer : struct
%   User inputs related to the power optimizer characteristics
% cable_losses : struct
%   User inputs related to the cable losses
% 
% Returns
% -------
% TOOLBOX_input : struct
%   User input parameters to the toolbox

TOOLBOX_input.runACConversionPart = true;
TOOLBOX_input.Conversion.plot_conversion = inverter.plot_graphs;
TOOLBOX_input.Conversion.ConversionType = inverter.type;


if ~isempty(power_optimizer)
    TOOLBOX_input.Conversion.PowerOpt.Model = power_optimizer.model;
    TOOLBOX_input.Conversion.PowerOpt.AddCentralInverter = power_optimizer.central_inverter;
    TOOLBOX_input.Conversion.PowerOpt.FixedVoltage = power_optimizer.fixed_voltage;
else
    power_optimizer = struct();
    power_optimizer.central_inverter = true;
end

if power_optimizer.central_inverter
    TOOLBOX_input.Conversion.Model = inverter.model;
    TOOLBOX_input.Conversion.Parallel_Modules = inverter.parallel;
    TOOLBOX_input.Conversion.Series_Modules = inverter.series;
    
    TOOLBOX_input.Conversion.CableLossesDetailedCalculation = cable_losses.detailed_calc;
    TOOLBOX_input.Conversion.CableLoss = cable_losses.fixed_percentage;
    if ~isempty(cable_losses.user_ans)
        TOOLBOX_input.Conversion.InverterCableResistivity = ...
            cable_losses.user_ans.inv.spec_res;
        TOOLBOX_input.Conversion.InverterCableLength = ...
            cable_losses.user_ans.inv.len;
        TOOLBOX_input.Conversion.InverterCableCrossSection = ...
            cable_losses.user_ans.inv.cross_sec;
        if ~isempty(cable_losses.user_ans.str)
            TOOLBOX_input.Conversion.StringCableResistivity = ...
                cable_losses.user_ans.str.spec_res;
            TOOLBOX_input.Conversion.StringCableLength = ...
                cable_losses.user_ans.str.len;
            TOOLBOX_input.Conversion.StringCableCrossSection = ...
                cable_losse