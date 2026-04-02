function Pdc = include_cable_losses(Pdc_raw,Idc_string,Idc_system,...
    cable_losses,inv_type)
%Modify the DC power to include cable losses
% These can be computed either with a detailed calculation considering the
% cable characteristics, or with a fixed percentage
%
% Parameters
% ----------
% Pdc_raw : double/cell
%   DC power at system, string (string inverter) or module (micro inverter)
%   level
% Idc_system : double/cell
%   DC current at system level
% Idc_string : double/cell
%   DC current at string level. For some inverters and simulation types,
%   Idc_string = Idc_system
% cable_losses : struct
%   User inputs related to the cable losses
% inv_type : char
%   Type of inverter considered in the simulations
%
% Returns
% -------
% Pdc : double
%   DC power after cable losses have been considered

if cable_losses.detailed_calc
    Pdc = cable_losses_detailed_calculation(Pdc_raw,Idc_string,Idc_system,...
        cable_losses,inv_type);
else
    Pdc = cable_losses_fixed_percentage(Pdc_raw,...
        cable_losses.fixed_percentage);
end

end

function Pdc = cable_losses_fixed_percentage(Pdc_raw,cable_loss)
%Compute the cable losses assuming a fixed percentage
%
% Parameters
% ----------
% Pdc_raw : double
%   DC power at system, string (string inverter) or module (micro inverter)
%   level
% cable_loss : double
%   Percentage of cable losses between 0 an 100 (usually 1)
%
% Returns
% -------
% Pdc : double
%   DC power after considering cable losses

cable_loss = cable_loss/100;
if iscell(Pdc_raw)
    Pdc = cellfun(@(x) x*(1-cable_loss),Pdc_raw,'UniformOutput',false);
else
    Pdc = Pdc_raw*(1-cable_loss);
end
end

function Pdc = cable_losses_detailed_calculation(Pdc_raw,Idc_string,...
    Idc_system,cable_losses,inv_type)
%Calculate the cable losses considering the cable characteristics
% Start by computing the total cable losses considering the resistances
% from the module to the inverter. Compute the percentual cable loss to
% assert a warning if these exceed 2%. Finally, modify the input DC power
% by substracting the losses due to cables
%
% Parameters
% ----------
% Pdc_raw : double
%   DC power at system, string (string inverter) or module (micro inverter)
%   level depending on the inverter and simulation type
% Idc_string : double/cell
%   DC current at string level. For some inverters and simulation types,
%   Idc_string = Idc_system
% Idc_system : double/cell
%   DC current at system level
% cable_losses : struct
%   User choices for the computation of cable losses
% inv_type : char
%   Type of inverter considered in the simulations
%
% Returns
% -------
% Pdc : double
%   DC power after cable losses have been considered

if contains(['MIC','STR'],inv_type)
    total_cable_loss = get_cable_loss(Idc_string, cable_losses.inv_resistance);
elseif contains(['POW-OPT','CEN'],inv_type)
    string_cable_loss = get_cable_loss(Idc_string, cable_losses.str_resistance);
    if iscell(string_cable_loss)
        string_cable_loss = sum([string_cable_loss{:}],2);
    end
    inv_cable_loss = get_cable_loss(Idc_system, cable_losses.inv_resistance);
    total_cable_loss = inv_cable_loss + string_cable_loss;
end

if iscell(Pdc_raw)
    percentual_power_loss = cellfun(@(x,y) x./y*100,total_cable_loss,Pdc_raw,...
        'UniformOutput',false);
    percentual_power_loss = cell2mat(percentual_power_loss);
else
    percentual_power_loss = total_cable_loss./Pdc_raw*100;
end
percentual_power_loss(isnan(percentual_power_loss)) = 0;
if max(percentual_power_loss,[],'all') > 2
    warning('conversion:HighPowerLoss', 'Cable losses are above 2%%. Consider a larger diameter cable.');
end

if iscell(Pdc_raw)
    Pdc = cellfun(@(x,y) x-y,Pdc_raw,total_cable_loss,...
        'UniformOutput',false);
else
    Pdc = Pdc_raw - total_cable_loss;
end

end

function cable_loss = get_cable_loss(Idc, resistance)
%Compute the cable losses considering current and resistance
%
% Parameters
% ----------
% Idc : cell/double
%   Electric current
% resistance : double
%   Electric resistance
%
% Returns
% -------
% cable_loss : cell/double
%   Losses in the cable
if iscell(Idc)
    cable_loss = cellfun(@(x) x.^2*resistance,Idc,'UniformOutput',false);
else
    cable_loss = Idc.^2*resistance;
end
end
