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
[all_constants, all_models] = xlsread('InverterConstants.xlsx');
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