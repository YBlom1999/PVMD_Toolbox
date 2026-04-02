function THERMAL_output = calculate_cell_temp(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output)
%CALCULATE_CELL_TEMP Calculate the cell temperature
%
% This function calculates the temperature of the cells in the system
% using the fluid dynamic model, therefore considering the weather 
% and the absorbed irradiance.
% 
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% WEATHER : struct
%   Results from the WEATHER module
% MODULE : struct
%   Results from the MODULE module
%
% Returns
% -------
% THERMAL_output : struct
%   Simulation results of the thermal module
% 
% Developed by unknown (A. Jamodkar? E. Garcia?). Improved by A. Nour.
% Commented by A. Alcaniz
% Refactored by M. Kok

%---- Create output structures
THERMAL_output = struct([]);

%---- Verify presence of temperature model
fields = {'Temperature_model'};
missing_fields = ~isfield(TOOLBOX_input.thermal, fields);
if any(missing_fields)
    missing_parameters = fields(missing_fields);
    fprintf('Aborting thermal calculation, missing parameters:\n')
    fprintf('    - %s\n', missing_parameters{:})
    return
end

%---- Verify presence of required input parameters
if strcmp(TOOLBOX_input.thermal.Temperature_model,'Fluid dynamic model')
    fields = {'cell_eff', 'temp_coeff', 'glass_thickness', 'plot_thermal'};
elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Duffie-Beckman model')
    fields = {'cell_eff', 'T_NOCT'};
elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Faiman model')
    fields = {'cell_eff', 'U0','U1','alpha'};
elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Sandia model')
    fields = {'a', 'b'};
elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Incropera model')
    fields = {'Nlayers','RearConvection', 'RearTemperature','Efficiency'};
    if ~check_input_Incropera(TOOLBOX_input.thermal)
        fprintf('Aborting thermal calculation, Input invalid\n')
        return
    end

end
missing_fields = ~isfield(TOOLBOX_input.thermal, fields);
if any(missing_fields)
    missing_parameters = fields(missing_fields);
    fprintf('Aborting thermal calculation, missing parameters:\n')
    fprintf('    - %s\n', missing_parameters{:})
    return
end

%---- Calculate cell temperatures
disp('Calculating Cell Temperatures. This may take a minute or two...');

%---- Find number of modules
if TOOLBOX_input.runPeriodic
    num_modules = 1;  
else
    num_modules = TOOLBOX_input.Scene.N_panels;
end

%---- Calculate plane of array irradiance
poa_cell = WEATHER_output.A;

%---- Calculate cell temperature using the selected thermal model
cell_temperature = cell(numel(poa_cell),1);
for i = 1:numel(poa_cell)
    if strcmp(TOOLBOX_input.thermal.Temperature_model,'Fluid dynamic model')
        if TOOLBOX_input.deviceOptic.exportAbsorptanceAll
            poa_cell{i} = poa_cell{i}(:,:,1);
        end
        cell_temperature{i} = individual_cell_FD(...
            WEATHER_output.ambient_temperature,...
            WEATHER_output.wind_speed,...
            poa_cell{i},...
            TOOLBOX_input,...
            MODULE_output,i);
    elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Duffie-Beckman model')
        cell_temperature{i} = cellTempCalcDB(WEATHER_output.ambient_temperature, ...
            WEATHER_output.wind_speed,poa_cell{i},TOOLBOX_input.thermal);

    elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Faiman model')
        cell_temperature{i} = cellTempCalcFaiman(WEATHER_output.ambient_temperature, ...
            WEATHER_output.wind_speed,poa_cell{i},TOOLBOX_input.thermal);

    elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Sandia model')
        cell_temperature{i} = cellTempCalcSandia(WEATHER_output.ambient_temperature, ...
            WEATHER_output.wind_speed,poa_cell{i},TOOLBOX_input.thermal);
    
    elseif strcmp(TOOLBOX_input.thermal.Temperature_model,'Incropera model')
        Temperature_Cells = zeros(length(WEATHER_output.day),MODULE_output.N);
        % Find index of solar cell
        layer_cell = TOOLBOX_input.thermal.assignment_layers(CELL_output.CELL_FRONT.Absmat_ind(2));
        ind_cell_Y = 0;
        for layer_i = 1:layer_cell-1
            ind_cell_Y = ind_cell_Y+TOOLBOX_input.thermal.layers(layer_i).nIntNodeY;
        end
        ind_cell_Y = ind_cell_Y + round(TOOLBOX_input.thermal.layers(layer_cell).nIntNodeY/2);


        for cell_i = 1:MODULE_output.N
            [pvModule,weather] = prepare_incropera(TOOLBOX_input, WEATHER_output,CELL_output,cell_i);
            [moduleTemp,~,~] = cellTemp_incropera2D(pvModule,weather,'rear convection',TOOLBOX_input.thermal.RearConvection,'rear temperature',TOOLBOX_input.thermal.RearTemperature);

            ind_cell_X = round((1 + round(pvModule.width/pvModule.layers(1).dx))/2);
            Temperature_Cells(:,cell_i) = moduleTemp(ind_cell_Y,ind_cell_X,:);
        end
        cell_temperature{i} = Temperature_Cells;
    end
end

% Add temperature to output structure
THERMAL_output(1).T = cell_temperature;

%---- Optional: plot figures
if TOOLBOX_input.thermal.plot_thermal

    if TOOLBOX_input.runPeriodic 
        figures_cell_temperature_per(...
            WEATHER_output.Irr,...
            THERMAL_output.T,...
            WEATHER_output.ambient_temperature,...
            WEATHER_output.day,...
            WEATHER_output.month, ...
            WEATHER_output.Period);
    elseif ~TOOLBOX_input.runPeriodic
        figures_cell_temperature_nonper(...
            WEATHER_output.Irr,...
            THERMAL_output.T,...
            WEATHER_output.ambient_temperature,...
            WEATHER_output,...
            TOOLBOX_input);
    end
end

disp('Thermal calculation finished.')
end
