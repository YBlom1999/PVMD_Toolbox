%new function for the script version of the PVMD TOOLBOX
%Co-developed by Karthik Ganapathi Subramanian & Malte Vogt
%Source code courtesy of Malte Vogt & Rudi Santbergen, PVMD Group, TU Delft
% Delft, the Netherlands.

% ©All rights reserved.

function [TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_output,DEGRADATION_output,LCOE_output,LOSS_ANALYSIS_output] = TB_script(input)

% Add functions to MATLAB path
here = pwd;
addpath(genpath(here));

%---- Check installed addons
required = {'Parallel Computing Toolbox', 'Symbolic Math Toolbox'};
if ~verify_installed_addons(required), return; end %if addons are missing, don't continue

%copy input to workspace
load(input,'TOOLBOX_input');

% Set up save folder
if TOOLBOX_input.save.enable
    TOOLBOX_input.save.folder = create_output_folder(...
        TOOLBOX_input.save.root_path,...
        TOOLBOX_input.save.project_name,...
        'user_name', TOOLBOX_input.save.user_name,...
        'description', TOOLBOX_input.save.description);
else
    TOOLBOX_input.save.folder = '';
end

%this needs to be true for the script version to run
TOOLBOX_input.script=true;

%% 1_CELL
if TOOLBOX_input.runDeviceOptic %running CELL simulations
    [CELL_output,TOOLBOX_input] = CELL_main(TOOLBOX_input);           %run the main file there
    disp('CELL Computations Completed!');
end

%% 2_MODULE
if TOOLBOX_input.runModulePart %running MODULE Simulations
    [MODULE_output, TOOLBOX_input] = MODULE_main(TOOLBOX_input, CELL_output);
    disp('MODULE Computations Completed!');
else
    MODULE_output.Output = nan;

end

%% 3_WEATHER
if TOOLBOX_input.runWeatherPart
        % Run weather module
        [WEATHER_output, TOOLBOX_input] = WEATHER_main(TOOLBOX_input, MODULE_output);
else
    WEATHER_output.Output = nan;
end
disp('WEATHER Computations Completed!');

%% 4_THERMAL
if TOOLBOX_input.runThermalPart
    
    % Run thermal module
    [THERMAL_output, TOOLBOX_input] = THERMAL_main(TOOLBOX_input, CELL_output, MODULE_output, WEATHER_output);
    

else
    THERMAL_output.Output = nan;

end

%% 5_ELECTRIC
if TOOLBOX_input.runElectricalPart
    [ELECTRIC_output, TOOLBOX_input] = ELECTRIC_main(TOOLBOX_input, CELL_output, MODULE_output, ...
                WEATHER_output,THERMAL_output); %running script for ELECTRIC calculations
    disp('ELECTRIC Computations Completed!');
else
    ELECTRIC_output.Output = nan;

end

%% 6_DEGRADATION
%running degradation calculations
if TOOLBOX_input.runDegradationPart   %for yearly simulations
    %calculating LCOE
    [DEGRADATION_output, TOOLBOX_input] = DEGRADATION_main(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output);
else %no LCOE calculation
    DEGRADATION_output.Output = nan;
%     assignin('base','DEGRADATION_output',DEGRADATION_output);
end

%% 7_CONVERSION
%running AC Conversion
if TOOLBOX_input.runACConversionPart
    % Run conversion module
    [CONVERSION_output,TOOLBOX_input] = CONVERSION_MAIN(TOOLBOX_input,ELECTRIC_output);

else
    CONVERSION_output.Output = nan;

end

%% 8_FINANCIALS
%running financial calculations
if TOOLBOX_input.runFinancialPart
    %calculating LCOE
    [LCOE_output, TOOLBOX_input] = LCOE_MAIN(TOOLBOX_input,ELECTRIC_output,DEGRADATION_output);
else %no LCOE calculation
    LCOE_output.Output = nan;
end

%% 9_LOSS ANALYSIS
%running loss analysis calculations
if TOOLBOX_input.runLA
    %calculating LCOE
    [LOSS_ANALYSIS_output, TOOLBOX_input] = LOSS_ANALYSIS_main(TOOLBOX_input,CELL_output,MODULE_output,WEATHER_output,THERMAL_output,ELECTRIC_output,CONVERSION_output);
else %no Loss analysis calculation
    LOSS_ANALYSIS_output.Output = nan;
end

%% Save simulation parameters
if TOOLBOX_input.save.enable
    disp('Saving Results...')
    file = fullfile(TOOLBOX_input.save.folder, 'Results.mat');
    save(file, 'TOOLBOX_input',"CELL_output","MODULE_output","WEATHER_output","THERMAL_output","ELECTRIC_output","DEGRADATION_output","CONVERSION_output","LCOE_output","LOSS_ANALYSIS_output");
end

end