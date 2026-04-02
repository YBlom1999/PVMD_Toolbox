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

if ~TOOLBOX_input.script
    inv_type = ask_user_inverter_type();
    if strcmp(inv_type,'POW-OPT')
        pow_opt_model = ask_user_converter_model('pow-opt');
        add_central_inv = ask_user_yes_no_question(...
            'Would you like to add a central inverter? [Y/N]: ');
        if add_central_inv
            [parallel,series] = ask_user_nr_mod_parallel_series(...
                    TOOLBOX_input.runPeriodic,nummod,inv_type);
            inv_model = ask_user_converter_model('inv');
            fixed_voltage = ask_user_yes_no_question(['Would you like ',...
                'to operate the inverter at a fixed voltage? [Y/N]: ']);
        else
            fixed_voltage = nan;
        end
    else
        [parallel,series] = ask_user_nr_mod_parallel_series(...
                TOOLBOX_input.runPeriodic,nummod,inv_type);
        Vmax = estimate_max_input_voltage(ELECTRIC_output.Vmpp,series,inv_type);
        inv_model = ask_user_inverter_model(Vmax);
        add_central_inv = true;
    end
    
    if add_central_inv
        [detailed_calc,str_resistance,inv_resistance,fixed_percentage,...
            user_cable_ans] = ask_user_cable_losses(inv_type);
    end
    plot_conversion = ask_user_yes_no_question(...
        'Plot graphs of results? [Y/N]: ');
else
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
end

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


function inv_type = ask_user_inverter_type()
%ASK_USER_INVERTER_TYPE Determine type of inverter depending on the user
%
% Returns
% -------
% inv_type : char
%   type of inverter

inverters = ["Central inverter";"String inverter";...
    "Micro inverter";"Power optimiser"];
[idx,~] = listdlg('PromptString','Select Inverter Type',...
    'ListSize',[200 100],...
    'ListString',inverters,'SelectionMode','single',...
    'CancelString','Cancel');
if isempty(idx)
    error('No inverter type has been selected');
end
inverter_types = {'CEN';'STR';'MIC';'POW-OPT'};
inv_type = inverter_types(idx);
inv_type = inv_type{1};
end

function converter_model = ask_user_converter_model(type)
%ASK_USER_CONVERTER_MODEL Ask the user for the converter model
% Valid for both inverter and power optimizer
% 
% Parameters
% ----------
% type : char
%   type of converter: 'inv' or 'pow-opt'
%
% Returns
% -------
% converter_model : char
%   model of converter
cd('Codes')
data_folder = get_data_path(type);
cd('..')
% Get info from data folders
file_list = dir(data_folder);
folders = file_list([file_list(:).isdir]==1);
folders = folders(~ismember({folders(:).name},{'.','..'}));
folders(2).name = folders(2).name + " " + "(NOT RECOMMENDED!)";

% User interaction
Q ='Select Manufacturer';
A = listdlg('PromptString',Q,'SelectionMode','single',...
    'ListString',{folders.name},'ListSize',[150 100]);
if isempty(A); error('Select a manufacturer'); end

if A == 1
    folder_name = 'Solar Edge';
elseif A == 2
    folder_name = 'User Specified';
end

solar_edge_data = fullfile(data_folder,folder_name);
list = dir(fullfile(solar_edge_data, '*.xlsx'));

Q ='Select Power Optimizer Type';
A = listdlg('PromptString',Q,'SelectionMode','single',...
    'ListString',{list.name},'ListSize',[150 200]);
if isempty(A); error('Select a power optimizer'); end

% Type of power optimizer or inverter
converter_model = split(list(A).name,'.xlsx');
converter_model = converter_model{1};
end


function V_max = estimate_max_input_voltage(V_mod,series,inv_type)
%Estimate the maximum system voltage to the inverter
% This method is accurate for periodic simulations but in nonperiodic the
% mismatch between cells need to be considerated
%
% Parameters
% ----------
% V_mod : double/cell
%   DC voltage at module level
% series : double
%   Number of PV modules in series
% inv_type : char
%   type of inverter
% 
% Returns
% -------
% V_max : double
%   Approximate maximum DC voltage in the time period

if strcmp(inv_type,'MIC')
    if iscell(V_mod)
        V_max = max([V_mod{:}],[],'all');
    else
        V_max = max(V_mod);
    end
else
    if iscell(V_mod)
        V_string = cellfun(@(x) x*series, V_mod, 'UniformOutput', false);
        V_max = max([V_string{:}],[],'all');
    else
        V_string = V_mod*series;
        V_max = max(V_string);
    end
end

end


function inv_model = ask_user_inverter_model(Vdc_max)
%ASK_USER_INVERTER_MODEL Ask the user the inverter model
%
% Parameters
% ----------
% Vdc_max : double
%   Approximate maximum DC voltage in the time period
%
% Returns
% -------
% inv_model : char
%   model of inverter


[~,~,data_folder] = get_folder_structure;
file = fullfile(data_folder,'Conversion','constants.xlsx');
[all_constants, all_models] = xlsread(file);
all_models = all_models(4:end,1);

%Remove incompatible inverter types
compatible = all_constants(:,13) > Vdc_max;
models = all_models(compatible);

Q ='Select Inverter Model';
model_idx = listdlg('PromptString',Q,'SelectionMode','single',...
        'ListString',models,...
        'ListSize',[350 500]);
if isempty(model_idx)
    error('No inverter model has been selected');
end
inv_model = models(model_idx);
inv_model = inv_model{1};
end

function [parallel,series] = ask_user_nr_mod_parallel_series(periodic, ...
    nummod, inv_type)
%Ask the user the number of modules in parallel and in series
%
% Parameters
% ----------
% periodic : logical
%   Whether the simulation is periodic
% nummod : double
%   Number of modules in the system (for non-periodic simulations)
% inv_type : char
%   type of inverter
%
% Returns
% -------
% parallel : double
%   Number of modules in parallel
% series : double
%   Number of modules in series

if periodic
    if contains(['CEN','POW-OPT'],inv_type)
        series = input('Specify number of modules in series: ');
        parallel = input('Specify number of modules in parallel: ');
    elseif strcmp(inv_type,'STR')
        series = input('Specify number of modules in the string: ');
        parallel = input('Specify number of strings in parallel: ');
    elseif strcmp(inv_type,'MIC')
        series = 1;
        parallel = input('Specify number of modules in the system: ');
    end
else
    if contains(['STR','CEN','POW-OPT'],inv_type)
        parallel = input("Enter Number of Strings: ");
        series = nummod/parallel;
    elseif strcmp(inv_type,'MIC')
        series = 1;
        parallel = nummod;
    end
end

end

function [detailed_calc,str_R,inv_R,fixed_perc, user_ans] = ask_user_cable_losses(inv_type)
% Ask if the user wants to include cable losses and how
%
% Parameters
% ----------
% inv_type : char
%   type of inverter
%
% Returns
% -------
% detailed_calc : logical
%     Whether the cable losses are calculated with detailed model or not
% str_R : double
%     electric resistance of the string cable
% inv_R : double
%     electric resistance of the inverter cable (if applicable)
% fixed_perc : double
%     fixed percentage of cable losses

info = input("Do you want to include Cable Losses? [Y/N] : ",'s');
if contains('Yy',info)
    choice_cl = {'Fixed Percentage for Cable Losses','Detailed Calculation'};
    [idx,~] = listdlg('PromptString','Select Type of Calculation',...
        'ListSize',[200 100],'ListString',choice_cl,'SelectionMode',...
        'single','CancelString','Cancel');
    if idx == 1
        Q = {'Specify overall cable losses [%]'};
        default = {'1'};
        detailed_calc = false;
        fixed_perc = (str2double(inputdlg(Q,'',1,default)))';
    elseif idx == 2
        detailed_calc = true;
        if contains(['CEN','POW-OPT'],inv_type)
            [str_R,str_ans] = ask_user_cable_properties('string');
        else
            str_R = nan;
            str_ans = struct([]);
        end
        [inv_R,inv_ans] = ask_user_cable_properties('inverter');
        fixed_perc = nan;
        user_ans.inv = inv_ans;
        user_ans.str = str_ans;
    else
        error('Specify the type of cable losses calculation')
    end
elseif contains('Nn',info)
    detailed_calc = false; fixed_perc = 0;
    str_R = nan; inv_R = nan; user_ans = struct([]);
else
    error('Specify the type of cable losses calculation')    
end

if ~detailed_calc; str_R = nan; inv_R = nan; user_ans = struct([]); end

end
    
function [resistance, user_ans] = ask_user_cable_properties(type)
% Obtain the properties of the cable when using the GUI option
%
% Parameters
% ----------
% type : char
%   options 'string' or 'inverter' indicates the cable type
%
% Returns
% -------
% resistance : double
%   Electric resistance of the cable
cd('Codes')
cable_properties_menu(type)
cd ..
try
    material = evalin('base', 'STRING_material');                    
    length = evalin('base', 'STRING_length');
    cross_sec_choice = evalin('base', 'STRING_cross_sec');

    evalin('base','clear STRING_material');
    evalin('base','clear STRING_length');
    evalin('base','clear STRING_cross_sec');
catch
    error('No inputs detected. Please provide proper input.')
end
length = str2double(length);

% Specific resistivity depends on material
% A try/catch condition should be included here!
spec_res_options = [0,0.0168,0.0282];
spec_res = spec_res_options(material);
if spec_res == 0
    error('Select either Copper or Aluminium.')
end

cross_sec_options = [0,0.5,0.75,1,1.5,2,2.5,4,6,10,16,25,35];
cross_sec = cross_sec_options(cross_sec_choice);
if cross_sec == 0
    error('Select a cross section value.')
end

resistance = spec_res*length/cross_sec;
user_ans.spec_res = spec_res;
user_ans.len = length;
user_ans.cross_sec = cross_sec;
end

function answer = ask_user_yes_no_question(ques)
% Ask the user a yes/no question through the command window
%
% Parameters
% ----------
% ques : char
%   Question to be asked
%
% Returns
% -------
% answer : logical
%   User answer

info = input(ques,'s');
answer = contains('Yy',info);
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