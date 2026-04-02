function TOOLBOX_input = set_config_thermal(TOOLBOX_input)
%SET_CONFIG_THERMAL Get user input if parameters are missing
%
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
%
% Returns
% -------
% TOOLBOX_input : struct
%


% Generate input parameters
if ~TOOLBOX_input.script
    TOOLBOX_input.runThermalPart = true; % FIXME: remove writing to workspace        
    THERMAL_type = ask_thermal_type;

    if THERMAL_type == 1
    
        THERMAL_config = get_user_input_thermal;
    
        % Abort if user presses 'cancel'
        if isempty(THERMAL_config), return, end     
        
        TOOLBOX_input.thermal = THERMAL_config; 
        TOOLBOX_input.thermal.Type = 1;   
        
        % Enable plotting of the thermal figures
        TOOLBOX_input.thermal.plot_thermal = true;
    else
        TOOLBOX_input.thermal.Type = 2;
        TOOLBOX_input.thermal.runpvt = get_user_input_PVT;
    end
   
end
end

function type = ask_thermal_type
%ask_thermal_type Generate dialog for user input on system type
%
% Returns
% -------
% type : double
%   Indication of system type (1 is PV module, 2 is PV Thermal module)

usrops = {'Normal solar cell','PV Thermal module'};
promptstring = 'Select Type';
type = listdlg('PromptString',promptstring...
    ,'ListString',usrops,'SelectionMode','single','ListSize',[200 100]);
end


function config = get_user_input_thermal
%get_user_input_thermal Generate dialog for user input on thermals
%
% Returns
% -------
% config : struct
%   Input parameters for thermal
Q ={'Select temperature model'};    %ask user
TempModel = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Fluid dynamic model','Duffie-Beckman model','Faiman model','Sandia model','Incropera model'},'ListSize',[150,100]); %user choice
if TempModel == 1
    Temperature_model = 'Fluid dynamic model';
    prompt = {
        'Efficiency of the PV module(Only used for temperature calculation)',...
        'Temperature coefficient of the PV module [1/K]',...
        'Thickness of the glass [m]',...
        };
    alias = {'cell_eff', 'temp_coeff', 'glass_thickness'};
    dlgtitle = 'Thermal module parameters';
    default_values = {'0.25','-0.0035','0.003'};

    config = inputdlg(prompt, dlgtitle,[1 35],default_values);
    config = struct(alias{1},str2double(config{1}),alias{2},str2double(config{2}),alias{3},str2double(config{3}));


elseif TempModel == 2
    Temperature_model = 'Duffie-Beckman model';
    prompt = {
        'Efficiency of the PV module(Only used for temperature calculation)',...
        'Nominal operating cell temperature (T_NOCT)',...
        };
    alias = {'cell_eff', 'T_NOCT'};
    dlgtitle = 'Thermal module parameters';
    default_values = {'0.25','60'};

    config = inputdlg(prompt, dlgtitle,[1 35],default_values);
    config = struct(alias{1},str2double(config{1}),alias{2},str2double(config{2}));
elseif TempModel == 3
    Temperature_model = 'Faiman model';
    prompt = {
        'Efficiency of the PV module(Only used for temperature calculation)',...
        'Constant heat transfer U0',...
        'Convective heat transfer U1', ...
        'Alpha'...
        };
    alias = {'cell_eff', 'U0','U1','alpha'};
    dlgtitle = 'Thermal module parameters';
    default_values = {'0.25','0','29','0.9'};

    config = inputdlg(prompt, dlgtitle,[1 35],default_values);
    config = struct(alias{1},str2double(config{1}),alias{2},str2double(config{2}),alias{3},str2double(config{3}),alias{4},str2double(config{4}));
elseif TempModel == 4
    Temperature_model = 'Sandia model';
    prompt = {
        'a',...
        'b',...
        };
    alias = {'a', 'b'};
    dlgtitle = 'Thermal module parameters';
    default_values = {'-3.47','-0.0594'};

    config = inputdlg(prompt, dlgtitle,[1 35],default_values);
    config = struct(alias{1},str2double(config{1}),alias{2},str2double(config{2}));
elseif TempModel == 5
    Temperature_model = 'Incropera model';
    prompt = {
        'Number of layers',...
        'Rear Convection',...
        'Rear temperature',...
        'Efficiency of the PV module(Only used for temperature calculation) [%]',...
        };
    alias = {'Nlayers','RearConvection', 'RearTemperature','Efficiency'};
    dlgtitle = 'Thermal module parameters';
    default_values = {'5','1','nan','20'};

    config = inputdlg(prompt, dlgtitle,[1 35],default_values);
    config = struct(alias{1},str2double(config{1}),alias{2},str2double(config{2}),alias{3},str2double(config{3}),alias{4},str2double(config{4}));
    layers = get_user_input_layers(config.Nlayers);
    assignment_layers = assign_layers_users(CELL_output.CELL_FRONT.lay,layers);
    config.layers = layers;
    config.assignment_layers = assignment_layers;
    
end
config.Temperature_model = Temperature_model;
end

function [layers] = get_user_input_layers(Nlayers)
%assign_layers_users aks the user for the parameters of all layers
%
% Returns
% -------
% layers : element
%   parameters of all layers
layers = repelem(PvLayer('name','','thickness',0,'k',0,'rho',0,'cp',0,'dx',0,'nIntNodeY',0,'emissivity',0),1,Nlayers);
for i = 1:Nlayers
    prompt = {
        'Name layer',...
        'thickness [m]',...
        'k (thermal conductivity) [W/m K]', ...
        'rho (density) [kg/m3]'...
        'cp (heat capacity) [J/kg K]'...
        'dx [m]'...
        'Number of nodes Y'...
        'emissivity [-]'...
        };
    dlgtitle = 'Thermal module parameters';
    default_values = {'Encapsulation','1e-3','0.24','1700','250','0.25','6','0.89'};
    answer =  inputdlg(prompt, dlgtitle,[1 35],default_values);
    layers(i) = PvLayer('name',answer{1},'thickness',str2double(answer{2}),'k',str2double(answer{3}),'rho',str2double(answer{4}),'cp',str2double(answer{5}),'dx',str2double(answer{6}),'nIntNodeY',str2double(answer{7}),'emissivity',str2double(answer{8}));
end
end

function [assignment_layers] = assign_layers_users(layers_GENPRO,layers_Incropera)
%assign_layers_users assigns the layer to the correct cel part
%
% Returns
% -------
% assignment_layers : double
%   overview of assignment of the layers
choice_layers = cell(1,length(layers_Incropera));
for layer_Incropera_i = 1:length(layers_Incropera)
    choice_layers(layer_Incropera_i) = {layers_Incropera(layer_Incropera_i).name};
end
assignment_layers = zeros(1,length(layers_GENPRO));
for layer_GENPRO_i = 2:length(layers_GENPRO)-1
    Q =append('Select layer in Incropera for ',layers_GENPRO(layer_GENPRO_i));    %ask user
    assignment_layers(layer_GENPRO_i) = listdlg('PromptString',Q,'SelectionMode','single','ListString',choice_layers,'ListSize',[250,100]); %user choice

end

end

function runpvt = get_user_input_PVT
%get_user_input_PVT Generate dialog for user input on PVT
%
% Returns
% -------
% runpvt : struct
%   function that should be executed
    Q_Q ={'Select option:'};    %ask user
    choice = listdlg('PromptString',Q_Q,'SelectionMode','single','ListString',{'Run PV-thermal Collector files', 'Run Thermal Collector files', 'Run Storage Tank files', 'Run Heat Pump files'},'ListSize',[150,80]); %user choice
    if isempty(choice), return, end
    if choice==1
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        types_folder = fullfile(pvmd_folder, '4_THERMAL','PV_Thermal','PV_Thermal_Collector');
        cd(types_folder);                               %go to the "PVthermal" folder
        list = dir('*.m');                              %list the file names there
        Q_Q ='Run saved simulations';                   %ask user
        A = listdlg('PromptString',Q_Q,'SelectionMode','single','ListString',{list.name}); %user choice
        cd(current_path);                               %return to the original folder

        if isempty(A), return, end                      %if user presses 'cancel'

    elseif choice == 2
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        types_folder = fullfile(pvmd_folder, '4_THERMAL','PV_Thermal','Thermal_Collector');
        cd(types_folder);                               %go to the "Thermal" folder
        list = dir('*.m');                              %list the file names there
        Q_Q ='Define thermal collector type';           %ask user
        A = listdlg('PromptString',Q_Q,'SelectionMode','single','ListString',{list.name}); %user choice
        cd(current_path);                               %return to the original folder

        if isempty(A), return, end                      %if user presses 'cancel'
    elseif choice == 3
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        types_folder = fullfile(pvmd_folder, '4_THERMAL','PV_Thermal','Storage_Tank');
        cd(types_folder);                               %go to the "Tank" folder
        list = dir('*.m');                              %list the file names there
        Q_Q ='Define collector with tank type';         %ask user
        A = listdlg('PromptString',Q_Q,'SelectionMode','single','ListString',{list.name}); %user choice
        cd(current_path);                               %return to the original folder

        if isempty(A), return, end                      %if user presses 'cancel'
    elseif choice == 4
        current_path = pwd;
        [pvmd_folder,~,~] = get_folder_structure;
        types_folder = fullfile(pvmd_folder, '4_THERMAL','PV_Thermal','Heat_Pump');
        cd(types_folder);                               %go to the "Tank" folder
        list = dir('*.m');                              %list the file names there
        Q_Q ='Define collector with HP type';           %ask user
        A = listdlg('PromptString',Q_Q,'SelectionMode','single','ListString',{list.name}); %user choice
        cd(current_path);                               %return to the original folder

        if isempty(A), return, end                      %if user presses 'cancel'
    end
    %copy choices to workspace
    runpvt=true;
    runpvt=list(A).name;
    cd(current_path)                                %return to parent folder
end