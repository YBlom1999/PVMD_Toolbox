function TOOLBOX_input = getDegradationInputs(TOOLBOX_input)
% getDegradationInputs Prepares the inputs to the main degradation function
%
% It asks the user questions if the script version is not run. These
% questions are required inputs for the degradation model
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
%
% Returns
% -------
% TOOLBOX_input : struct
%   Update of TOOLBOX_input by adding simulation parameters of degradation module
%
% Author: Youri Blom

%Ask for number of working years
N_years = input('Number of working years: ');
TOOLBOX_input.Degradation.N_years = N_years;

%Ask for RMC model
Options_RMC_model = ["Load COMSOL result";"Use Simple Equation"];
[idx,~] = listdlg('PromptString','Select Inverter Type',...
    'ListSize',[200 100],...
    'ListString',Options_RMC_model,'SelectionMode','single',...
    'CancelString','Cancel');

if idx == 1
    TOOLBOX_input.Degradation.loadCOMSOL = true;
    TOOLBOX_input.Degradation.simpleEquation = false;
    TOOLBOX_input.Degradation.COMSOL_file = askFile('RMC');
    TOOLBOX_input.Degradation.RMC_sat = nan;
    TOOLBOX_input.Degradation.T_95 = nan;
elseif idx == 2
    TOOLBOX_input.Degradation.loadCOMSOL = false;
    TOOLBOX_input.Degradation.simpleEquation = true;
    TOOLBOX_input.Degradation.COMSOL_file = ''; 

    %Ask for parameters
    [RMC_sat,T_95] = askParametersSimpleEquation();
    TOOLBOX_input.Degradation.RMC_sat = RMC_sat;
    TOOLBOX_input.Degradation.T_95 = T_95;

end

%Ask for Moisture ingress model
TOOLBOX_input.Degradation.MoistureDeg_file = askFile('MoistureIngress');

%Ask for Discoloration model
TOOLBOX_input.Degradation.Discoloration_file = askFile('Discoloration');
end


function sel_file = askFile(folder_name)
% askFile Asks the user which file it wants to use
%
% It shows a list of all possible simulations. It will return the filename
% of the selected simulation
%
% Parameters
% ----------
%
% Returns
% -------
% sel_file : string
%   Name of the selected result
%
% Author: Youri Blom
[~,~,data_folder] = get_folder_structure;
files_dir = fullfile(data_folder,'Degradation',folder_name);
prompt = 'Select the file';
list = dir([files_dir,filesep,'*.mat']);
A = listdlg('PromptString',prompt,'SelectionMode','single',...
    'ListString',{list.name},'ListSize',[200,100]);
if isempty(A), return, end

sel_file = list(A).name;

end

function [RMC_sat,T_95] = askParametersSimpleEquation()
% askParametersSimpleEquation Asks the user for the parameters of the
% simple equation model
%
% It asks for the saturation value and the time at which 95% of this value
% is reached.
%
% Parameters
% ----------
%
% Returns
% -------
% COMSOL_file : string
%   Name of the selected Comsol result
%
% Author: Youri Blom
prompt1 = {'Saturation value of RMC',...
    'Time at which the 95% saturation is reached [years]'};
dlgtitle = 'Data for simple RMC model';
dims = [1 40];
definput = {'0.5','8'};
Answers = inputdlg(prompt1,dlgtitle,dims,definput);

RMC_sat = str2double(Answers{1});
T_95 = str2double(Answers{2});

end