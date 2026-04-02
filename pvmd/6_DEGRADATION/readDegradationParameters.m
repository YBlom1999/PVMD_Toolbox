function Parameters = readDegradationParameters(TOOLBOX_input)
% readDegradationParamers reads the values of the different parameters
% involved in the degradation simulation
%
% It reads the filename with the requisted parameters and writes them to
% the structure parameters
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   User input parameters to the toolbox
%
% Returns
% -------
% Parameters : struct
%   Update of TOOLBOX_input by adding simulation parameters of degradation module
%
% Author: Youri Blom

filename = TOOLBOX_input.Degradation.Parameters_file;
[~,~,data_folder] = get_folder_structure;
files_dir = fullfile(data_folder,'Degradation','Parameters',filename);
load(files_dir,'A_mois','C_mois','Ea_mois','A_dis','Ea_dis','factor_LID','C_LID','Ea_TC','A_TC','n_TC','b_TC');

Parameters.A_mois = A_mois;
Parameters.C_mois = C_mois;
Parameters.Ea_mois = Ea_mois;
Parameters.A_dis = A_dis;
Parameters.Ea_dis = Ea_dis;
Parameters.factor_max = factor_LID;
Parameters.C_LID = C_LID;
Parameters.A_TC = A_TC;
Parameters.Ea_TC = Ea_TC;
Parameters.n_TC = n_TC;
Parameters.b_TC = b_TC;

end