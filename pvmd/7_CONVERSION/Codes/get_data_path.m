function data_folder = get_data_path(type_folder)
%Get the relative path to conversion data folders
% Some sections are commented since the idea is to move all data folders
% out of the main code
%
% Parameters
% ----------
% type_folder : char
%   Type of data folder to be accessed: 'pow-opt', 'inv', 'conv'
%
% Returns
% -------
% data_folder : char
%   Relative path to the data folder

% data_folder = strrep(pwd,'pvmd\7_CONVERSION\Codes','data');
[~,~,data_folder] = get_folder_structure;
conversion_data_folder = fullfile(data_folder,'Conversion');
if strcmp(type_folder,'pow-opt')
    data_folder = fullfile(conversion_data_folder,...
        'Power Optimizer Manufacturers');
elseif strcmp(type_folder,'inv')
    data_folder = fullfile(conversion_data_folder,...
        'Inverter Manufacturers');
elseif strcmp(type_folder,'conv')
    data_folder = conversion_data_folder;
end

end

