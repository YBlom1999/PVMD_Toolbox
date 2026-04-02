function [src_folder,test_folder,data_folder] = get_folder_structure
%GET_FOLDER_STRUCTURE Get src and tests folders
%   
% This file needs to reside in `/pvmd/utils` to produce the correct paths
%
% Returns
% -------
% src_folder : char
%   path to the pvmd folder
% test_folder : char
%   path to the tests folder
% data_folder : char
%   path to the data folder

% Define folder structure
here = mfilename('fullpath');
idx = strfind(here, filesep);
root = here(1:idx(end-2)-1);

src_folder = fullfile(root, 'pvmd');
test_folder = fullfile(root, 'tests');
data_folder = fullfile(root, 'data');
end

