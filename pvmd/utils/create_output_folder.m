function folder = create_output_folder(root_folder, project_name, varargin)
%CREATE_OUTPUT_FOLDER
%
% The folder structure used to store the results follows this schema:
%
% - results
%  |- user_name
%    |- project_name / thesis chapter
%      |- datestamp + [description] 
%        |- MODULE_1_config.mat
%        |- MODULE_1_output.mat
%        |- MODULE_2_config.mat
%        |- MODULE_2_output.mat
%        |- TOOLBOX_input.mat
%        \- README.txt
%
%
% Parameters
% ----------
% root_folder : char
%   Absolute or relative path to the output folder, i.e. pointing to top-level
%   /results folder
% project_name : char
%   Name of the project or thesis chapter for students
% user_name : char, optional
%   Name of the researcher. Default is the Matlab username.
% description : char optional
%   Additional identifiers. Default is empty string.
%
%
% Returns
% -------
% folder : char
%   Absolute path of the simulation timestamped output folder
%
%
% Examples
% --------
%
%   folder = create_output_folder('results', 'project_name')
%   folder = create_output_folder('results', 'project_name', ...
%            'user_name', 'John Doe')
%   folder = create_output_folder('results', 'project_name', ...
%            'user_name', 'John Doe', 'description', 'Thermal')
%   

% Set default values
default_user_name = getenv('username');
default_description = '';

% Create parser
p = inputParser;
addRequired(p, 'root_folder', @ischar);
addRequired(p, 'project_name', @ischar);
addParameter(p, 'user_name', default_user_name, @ischar);
addParameter(p, 'description', default_description, @ischar);

% Parse inputs
parse(p, root_folder, project_name, varargin{:});

% Create timestamp with description
formatOut = 'yyyy-mm-dd__HH-MM-SS_';
timestamp = datestr(now, formatOut);
descriptor = strcat(timestamp, p.Results.description);

% Create output folder
folder = fullfile(...
    p.Results.root_folder,...
    p.Results.user_name,...
    p.Results.project_name,...
    descriptor);
mkdir(folder)
end