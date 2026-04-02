function config = get_user_input(prompt, field_names, varargin)
%GET_USER_INPUT Get user input through dialog box
%
%
% Parameters
% ----------
% prompt : cell
%   Questions for inpu, e.g. prompt={'What is your name?', 'How old are you?'}
% field_names : cell
%   Cell with the names for the fields of the returned structure. 
%   Each field name corresponds to a question. Must be of equal length to prompt.
% default_input : cell, optional
%   Optional default values, each value must be of type char, 
%   i.e. {'1', '2', '3'}. Default is empty cell.
% title : char, optional
%   Optional dialog box title. Default is '' 
%
%
% Returns
% -------
% struct
%   config parameters of type double for each field name
%
%
% Examples
% --------
%   config = get_user_input({'Input 1?', 'Input 2?'}, {'in_1', 'in_2'})
%   config = get_user_input({'Input 1?', 'Input 2?'}, {'in_1', 'in_2'},...
%            {'10', '20'})
%   config = get_user_input({'Input 1?', 'Input 2?'}, {'in_1', 'in_2'},...
%            'title', 'Input parameters')
%   config = get_user_input({'Input 1?', 'Input 2?'}, {'in_1', 'in_2'},...
%            {'10', '20'}, 'title', 'Input parameters')

% Set defaults
default_title = '';
default_input = repmat({''}, length(prompt), 1); % cell with empty strings of length prompt

% Parse arguments
p = inputParser;
addRequired(p, 'prompt', @iscell);
addRequired(p, 'field_names', @iscell);
addOptional(p, 'default_input', default_input, @iscell);
addParameter(p, 'title', default_title, @ischar);
parse(p, prompt, field_names, varargin{:});

% Execute input dialog
opts.Resize = 'off';
answers = inputdlg(p.Results.prompt, p.Results.title, 1,...
    p.Results.default_input, opts);

% Convert answers to double
answers_double = cellfun(@str2double, answers);

% Store input in config structure
config = struct([]);
if ~isempty(answers_double)    
    for i = 1 : length(field_names)
       config(1).(field_names{i}) = answers_double(i); 
    end
end
end