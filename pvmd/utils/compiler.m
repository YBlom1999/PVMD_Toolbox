function compiler
% P-code compiler
%
% Create exact copy of matlab source code and replace all matlab files with
% pcode

% Get directory structure
here = pwd;
root = fullfile(here, '../../');
src_matlab = fullfile(root, 'pvmd');

% Create output fodler for pcode
output_folder = fullfile(root, 'pcode');
if ~isfolder(output_folder)
    mkdir(output_folder)
end

% Copy Matlab code base
copyfile(src_matlab, output_folder)

% Replace .m files with .p files
list = dir(fullfile(output_folder, '**/*.m'));
for i = 1:length(list)
    file = fullfile(list(i).folder, list(i).name);
    pcode(file,'-inplace');
    delete(file);
end

end
