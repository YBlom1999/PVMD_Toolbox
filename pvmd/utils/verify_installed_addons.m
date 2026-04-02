function out = verify_installed_addons(required)
%VERIFY_INSTALLED_ADDONS 
%
% Parameters
% ----------
% required : cell
%   Cell array with names of required installed addons
%
% Returns
% -------
% logic
%

out = 1;
info = ver;
installed_addons = {info.Name};

for i = 1:length(required)
   if ~any(strcmp(installed_addons, required{i}))
      fprintf('Missing dependency, please install addon %s', required{i})
      out = 0;
   end
end
end