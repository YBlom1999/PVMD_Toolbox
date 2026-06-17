function Settings = read_settingsBackwardTracer(SettingsFilePath,Settings)
%reads settings from Excel file to struct

set = readcell(SettingsFilePath,"Sheet","Module","Range","A3:B4");    %load settings file

set(cellfun( @(c) isa(c,'missing'),set(:,1)),:) = [];  %remove empty lines

%---convert cell to struct--- (more convenient)
for line = 1:size(set,1)
    Settings.(set{line,1}) = set{line,2};
end
%---
end