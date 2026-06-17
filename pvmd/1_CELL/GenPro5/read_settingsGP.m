function SET = read_settingsGP(SettingsFilePath)
%reads settings from Excel file to struct

set = readcell(SettingsFilePath,"Sheet","GENPRO","Range","A2:B30");    %load settings file

set(cellfun( @(c) isa(c,'missing'),set(:,1)),:) = [];  %remove empty lines

%---convert cell to struct--- (more convenient)
SET = struct();                %note: capital letters indicates struct
for line = 1:size(set,1)
    SET.(set{line,1}) = set{line,2};
end
SET.settingsFile= SettingsFilePath;
%---
end