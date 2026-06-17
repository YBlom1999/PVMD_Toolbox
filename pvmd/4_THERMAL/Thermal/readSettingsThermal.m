function thermalSettings = readSettingsThermal(SettingsFilePath,thermalSettings)
%reads settings from Excel file to struct


%Read specific settings
set = readcell(SettingsFilePath,"Sheet","Thermal","Range","A3:B5");    %load settings file


set(cellfun( @(c) isa(c,'missing'),set(:,1)),:) = [];  %remove empty lines

%---convert cell to struct--- (more convenient)
for line = 1:size(set,1)
    thermalSettings.(set{line,1}) = set{line,2};
end
%---


end