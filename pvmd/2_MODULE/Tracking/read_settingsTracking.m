function module_mounting = read_settingsTracking(SettingsFilePath,module_mounting,trackingType)
%reads settings from Excel file to struct


%Read specific settings
if strcmp(trackingType,'HSATS')
    set = readcell(SettingsFilePath,"Sheet","Module","Range","A16:B21");    %load settings file
    set{size(set,1)+1,1} = 'ModAzimuth';
    set{size(set,1),2} = mod(set{6,2}-90,360);
    set{size(set,1)+1,1} = 'ModTilt';
    set{size(set,1),2} = 0;
elseif strcmp(trackingType,'TSATS')
    set = readcell(SettingsFilePath,"Sheet","Module","Range","A25:B31");    %load settings file
    set{size(set,1)+1,1} = 'ModAzimuth';
    set{size(set,1),2} = mod(set{7,2}-90,360);
    set{size(set,1)+1,1} = 'ModTilt';
    set{size(set,1),2} = 0;
elseif strcmp(trackingType,'DATS')
    set = readcell(SettingsFilePath,"Sheet","Module","Range","A35:B40");
    set{size(set,1)+1,1} = 'ModAzimuth';
    set{size(set,1),2} = 0;
    set{size(set,1)+1,1} = 'ModTilt';
    set{size(set,1),2} = 0;
end

set(cellfun( @(c) isa(c,'missing'),set(:,1)),:) = [];  %remove empty lines

%---convert cell to struct--- (more convenient)
for line = 1:size(set,1)
    module_mounting.(set{line,1}) = set{line,2};
end
%---

%Read settings for tracker
set = readcell(SettingsFilePath,"Sheet","Module","Range","A44:B46");
set(cellfun( @(c) isa(c,'missing'),set(:,1)),:) = [];  %remove empty lines

%---convert cell to struct--- (more convenient)
for line = 1:size(set,1)
    module_mounting.(set{line,1}) = set{line,2};
end
%---

end