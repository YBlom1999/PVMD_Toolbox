function [rgb,a] = read_rgb(list,SET)
% Reads material rgb values and abs status from settings file (for plot).

% ---INPUT---
% list   cell with list of material names (L = nr. of materials on list)
% root   file path
% ---OUTPUT---
% rgb    cell with rgb triplets of each material (values 0 - 1)
% a      absorber status (1=absorber, 0=not)
%
% note: assigns random but reproducible rgb value if material not in file.

L   = length(list);         %nr of materials in device list
rgb = cell(L,1);            %initialize rgb cell
a   = zeros(L,1);           %initialize absorber status vector

ColorTable = readtable(SET.settingsFile,'Sheet','materials','ReadRowNames',true);

%list

for l = 1:L                 % for every layer in device
    if ismissing(list(l)) || isempty(char(list(l)))   % if layer name is missing ('blank' layer)

        rgb{l} = [1,1,1];   % assign white color
        
    elseif ~any(strcmp(list{l},ColorTable.Properties.RowNames)) 
        %if layer name is not in settings file, assign 'random' color

        rng(sum(list{l}))   %string ascii sum as seed for rng
        rgb{l} = rand(1,3); %random (but reproducable) rgb

        %absorber status remains zero (not absorber)

    else                    %if layer name is in settings file

        rgb{l} = ColorTable{list{l},1:3};   %lookup rgb values
        a(l)   = ColorTable{list{l},4};     %lookup absorber status

    end

end
a = logical(a);

end
