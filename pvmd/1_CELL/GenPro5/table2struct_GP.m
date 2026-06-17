function [Lay,Int,ID,absMat] = table2struct_GP(DeviceTable,SET)
% Convert 3-column device table input from GUI or file to
% Lay and Int struct format to be passed on to GenPro_core.

% note: table is most user friendly input and output format while
% struct is most convenient format to use insinde GENPRO_core
% (e.g. for saving and reusing SM and AT data). Therefore,
% conversion table-->struct-->table is unavoidable.

% Assuming that the table does NOT contain in/outgoing medium.
% Not going to add it here, this will be done inside GENPRO_core.

% ---INPUT---
% DeviceTable   three-column table describing device
%               1 = material/texture,
%               2 = thickness,
%               3 = unit
% ---OUTPUT---
% Lay, Int      default input structures for GENPRO_core
% ID            matrix with indices of flagged layers or coatings
%               1st column: is absorber 1=yes,0=no
%               2nd column: is thickness vector 1=yes,0=no
%               3rd column: index (which Lay or Int)
%               4th column: sub-index (which coat, is 0 if Lay)
% absMat        Index of the absorber layers


% convert Thickness column from string to double (scalar or vector)
DeviceTable = string2thickness(DeviceTable);

% initialize Lay,Int structs (gives warning if not done)
nrJ = sum(~startsWith(DeviceTable.Layer,'*') & ...
    ~strcmp(DeviceTable.Unit,"nm"));     % nr incoherent layers

Lay(nrJ).med = [];        % because in/outgoing media are not speci-
Int(nrJ+1).tex = [];      % fied, nr of interfaces will nr layers +1.

[~,IsAbs] = read_rgb(DeviceTable.Layer,SET);
absMat = find(IsAbs);

% initialize ID matrix
R = size(DeviceTable,1);             % nr of table rows
ID = zeros(R,4);

lx = 1;               % current layer/coating index

for row = 1:R         % for every table row

    if startsWith(DeviceTable.Layer(row),'*')        % if texture
        Int(lx).tex = DeviceTable.Layer{row}(2:end); % set name
        % note: if the user incorrectly set multiple textures at
        % the same interface, the LAST one will overwrite any
        % previous ones. So only last texture will be used then.

        % note: not opening the texture file here to load model,xy,
        % Z-data. This will be done by GENPRO_core.

        % reduce the absMat index so that texturing is not counted as layer
        absMat(absMat > row) = absMat(absMat>row) - 1;

    else       %if is is not texture, it is real layer/coating

        % ---check ID column, extract data for ID matrix---
        id = zeros(1,2);                 % matches nr of flag types
        if IsAbs(row), id(1) = 1; end

        % check if Thickness is a row (range) or column (slicing)
        if size(DeviceTable.Thickness{row},2)>1         % if row

            id(2) = 1;  % identifies this layer thickness as range
            thickness = DeviceTable.Thickness{row};
            % thickness vector will be interpreted as range

        elseif size(DeviceTable.Thickness{row},1)>1     % if column

            thickness = slicer(DeviceTable.Thickness{row});
            % thickness vector will be interpreted as slicing

        else                                            % if scalar
            thickness = DeviceTable.Thickness{row};
        end



        % ---check unit (nm = coating, um/mm = layer)---
        % note: in GENPRO_core a 'coating' is treated coherent and
        % conformal, 'layer' is treated incoherent an nonconformal.
        % All thickness should be converted to um.

        if (DeviceTable.Unit{row} == "nm")  % treat as 'coating'

            if ~isfield(Int,'coat')         % if no coating yet
                C = 1;                      % will be first.
            else                            % if coating exists
                C = length(Int(lx).coat)+1;  % add below previous.
            end

            Int(lx).coat(C).med = DeviceTable.Layer{row}; % medium
            Int(lx).coat(C).thi = thickness/1000;       % nm to um
            Int(lx).coat(C).isabsorber = id(1);       % for a_plot
            ID(row,:) = [id,lx,C];          % add row to ID matrix

        else                                 % treat as 'layer'

            if DeviceTable.Unit{row} == "mm" % if in mm...
                Lay(lx).thi = 1000*thickness; % convert to um
            else
                Lay(lx).thi = thickness;      % already in um
            end
            Lay(lx).med = DeviceTable.Layer{row}; % material
            Lay(lx).isabsorber = id(1);       % for a_plot
            ID(row,:) = [id,lx,0];            % add row ID matrix

            lx=lx+1;       % only increment after 'layer' was added
        end
    end
end
%..................................................................
    function slicing = slicer(thickness)
        %Creates slicing vector for depth profile (GENPRO_core input).

        %---INPUT---
        % thickness   step;stop column vector (e.g. [20;180])[nm,um,mm]
        %             note: final value (180) indicates layer thickness

        %---OUTPUT---
        %thickness_vector from 0 to thickness in steps (0,20,40,...180)
        %             note: format accepted by GENPRO_core. [same unit]

        v = [0;thickness];        %add first (0) and final (thickness)

        slicing = 0;       %inital value always zero

        for i = 1:(length(v)-1)/2   %for every number pair
            sub_vector = v(2*i-1):v(2*i):v(2*i+1);  %create subvector
            if sub_vector(end)<v(2*i+1)
                sub_vector(end+1) = v(2*i+1); %#ok<AGROW>
            end
            slicing = [slicing,sub_vector(2:end)]; %#ok<AGROW>
            %do not add first value of subvector to avoid duplication
        end

    end
end


function tab = string2thickness(tab)
        % DeviceTable.Thickness is given as string. This can be just a
        % number, or contain ':' to indicate a thickness range, or contain
        % '/' to indicate slicing. This function converts string to double,
        % i.e. numerical value. Because this can be a scalar or vector, the
        % 'Thickness' variable is converted to cell array.
        % A range '40:10:90' is converted to a ROW vector [40,10,90]
        % A slicing '20/180' is converted to a COLUMN vector [20;180]


        % Note: a similar function is used by DevicePlot2D, but this does
        % not output a thickness vectors, just a scalar used for the plot.

        %---INPUT---
        % tab           table (should contain the variable 'Thickness')

        % convert 'Thickness' variable from string to cell, i.e. {string}
        tab.Thickness = arrayfun(@(dIn) {dIn},tab.Thickness);

        R = size(tab,1);             % nr of rows

        for rr = 1:R
            if contains(tab.Thickness{rr},':')       % if it is a range

                % extract start, step, stop numbers
                nrs = str2double(split(tab.Thickness{rr},':'))';  
                tab.Thickness{rr} = nrs;             % store vector in cell

            elseif contains(tab.Thickness{rr},'/')   % if it is a slicing

                % extract step, stop numbers
                nrs = str2double(split(tab.Thickness{rr},'/'));
                % set thickness to final stop
                tab.Thickness{rr} = nrs;

            else                             % else assume it's a scalar
                % scalar in cell
                tab.Thickness{rr} = str2double(tab.Thickness{rr}); 
            end
        end
    end