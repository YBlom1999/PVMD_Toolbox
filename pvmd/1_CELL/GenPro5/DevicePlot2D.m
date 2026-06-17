function DevicePlot2D(tab,varargin)
% Plot 2D device cross-section, including texture, based on input table.
%
% Plot is created using a stacked area plot.Thicknesses are on log-scale 
% to make even thinnest layers visible.
% This function is called by GENPRO_core (in Matlab) and by GENPRO (GUI).

%---INPUT---
% tab        device structure in table format or name of xlsx file
% ax         optional: axes handle (where figure will be plotted)
% root       optional: root folder, where to find settings file etc.

%---settings that affect how the plot looks---
T_min = 1e-3;      %min. thickness [um] below which layer is invisibly thin
Z_min = 1e-2;      %min. texture amplitude [um] below which invisibly small
X_res = 128;       %texture resolution in x-direction
X = 40;            %x-range in plot units (not um)
%---

%---get axis handle (indicating where the plot should go)---
if nargin < 2           % if ax not given as input
    ax = 1;             % set default Figure(1) for plot
else
    ax = varargin{1};   % else use input value
end

if isnumeric(ax)        % if ax is not axis handle provided by GUI
    if ax==0            % it can either be 0
        return          % --> don't plot (silent mode) skip whole function
    else                % or a positive number
        figure(ax)      % --> indicating figure nr where to plot
        ax = gca;       % get corresponding axis handle
    end
end
cla(ax);                %clear any previous plot

%---get root (indicating where setting and texture files can be found)
if nargin < 3           % if root not defined
    root = [];          % assume root is current folder
else                    % if root is given (usually by GUI)
    root = varargin{2}; % set root 
end

%---check input file name---
if ~istable(tab)       % if table not input directly
    tab = read_device(tab,root);    % assume name of file containing table
end

%---adjusting the input table---
% GUI uses 'categoricals' that have to be converted to string
tab.Layer = string(tab.Layer);            % convert categorical to string
tab.Unit = string(tab.Unit);              % convert categorical to string

tab = fillmissing(tab,"constant",["air","0","um"]);  %fill missing data
tab = string2thickness(tab);              % convert string to number
% note: input thickness can be range, e.g. 40:5:90. In that case plot will
% be based on average value 65.

% Layers of infinite, zero or negative thickness will be removed before
% plotting. Textures may have 0 thickness but are NOT removed here.
remove_row = (tab.Thickness <= 0 | ~isfinite(tab.Thickness)) ...
    & ~startsWith(tab.Layer,'*');  % identify

tab(remove_row,:) = [];     % remove

tab.Thickness(tab.Thickness<T_min) = T_min;  %overwrites below min with min

% Incident medium is not in input table and is not plotted.
% However, if the first interface is textured, a 'blank' (unlabelled) layer
% needs to be inserted to make plotting method work.

first_tex = find(startsWith(tab.Layer,'*'),1,"first"); %1st texture
first_lay = find(strcmp(tab.Unit,'um')|...
    strcmp(tab.Unit,'mm'),1,"first");     %1st layer

if first_tex<=first_lay   %if 1st texture above 1st non-conformal layer
    tab(2:end+1,:) = tab(1:end,:);   % move all rows one down
    tab(1,:) = {'',1,'um'};          % add 'blank' layer (1 um thick)
end
%---

isTexture = startsWith(tab.Layer,'*');       % texture rows

%---thickness transformation to log scale---

% convert all thicknesses to um
tab.Thickness(tab.Unit == 'nm') = tab.Thickness(tab.Unit == 'nm')/1000;
tab.Thickness(tab.Unit == 'mm') = tab.Thickness(tab.Unit == 'mm')*1000;
%unit is NOT updated because nm indicates conformality of layer!

tab.Thick4plot = log(tab.Thickness)-log(T_min);    % TRANSFORMATION

% texture should not affect layer positioning, regardless of 'thickness'.
tab.Thick4plot(isTexture) = 0;                     % set tex thickness to 0
Zmid = cumsum(tab.Thick4plot)-0.5*tab.Thick4plot;  % mid-point for label
X_plot = linspace(0,X,X_res);                      % x domain in plot
TT_plot = tab.Thick4plot*ones(1,X_res);            % make it a matrix
%TT_plot is the matrix that will be area plotted. At this point it
%represents a flat solar cell, without any textures.
%note: coordinates are stored in row format (not column)

%---TEXTURE (load actual texture from file, display cross section)---
for i = find(isTexture)'                           % for every texture
    % note: if one interface has multiple textures, the plot will show
    % their superposition. GENPRO can use only one and will take the final 
    % one. So plot does not correspond to what will be simulated.

    tex_name = char(tab.Layer(i));                 % texture filename
    [xy,Z,~] = read_tex(tex_name(2:end),root);     % read xy and Z data
    %TODO: use unit data (now assuming um?)

    %---extract texture cross-section [um]---
    z_phys = diag(Z)';               %take 1D cross section at diagonal
    z_phys = z_phys - z_phys(1);     %shift 1st point to zero (for label)
    x_res = length(z_phys);          %resolution of texture data

    %---texture size transformation from um to plot units---
    %note: tiny textures need to be visible. Scale similar to layers.
    amp_phys = max(z_phys)-min(z_phys);    % amplitude of texture max-min
    amp_plot = log(amp_phys) - log(Z_min); % texture amplitude in plot
    sf = amp_plot / amp_phys;              % scale factor (plot/phys)

    Z_plot1 = z_phys*sf;                   % scaled z dimension
    x_plot1 = sqrt(2)*xy(1)*sf;            % scaled range
    % this is now in plot units, but still only a single period at the 
    % origninal resolution.

    %---texture replication/interpolation---
    % note: every texture potentially has a different range and resolution.
    % replication/interpolation is needed to make them all match. When
    % textures are 'tilable' (or periodic) this will not cause a jump.

    row = ceil(X/x_plot1);             % nr of times fits in larger domain
    Z_plot2 = repmat(Z_plot1,1,row);   % replicate texture r times
    X_plot2 = linspace(0,x_plot1*row,x_res*row);
    Z_plot = interp1(X_plot2,Z_plot2,X_plot,'linear',0); %interpolate
    %this sets correct range and resolution
    %TODO: if amplitude below minimum, ignore

    %---add texture info to TT_plot matrix---

    layer_above = find(~strcmp(tab.Unit(1:i-1),'nm') &...
            ~isTexture(1:i-1),1,'last');    %find last layer BEFORE texture

    layer_below = i + find(~strcmp(tab.Unit(i+1:end),'nm') &...
            ~isTexture(i+1:end),1,'first'); %find first layer AFTER texture
    % note: texture can have units 'um', but does not count as layer.

    % note: in a stacked area plot, adding a texture to one layer, will 
    % affect all layers on top (conformally) until at one layer that same 
    % texture is extracted. That layer will seem to absorb the texture
    % (non-conformal).

    % add texture to layer below
    TT_plot(layer_below,:) = TT_plot(layer_below,:) + Z_plot;

    % subtract texture to layer above
    TT_plot(layer_above,:) = TT_plot(layer_above,:) - Z_plot;
    %TODO: double check that texture is not upside down
end

TT_plot(isTexture,:) = [];            % remove texture lines
Zmid(isTexture) = [];

%=== PLOTTING===
h = area(ax,X_plot,TT_plot');         % plot layers bottom to top

%---assign colors---
[rgb,isAbsorber] = read_rgb(tab.Layer); % get rgb value of every layer
set(h,{'FaceColor'},rgb(~isTexture))      % set rgb color (except textures)

% ---add medium as y-labels---
label = strings(sum(~isTexture),1);     % initialize string array
l = 1;                                  % label index

for r = find(~isTexture).'              % for every layer (not texture)

    label(l) = tab.Layer{r};            % label is layer name

    % make label bold font \bf for absorbers
    if isAbsorber(r), label(l) = strcat("\bf ",label(l));end

    % make label blue for coatings
    if strcmp(tab.Unit{r},"nm"), label(l) =...
            strcat("\color[rgb]{0,0,0.5} ",label(l)); end
    % note: italic \it or \textit doesn't work in GUI.

    l = l+1;                            % increment label index
end

% two zero-thickness layers in a row will cause overlapping labels
incr = Zmid(2:end)-Zmid(1:end-1);   % check increment
Zmid(incr<=0) = [];                 % zero incr means remove y-tick
label(incr<=0) = [];                % and corresponding label

ax.YTickLabel = label;              % set y-label
%---

ax.XTick = [];
ax.YTick = Zmid;
set(ax,'TickDir','out');
axis(ax,'equal')
axis(ax,'tight')
axis(ax,'ij')                     % flips y-axis! Area stacks top to bot.
drawnow

%TODO: check thickness value is present and positive nr
%TODO: flip texture when at rear

%--------------------------------------------------------------------------
    function tab = string2thickness(tab)
        % DeviceTable.Thickness is given as string. This can be just a
        % number, or contain ':' to indicate a thickness range, or contain
        % '/' to indicate slicing. This function converts the string to a 
        % single number, which is needed to create the plot.

        %---INPUT/OUTPUT---
        % tab           table (only 'Thickness' column is changed)

        R = size(tab,1);             % nr of rows

        for rr = 1:R                 % for every table row
            if contains(tab.Thickness(rr),':')   % if it is a range
                % extract start, step, stop numbers
                nrs = str2double(split(tab.Thickness(rr),':')); 
                % set thickness to AVERAGE value of start and stop.
                tab.Thickness2(rr) = (nrs(1)+nrs(end))/2;   

            elseif contains(tab.Thickness(rr),'/') % if it is a slicing
                % extract step, stop numbers
                nrs = str2double(split(tab.Thickness(rr),'/'));
                % set thickness to final STOP value
                tab.Thickness2(rr) = nrs(end);

            else     % if not range or slicing, assume it's a single nr

                tab.Thickness2(rr) = str2double(tab.Thickness{rr});

            end
        end
        % extra column 'Thickness2' has been created
        tab.Thickness = tab.Thickness2; % this replaces the original
        tab.Thickness2 = [];            % and is then removed
    end
end