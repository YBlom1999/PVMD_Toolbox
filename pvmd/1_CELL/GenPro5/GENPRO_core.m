function [Lay,Int,out]=GENPRO_core(Lay,varargin)
% 1D optical model for solar cells with coatings and textures.
%
% This is the 'core' function that performs the optical calucations. It is 
% called by GENPRO_shell, which in turn is called by the graphical user 
% interface. It can also be called directly from the Matlab command line,
% as shown in example1.m, example2.m etc.
%
% The input is given in the form of layer and interface structures Lay and
% Int. Layers are treated incoherently (without interference) and coatings
% are treated coherently. Coatings and textures are part of the interface. 
%
% The extended net-radiation method is used to calculate the absorptance in
% each layer, as a function of wavelength and depth in the device. By 
% combining this with the solar spectrum, the photocurrents and generation
% profile are obtained.
%
% 3D scatter models 'flat', 'ray',and 'wave' are used to calculate the
% scatter matrices SM (= bidirectional scatter distribution function) for
% each interface. The'ray' model also allows saving of angle tree AT. 
% SM and AT data can be reused to reduce computation time in a for-loop,
% making GENPRO highly suitable for fast optimization.
% 
% For more information: https://doi.org/10.1109/JPHOTOV.2017.2669640
%
% ===INPUT===
% Lay           layer struct (length = nr layers), with the fields:
% .med          medium (materials file with nk data, or fixed 'N' value)
% .thi          thickness [um] (can be depth profile vector)
%
% Int           interface struct, with the fields:
% .coat         coating (with sub-fields: med/thi see Lay) empty = bare
% .tex          texture (textures file with model,xy,Z data) empty = flat    
%
% S             settings (when not given, default settings are used)
%
%===OUTPUT===
% Lay/Int.coat  same as input but with additional fields
% .abs          absorptance (%) at every wavelength (first=R, final=T)
% .cur          implied photocurrent (integrated over spectrum)
%
% out           same absorptance and current data in table form
%
%===precalculated output data that can be reused as input===
%Lay/Int.coat
% .N            refractive index for every wavelength (overrules .med)
% .isabsorber   is it an absorber layer? 1=yes, 0=no (for plotting)
%
% Int
% .SM           scatter matrix coat & tex (overrules .AT/.model/.tex/.coat)
% .AT           angle tree data (ray model only) (overrules .model/.tex)
% .model        'flat', 'ray' or 'wave' scatter model   (overrules .tex)
% .xy           lateral size of Z-matrix [um] (required for ray/wave) 
% .Z            texture height Z-matrix [um] (required for ray/wave) 
%
% ===EXAMPLE===
% clear Lay Int                    %clear variables
% 
% %---LAYERS---
% Lay(1).med = 'c-Si';             Lay(1).thi = 0:20:200;
% 
% %---INTERFACES---
% Int(1).tex = 'pyramids';      %pyramid texture
% Int(1).coat(1).med = 'ITO';      Int(1).coat(1).thi = 0.080;
% Int(1).coat(2).med = 'a-Si(p)';  Int(1).coat(2).thi = 0.005;
% Int(1).coat(3).med = 'a-Si(i)';  Int(1).coat(3).thi = 0.005;
% %---
% 
% Int(2).tex = 'pyramids';     %same pyramid texture
% Int(2).coat(1).med = 'a-Si(i)';  Int(2).coat(1).thi = 0.005;
% Int(2).coat(2).med = 'a-Si(n)';  Int(2).coat(2).thi = 0.005;
% Int(2).coat(3).med = 'Ag';       Int(2).coat(3).thi = 0.300;
% %---
% 
% S.SM = [1,1];                %save SM (of both interfaces)
% 
% [Lay,Int,tab] = GENPRO_core(Lay,Int,S);
   

%TODO: check consistency units xy, Z different models
%TODO: geo_p and geo_w data is it possible/desirable to reuse data?
%TODO: simulation log where details can be found, e.g. in 'out' variable.
%TODO: show device plot (as in GUI)
%TODO: WAVE model has been disabled, but some people might need it. Put
%back?


ok = LICENCE(nargin);                %check GenPro license

if ~ok, return; end                  %if license not ok, don't continue

[Lay,Int,set] = INPUT(Lay,varargin); %check and complete input

Int           = TEX01(Lay,Int,set);  %reads texture files, appends data

[Lay,Int]     = NK(Lay,Int,set);     %reads nk files, appends data

[out,Lay,Int] = TABULATE1(Lay,Int,set);  %creates empty output table

DevicePlot2D(unslice(out,2),set.plot_fig_input); % create plot of device

%---initialize loop variables---
W = length(set.wavelength_um);       %nr of wavelengths

QQ = zeros(length(set.incident_light),W);    %initialize fluxes matrix QQ
SM = cell(set.nr_interf,W);            %initialize scatter matrix
AT = cell(set.nr_interf,W);          %initialize AngleTree
logbook = cell(2,W);                 %initialize log (contains warnings)
%TODO: don't write log to xls file if there are no mesages

%if set.text_display==1, fprintf([repmat('.',1,W) '\n\n']);end %progress bar

if set.text_display==1      % if output to Matlab window
    fprintf([repmat('.',1,W) '\n\n'])
elseif set.text_display~=0        % else output to App command dialog
    set.text_display.Message = ...
        ['Simulating ',num2str(W),' wavelengths. This could take a minute.'];
end  

%===calculate scatter matrices and solve set of eq. to get fluxes===
if set.parallel_compute && set.pc_toolbox      %if willing and able

    parfor w=1:W                    %parallel processing wavelength loop
        [QQ(:,w),SM(:,w),AT(:,w),logbook(:,w)]=FLUX06(Lay,Int,set,w);
    end
    
else                                %if not willing or able

    for w=1:W                       %normal wavelength loop (not parallel)
        [QQ(:,w),SM(:,w),AT(:,w),logbook(:,w)]=FLUX06(Lay,Int,set,w);
    end

end
%===

Int = STORAGE(Lay,Int,SM,AT,set);   % store SM and AT data for reuse

Int = AID_COATABS(QQ,Int,set);      % split QQ into Int.AID and Int.coat.abs

[Lay,Int] = LAYABS(Lay,Int,set);    % convert AID to Lay.abs

out = TABULATE2(out,Lay,Int,set);   % add .abs data to table

[out,Lay,Int] = SPECTRUM(out,Lay,Int,set);  % calculate J, add to table
                                    
a_plot2(out,set.plot_fig_output,set.root);           % plot output in area plot

end

%==========================================================================

function ok = LICENCE(nr_inputs)
% Checks if licence is still valid. If ok, ok = 1, else ok = 0.

% Compares current date from datetime function to expiration date set here.
% Display reminder 2 weeks before expiration, message after expiration and 
% any other time when GenPro runs without input.

% TODO: If run from GUI message should go to GUI. 

%===Set customer name, affiliation, and email. Set licence end date===
licencee = 'Photovoltaic Materials and Devices group, TU Delft';
end_date = '1-Jan-2027';
%Convert m-files to p-code so customer can not see code or change date!
%===

%---message lines given as output---
line1  = 'GENPRO, by Rudi Santbergen, Delft University of Technology';
line2  = ['Copy of ',licencee];
line3a = ['Licence valid until: ',end_date];
line3b = ['Reminder: Your license will expire on ',end_date];
line3c = ['Your license expired on ',end_date];
line4  = ['To extend your license please contact ',...
                    'Rudi Santbergen (r.santbergen@tudelft.nl)'];
n = newline;
%---

%---exact message depends on licence time remaining---
remain = datetime(end_date,'Locale','en_US') - datetime('now');  %hours of licence remaining
two_weeks = duration('14:00:00:00');            %set duration for warning

if remain > two_weeks                  %if licence valid > 14 days

    if nr_inputs==0                    %message only if run without input
        ok = 0;                        %don't continue without input
        disp([line1 n line2 n line3a])
    else
        ok = 1;                        %ok, no message, just continue
    end

elseif datetime(end_date) >= datetime('now') %else if still valid

    ok = 1;                            %ok, but display reminder
    disp([line1 n line2 n line3b n line4]);

else                                   %licence expired,

    ok = 0;                            %not ok, message of expiration
    disp([line1 n line2 n line3c n line4]);

end

end

%--------------------------------------------------------------------------
function [Lay,Int,set] = INPUT(Lay,additional_input)
% Checks input variables and completes missing data. Data that is needed by
% several sub-functions is stored in settings struct set.

% ---INPUT---
% Lay                       mandatory input 1 = Lay
% additional_input          optional input 2 = Int, 3 = set

% ---OUTPUT---
% Lay                       layer struct
% Int                       interface struct
% set                       settings struct (combined from input and file)
% with the following fields:                        
% start_wavelength_nm       initial wavelength [nm]
% step_wavelength_nm        wavelength step [nm]
% stop_wavelength_nm        final wavelength [nm]
% nr_ang_interval           number of angular intervals
% incident_ang_interval     incident light angular interval
% parallel_compute          parallel computation for speed (1 =yes, 0 =no)
% plot_fig_input            plot results (0 = no, >0 = yes, axis_handle =
%                           yes, plot on that specific axis (e.g. in GUI))
% plot_fig_output           plot results (0 = no, >0 = yes, axis_handle =
%                           yes, plot on that specific axis (e.g. in GUI))
% spectrum                  name of spectrum in spectrum.xlsx file
% energy_conserv            warning threshold energy non-conservation
% ray_nr                    nr of rays (10 = fast, 1000 = accurate)
% ray_energy_conserv        warning threshold total terminated ray int.
% out_scatter_matrix        give 'ScatterMatrix' as output (1 =yes, 0 =no)
% out_ang_tree              'AngleTree' ray-data as output (1 =yes, 0 =no)
%===created here===
% pc_toolbox                parallel computing toolbox installed
% bound_ang_interval        angular interval boundaries
% coat_slices               coating cell aray. For each interface there is 
%                           a cell containing a vector. For each coating 
%                           this vector has one element indicating
%                           the number of subcoatings (for depth profiling)
% lay_slices

%TODO: if run from GUI, error messages should go to GUI
%TODO: check which layer is absorber

%---check 1st input argument = layer struct LAY---
if ~isfield(Lay,'med'),error('NO LAYER MEDIUM'); end        %mandatory
if ~isfield(Lay,'thi'),error('NO LAYER THICKNESS'); end     %mandatory

nr_input = length(additional_input);     %nr of ADDITIONAL input arguments

%---check 2nd input argument = interface struct INT---
nr_interf = length(Lay)+1;               %nr of interfaces there SHOULD be
% assuming incident medium is not given as input (will be checked later).

if nr_input == 0 || isempty(additional_input{1})   %2nd argument NOT given

    %make all interfaces bare (=no coating) and flat (=no texture)
    Int(nr_interf).coat = [];            %create empty coatings (bare)
    Int(nr_interf).tex = [];             %create empty textures (flat)

else                                               %2nd argument IS given

    Int = additional_input{1};

    %if coating not defined, add all bare (= no coating) interfaces
    if ~isfield(Int,'coat'),  Int(nr_interf).coat = []; end
    
    %if texture not defined, add all flat (= no texture) interfaces
    if ~isfield(Int,'tex'),   Int(nr_interf).tex = []; end

    %if coatings and textures are defined, but not for final interface(s)
    if length(Int)<nr_interf, Int(nr_interf).tex = []; end   %append 

end

%---check 3rd input argument = settings struct set---    
if nr_input > 1                                     % if settings input
    set = additional_input{2};                    % but no root
else                                   % if no input, all setting from file
    set = read_settings('settings.xlsx');           % setting file in current folder
    set.root = '';                     % create empty root
end
% set.root comes from GUI such that when deployed exe can find files. If
% not given (when run without GUI) root is empty means settings file
% (and others) can be found in current directory.

%---add incident and outgoing medium (as specified in settings)---

if ~strcmp(Lay(1).med,set.incident_medium)       %if first layer is not air

    fn = fieldnames(Lay);                          % get each field name
    for f = 1:length(fn), iLay.(fn{f}) = []; end   % initialize empty
    iLay.med = set.incident_medium;                % set inc. medium
    Lay = [iLay,Lay];                              % concatenate struct

end
Lay(1).thi = inf;           % set first layer to infinite thickness
Lay(1).isabsorber = 0;      % air is not absorber layer

if ~strcmp(Lay(end).med,set.incident_medium)     %if final layer is not air
    Lay(end+1).med = set.incident_medium;        %add air as final layer
end
Lay(end).thi = inf;           %set final layer to infinite thickness
Lay(end).isabsorber = 0;
%---

%---number of interfaces (after adding incident medium)---
set.nr_interf = length(Lay)-1;          %nr of interfaces goes in settings
Int(set.nr_interf+1:end) = [];         
% remove any excess interfaces in case assumption that no incident
% medium was input is incorrect.

%---make sure there is an SM and AT setting for each interface----
s = length(set.out_scatter_matrix);
if s < set.nr_interf
    if s==1          %if scalar, apply to all
        set.out_scatter_matrix = ...
            ones(1,set.nr_interf)*set.out_scatter_matrix;
    else            %else append zeros
        set.out_scatter_matrix = ...
            [set.out_scatter_matrix,zeros(1,set.nr_interf-s)];   
    end
end

a = length(set.out_ang_tree);
if a < set.nr_interf
    if a==1          %if scalar, apply to all
        set.out_ang_tree = ones(1,set.nr_interf)*set.out_ang_tree;
    else            %else append zeros
        set.out_ang_tree = [set.out_ang_tree,zeros(1,set.nr_interf-a)];   
    end
end

%===other useful data (also included in settings struct)===
%this is internal data, not from settings file or user input

%---wavelength settings---
set.wavelength_um = (set.start_wavelength_nm:set.step_wavelength_nm:...
                                            set.stop_wavelength_nm)/1000;

%---toolbox settings---
% v = ver;                              %gives version of INSTALLED toolboxes
% set.pc_toolbox = sum(strcmp({v.Name},'Parallel Computing Toolbox'));
set.pc_toolbox = contains([ver().Name],"Parallel Computing Toolbox");

if set.pc_toolbox                      %if toolbox available

    %if parallel pool does not exist, but is desired... create
    if isempty(gcp('nocreate')) && set.parallel_compute, parpool;end

    %if parellel pool exists, but is not desired... delete
    if ~isempty(gcp('nocreate')) && ~set.parallel_compute
        delete(gcp('nocreate'));
    end   
end

%---angular interval settings---
set.bound_ang_interval = linspace(0,pi/2,set.nr_ang_interval+1); 

%---coating vector and slicing settings---
% Note: layers and coating can be sliced to obtain depth-resolved profiles.
% Thickness scalar indicates no slicing, thickness vector indicates slice 
% depths (plus first value is 0 and final value is thickness)

for c1=1:set.nr_interf                          %for every interface

    l = length(Lay(c1).thi);                    %thickness vector length
    set.lay_slices(c1) = max(1,l-1);            %nr of slices in layer on top

    if isfield(Int(c1).coat,'thi')              %if coating thickness given
        %for all coatings at this interface
        l = cellfun('length',{Int(c1).coat.thi});  %thickness vector lenght
        set.coat_slices{c1} = max(1,l-1);       %nr of coating slices
    else
        set.coat_slices{c1} = [];   %no coating (or undefined thickness)
    end
end

set.lay_slices(set.nr_interf+1) = 1;            %final layer is not sliced

%---illumination settings---
tot_coatslices = sum([set.coat_slices{:}]);

%mega matrix is the main net-radiation matrix holding all linear equations.
%mega-matrix size is total number of subfluxes and coating slices.
mms = 4 * set.nr_interf * set.nr_ang_interval + tot_coatslices;  

set.incident_light = zeros(mms,1);              %vector size matches mms

if set.incident_ang_interval > 0                %if front-side illumination

    interval = set.incident_ang_interval;       %illumination interval
    set.incident_light(interval)=1;             %incident light vector

else                                            %rear-side illumination

    interval = set.nr_ang_interval * (4 * set.nr_interf - 2)...
              - set.incident_ang_interval;      %counting back from final 

    set.incident_light(interval)=1;             %incident light vector

    set.skip_dark = 0;       %sd acceleration does not work, turn off
end

    %......................................................................
    function set_f = override(set_f,set_i)
        %override default settings (from settings.xlsx file) with any 
        %settings given by user as input

        setting = fieldnames(set_i);           %get input settings names

        for c = 1:length(setting)              %for each input setting

            %if isfield(set_f,setting{c})       %check if name matches

            set_f.(setting{c}) = set_i.(setting{c}); %if yes, overwrite

            %end
        end
    end
    %......................................................................
end

%--------------------------------------------------------------------------

function Int = TEX01(Lay,Int,set)
% If a texture is specified, take data (model/xy/Z) from file. If data 
% already exists, this will be used (not the data from the file).

for i = 1:set.nr_interf                          %for every interface

    if ~isfield(Int,'SM') || isempty(Int(i).SM)
        %if SM exists, no texture data needed, ALL below will be skipped
        if ~isfield(Int,'model') || isempty(Int(i).model) 
            %if model/xy/Z data already exists, skip below 
            if isempty(Int(i).tex)
                %if also no tex-file specified, assume flat interface
                Int(i).model = 'flat';
            else                            
                %if tex-file IS specified, read its data
                [xy,Z,~,model,substrate] = read_tex(Int(i).tex,set.root);
                %TODO: use unit information

                Int(i).xy = xy;             % add xy data to Int(i) struct
                Int(i).model = model;       % add model to Int(i) struct
                
                % add Z data to Int(i) struct, accounting for orientation
                if strcmp(substrate,Lay(i).med) % if substrate is ABOVE
                    Int(i).Z = -Z;          % flip texture upside down
                else
                    Int(i).Z = Z;           % otherwise normal upright
                end

            end
        end

        % at this point we know the model for interface i. Just in case it
        % is the 'ray' model some extra steps are needed

        if strcmp(Int(i).model,'ray')

            % if AT data is available, no ray-tracing is required
            if ~isfield(Int,'AT') || isempty(Int(i).AT)    %if AT not input

                % name of file where AT is potentially stored
                %AT_file_name = ['textures\ray_data\',Lay(i).med,'_',...
                %                       Int(i).tex,'_',Lay(i+1).med,'.mat'];

                AT_file_name = fullfile("textures","ray_data",...
                    [Lay(i).med,'_',Int(i).tex,'_',Lay(i+1).med,'.mat']);

                if set.load_ang_tree && isfile(AT_file_name)    
                    % if AT found in file
                    load(AT_file_name,'ATi');        % load AT data
                    Int(i).AT = ATi;                 % add data to Int(i)
                    % TODO: check that wl range is the same

                else  % if AT (input or file) not available, do ray-tracing

                    % get geometry data (w=whole, p=parts)
                    [Int(i).geo_w,Int(i).geo_p]=...
                        pRAY02(Int(i).Z,Int(i).xy); 
                    % Note: the presence of this geometry data indicates
                    % fresh ray-tracing will be / has been done. Can be
                    % checked when saving AT data
                end
            end
        end
    end
end

end
%--------------------------------------------------------------------------

function [Lay,Int] = NK(Lay,Int,set)
% If nk-data not included, it is read from nk files and included.

% INPUT (see above)
% OUTPUT (added to input struct)
% Int.coat.N complex refractive index of coating (for every wl)
% Lay.N      complex refractive index of layer (for every wl)

I=set.nr_interf;                 % nr. of interfaces

%TODO: if interface has SM data, N won't be used
for i=1:I             %for every interface...
    for c=1:length(Int(i).coat)   %...for every coating...
        %if N not included, get N from file
        if ~isfield(Int(i).coat,'N') || isempty(Int(i).coat(c).N) 

            Int(i).coat(c).N = ...
                read_nk(Int(i).coat(c).med,set.wavelength_um);

        end
    end
end

for l=1:I+1                 %for every layer...
    if ~isfield(Lay,'N') || isempty(Lay(l).N)   %if N not included 

        Lay(l).N= read_nk(Lay(l).med,set.wavelength_um);   %get N

    end
end

end

%--------------------------------------------------------------------------
function [table_out,Lay,Int] = TABULATE1(Lay,Int,set)
% Put INPUT data from Lay,Int struct into 3 column table (Layer, Thickness,
% Unit). If a layer is divided into slices, then each slice has its own
% row.
% 
% Note: This table is input for deviceplot (which plots the device
% structure). After the simulation, the results will be added and then the
% table is input for aplot (which plots the spectral absorptance).
% The table output from GENRPO_core should be useable as input for
% GENPRO_shell and GENRPRO_gui, as first 3 layers are un
% Note: to every Lay/Int entry a row number is added, to keep track of
% where this data can be found in the table (which will be useful later).
%
% TODO: in GENPRO_shell table is converted to struct and here it is
% converted back to table again. GENPRO_core works very conveniently with
% struct (e.g. when loading and saving SM) so this approach is OK. But
% perhaps more elegant approaches can be developed. 
% TODO: when SM is given as input, that data will be used. Other input (the
% data that goes into the table) will be ignored. In principle there can be
% a mismatch between table input and results output. Solution: interface
% data should somehow be added to SM.
%
% ---INPUT---
% Lay,Int,set   layer, interface and settings structs

% ---OUTPUT---
% table_out     table with row for R, A (every slice), T and texture
% Lay,Int       table row added to struct (to later find corresponding row)

%---Determine the nr of rows for output table---
I = set.nr_interf;

has_texture = zeros(I,1);      % initialize 0

    for i = 1:I             % for every interface

        %check if it has texture from tex field
        if isfield(Int,'tex') && ~isempty(Int(i).tex)
            has_texture(i) = 1;
        end
        %check if it has texture Z-data
        if isfield (Int,'Z') && ~isempty(Int(i).Z)
            has_texture(i) = 1; 
        end

    end 

%nr of rows in table = nr of slices + nr of textures
rows = sum(set.lay_slices) + sum([set.coat_slices{:}]) +...
           sum(has_texture);

%---initialize table---
sz = [rows,3];              % three columns
varNames = {'Layer' ,'Thickness','Unit'  };   % names of columns
varTypes = {'string','string'   ,'string'};   % data types

table_out = table('Size',sz,'VariableTypes',varTypes,...
    'VariableNames',varNames);  % initialize table

row = 0;                       % row index

for i = 1:I                         % for every interface
    
    row = add_row(Lay(i),row,0);    % add layer above interface to table
    Lay(i).row = row;               % add row number to struct

    %---include row of texture data---
    if has_texture(i)               % if has texture
        row = row(end) + 1;

        if isfield(Int,'tex') && ~isempty(Int(i).tex)   % from tex file?

            table_out{row,1} = string(['*',Int(i).tex]);% use tex file name

        else

            table_out{row,1} = "*texture";       % no, unnamed texture
            % TODO: in GENPRO_core it is possible to input texture data 
            % directly, without having a texture file. Such a file is
            % needed for Deviceplot (and for reproducability).
            % Create texture.xlsx file? And include sheet in output xlsx?

        end
        Int(i).row = row;           % add row number to Int (not Int.coat) 
    end
    %---
    
    for j = 1:length(Int(i).coat)   % for every coating at this interface
        
        row = add_row(Int(i).coat(j),row,1);       % add
        Int(i).coat(j).row = row;       % add row number to struct

    end

end

row = add_row(Lay(I+1),row,0);          % final layer
Lay(I+1).row = row;                     % add row number to struct

%abs in incident medium should be interpreted as reflected/transmitted
if set.incident_ang_interval > 0        % if light incident from front
    table_out{1,1}   = "reflected";
    table_out{end,1} = "transmitted";
else                                    % if light incident from rear
    table_out{1,1}   = "transmitted";        
    table_out{end,1} = "reflected";
end

    %......................................................................
    function row = add_row(struc,row,coh)
    %add single row, from Lay/Int.coat struct, to output table.

    %---INPUT---
    % struc          Lay or Int.coat struct
    % row            previous row
    % coh            is coherent
    %
    %---OUTPUT---
    % row            current row (for sliced layers this is a vector)

    %---slicing---
    if length(struc.thi)==1                % input scalar = thickness value
        slice_thick = struc.thi;           % slice is complete layer
    else                                   % input vector = slice depth
        %convert to slice thicknesses
        slice_thick = struc.thi(2:end)-struc.thi(1:end-1);
    end

    nr = length(slice_thick);              % nr of layer/coating slices
    row = row(end) + (1:nr);               % row nr (range)

    %---Layer column---
    C = cell(nr,1);                        % in case of n slices
    C(:) = {struc.med};                    % repeat material name n times
    table_out(row,1) = C;                  % add material name to table

    %---Thickness and Unit column---
    if coh            % coherent layers can be recognized by unit being nm

        table_out{row,2} = string(slice_thick'*1000); % layer/slice thickness in nm
        table_out{row,3} = "nm";             % units nm

    elseif struc.thi(end)>1000 % incoherent and very thick
        % Layers thicker than 1000 um will be converted to mm

        if isfinite(slice_thick)
            table_out{row,2} = string(slice_thick'/1000); % layer/slice thickness in mm
        else
            table_out{row,2} = "Inf";
        end
        table_out{row,3} = "mm";

    else % incoherent and not very thick, keep unit in um

        table_out{row,2} = string(slice_thick');      % layer/slice thickness in um
        table_out{row,3} = "um";             % units um
        % TODO: why not remove the mm option completely?
    end

    %---ID column---   
%     if isfield(struc,'isabsorber') && any(struc.isabsorber) % if absorber layer
%         C = cell(nr,1);                      % in case of n slices
%         C{:} = 'a';                          % repeat 'a' n times
%         table_out(row,4) = C;                % add 'a' to 'ID' column
%     end

    %simply pass on input absorber status (usually present when input from
    %GUI, not when input from m-file). Reading missing data (and rgb) from
    %settings.xlsx file not done here, but in aplot (so only reading file 
    %when plotting.)
    end

    %......................................................................

end

%--------------------------------------------------------------------------
%CODE BELOW IS CALLED FROM PARFOR LOOP (RUNS FOR EVERY WAVELENGTH)
%--------------------------------------------------------------------------

function [Q,SMi,ATi,logbook]=FLUX06(Lay,Int,set,w)
%single wavelength complete calculation of fluxes
%
%INPUT
%Lay, Int       struct with layer and interface information
%S              struct with settings
%w              nr of current wavelength
%OUTPUT
%Q              flux vector
%SMi            scatter matrix for every interface
%ATi            angle tree for every interface

set.w=w;                                        % put current wavelength in S
MMp = diag(ones(1,length(set.incident_light)),0); MMs=MMp; %initialize mega-matrix
SMi = cell(set.nr_interf,1);
ATi = cell(set.nr_interf,1);
logbook = cell(2,1);

logbook{1} = [num2str(1000*set.wavelength_um(w)),' nm'];
%logbook{2} = 'hello';

cad_front = 0;                       %initialize cumulative absorption coef. thickness product
dum = cell(3,2,2);                   %if no light reaches and use this dummy SM instead.

for i=1:set.nr_interf                          %for every interface
    Int(i).skip_rear = 0;              %assume rear illum cannot be skipped
    if i>1
        cad_front = cad_front + 4*pi*imag(Lay(i).N(w))/set.wavelength_um(w) * Lay(i).thi(end); 
    end %unitless [1/um * um]
    
    if i==set.nr_interf    %if final interface
        cad_rear = 99;
    else
        cad_rear = cad_front + 4*pi*imag(Lay(i+1).N(w))/set.wavelength_um(w) * 2 * Lay(i+1).thi(end);
        %ad product for rear side is ad front plus twice (down & up) ad of
        %the layer underneath
    end
    %note: transmittance is checked only for layers, not coatings.
    if cad_front > 5 && set.skip_dark           %cad>5 means less than 1% of the light reaches
        %no light reaches interface, don't use models, dummy SM instead
        %TODO: reusing SM in thickness variation (example4.m) will 'expose'
        %artificial SM. Solution: generate SM with thinnest layer. Explain
        %in manual?
        dum([1,3],:,:) = {diag(0.5*ones(set.nr_ang_interval,1))};    %dimensions of rt-matrices in cells should match.
        dum(2,:,:) = {zeros(sum(set.coat_slices{i}),set.nr_ang_interval)};   %same for a-matrices
        SMi{i} = dum;
        ATi{i} = [];                           %empty AT seems OK
    else
        %light does reach interface, use models (ray, flat)
        if cad_rear > 5 && set.skip_dark %but if light does not reach rear side
            %after reflection from lower interface
            Int(i).skip_rear = 1;           %take note of that, such that models can skip rear illum
        end

        [SMi{i},ATi{i},logbook{2}] = CSM(Lay(i).N(set.w),Lay(i+1).N(set.w),Int(i),set,i); %calculate SM
       
        %TODO: interface is not illuminated from rear if layer underneath
        %is not transparent --> half of SM not needed
    end

    [TM] = CTM(Lay(i+1),set);                   %get layer transmittance matrix
    [MMs,MMp] = ATM(MMs,MMp,SMi{i},TM,set,i);   %add to mega-matrix
end
%---SOLVE SET OF EQUATIONS---
Qp = MMp\set.incident_light;
Qs = MMs\set.incident_light;
Q = (Qp+Qs)/2;      %invert matrix
%if matrix is close to singular, increase nr. of rays


% once Q is obtained, update progress

if set.text_display==1      % if output to Matlab window
    fprintf('\b|\n');       % add vertical line
elseif set.text_display~=0  % else output to App command dialog
    set.text_display.Message = ...
        ['Simulating ',num2str(1000*set.wavelength_um(w)),' nm.'];
    %TODO: this approach does NOT work when parallel compute is enabled.
    %Looked into workarounds, but is not worth the effort.
end                          

end

%==========================================================================

function [TM]=CTM(Lay,SET)
%calculate (layer) transmittance matrix

%INPUT
%Lay    struct with 1 layer info (refr. index, thickness)
%S      struct with general settings
%OUTPUT
%TM     layer transmittance matrix

aic=(SET.bound_ang_interval(2:end)+SET.bound_ang_interval(1:end-1))/2;     %angular interval centres

k=imag(Lay.N(SET.w));      %extinction coefficient
alpha=4*pi*k/SET.wavelength_um(SET.w); %absorption coefficient
t=Lay.thi(end)./cos(aic);     %effective thickness for each angle
tau=exp(-alpha*t);       %transmittance for each angle
TM=diag(tau,0);          %layer transmittance matrix
end

%--------------------------------------------------------------------------

function [SM12,AT,logbook]=CSM(N1,N2,Int,SET,i)
%calculate scatter-matrics
%INPUT
%N1,N2      refractive index on either side
%Int        struct with interface info
%S          struct with general settings
%nra,nrb    nr. of rays vector above/below (for every angle)
%OUTPUT
%sm12       3x2x2 cell containing sub-matrices
%dim1       reflectance, absorptance, transmittance
%dim2       incidence from above, below
%dim3       p-pol, a-pol


AT = [];
logbook = [];
%---select model to calculate scatter matrices---

if isfield(Int,'SM') && ~isempty(Int.SM)            %if SM given as input
    SM12=Int.SM{SET.w};                               %use that
else                                                %else, calculate
    switch Int.model
        case 'ray'          %RAY-TRACING MODEL
            [Rpa,Rsa,Rpb,Rsb,Apa,Asa,Apb,Asb,Tpa,Tsa,Tpb,Tsb,AT,logbook]=...
                RAY35at(N1,N2,Int,SET);                 %<---- using new 2step raytracing RAY35at AngleTree version !!!!!!!!
        case 'flat'         %FLAT INTERFACE
            [Rpa,Rsa,Rpb,Rsb,Apa,Asa,Apb,Asb,Tpa,Tsa,Tpb,Tsb]=...
                FLAT17(N1,N2,Int.coat,SET);
        otherwise           %error if none of the above
            error('CSM:model',[Int.model,' unknown scatter model'])
    end

    %---calculate energy loss, i.e. check COLUMN sums---
    ELpa=sum(Rpa,1)+sum(Tpa,1)+sum(Apa,1)-1;    %above, p
    ELsa=sum(Rsa,1)+sum(Tsa,1)+sum(Asa,1)-1;    %above, s
    ELpb=sum(Rpb,1)+sum(Tpb,1)+sum(Apb,1)-1;    %below, p
    ELsb=sum(Rsb,1)+sum(Tsb,1)+sum(Asb,1)-1;    %below, s

    if max(ELpa>SET.energy_conserv) || max(ELsa>SET.energy_conserv) || max(ELpb>SET.energy_conserv) || max(ELsb>SET.energy_conserv) %if too much
        logbook = [logbook,' | ',Int.model,'-model does not conserve energy!'];
        %disp([Int.model,'-model does not conserve energy!'])
        %TODO: could add energy loss values to log
        %disp(ELpa);disp(ELsa);disp(ELpb);disp(ELsb);
    end

    %---put all matrices in struct---
    SM12{1,1,1}=Rpa;  SM12{2,1,1}=Apa;  SM12{3,1,1}=Tpa;
    SM12{1,1,2}=Rsa;  SM12{2,1,2}=Asa;  SM12{3,1,2}=Tsa;
    SM12{1,2,1}=Rpb;  SM12{2,2,1}=Apb;  SM12{3,2,1}=Tpb;
    SM12{1,2,2}=Rsb;  SM12{2,2,2}=Asb;  SM12{3,2,2}=Tsb;
end

% if SET.plot_fig>2
%     %---plot matrices for debugging
%     h=figure(99);
%     set(h,'Name',['wavelength=',num2str(1000*SET.wavelength_um(SET.w)),' nm'])
%     xy_tile=ceil(sqrt(SET.nr_interf));
%     subplot(xy_tile,xy_tile,i)
% 
%     RT4p=[SM12{1,1,1},fliplr(SM12{3,1,1});...
%         flipud(SM12{3,2,1}),rot90(SM12{1,2,1},2)];    %join matrices
%     RT4s=[SM12{1,1,2},fliplr(SM12{3,1,2});...
%         flipud(SM12{3,2,2}),rot90(SM12{1,2,2},2)];    %join matrices
%     imagesc(RT4p+RT4s)                                  %average p&s pol
%     title(['i.',num2str(i),':',Int.model],'Color','red');
%     axis equal tight
%     axis off
%     if i>1, ii=i-1; elseif SET.w>1, ii=SET.nr_interf; else, ii=i; end
%     subplot(xy_tile,xy_tile,ii);
%     set(get(gca,'Title'),'Color','black')
% 
%     drawnow
% end

end

%--------------------------------------------------------------------------

function [MMs,MMp]=ATM(MMs,MMp,SM,TM,set,i)
%add sub-matrices (scatter matrices and layer transmittance matrix)
%to mega matrix
%INPUT
%MMp,MMs    mega-matrix (p&s-pol)
%SM         struct with scatter matrices
%TM         layer transmittance matrix
%S          struct with settings
%i          current interface
%OUPUT
%MMp,MMs    mega-matrix (with matrices for 1 interface added)

c = 4 * (i - 1);                 %offset constant
n = set.nr_ang_interval;         %nr of angular intervals

MMp(block(c+2),block(c+1)) = -SM{1,1,1}; %rp1
MMp(block(c+4),block(c+1)) = -SM{3,1,1}; %tp1
MMp(block(c+4),block(c+3)) = -SM{1,2,1}; %rp2
MMp(block(c+2),block(c+3)) = -SM{3,2,1}; %tp2
MMs(block(c+2),block(c+1)) = -SM{1,1,2}; %rs1
MMs(block(c+4),block(c+1)) = -SM{3,1,2}; %ts1
MMs(block(c+4),block(c+3)) = -SM{1,2,2}; %rs2
MMs(block(c+2),block(c+3)) = -SM{3,2,2}; %ts2

if ~isempty(set.coat_slices{i})      %if current interface has a coating
    if i==1                           %if it is the first interface
        y=1:sum(set.coat_slices{1});
    else
        y=1+sum([set.coat_slices{1:i-1}]):sum([set.coat_slices{1:i}]); 
    end

    y= y + 4 * set.nr_interf * n;        %absorptance matrices are appended to bottom
    %---add absorption matrices---
    MMp(y,block(c+1)) = -SM{2,1,1};  %a1p
    MMp(y,block(c+3)) = -SM{2,2,1};  %a2p
    MMs(y,block(c+1)) = -SM{2,1,2};  %a1s
    MMs(y,block(c+3)) = -SM{2,2,2};  %a2s
    %transposed a, OK? How about r and t?
end

if i < set.nr_interf      %if not final interface
    %---add next-layer transmission matrices---
    MMp(block(c+3),block(c+6)) = -TM;
    MMp(block(c+5),block(c+4)) = -TM;
    MMs(block(c+3),block(c+6)) = -TM;
    MMs(block(c+5),block(c+4)) = -TM;
    %---
end
%..................................................................
    function [x]=block(b)
        %convert block coordinates to matrix coordinates.
        x=(1:n)+n*(b-1);
    end
end

%--------------------------------------------------------------------------
%CODE BELOW IS OUTSIDE PARFOR LOOP (runs once)
%--------------------------------------------------------------------------

function Int = STORAGE(Lay,Int,SM,AT,set)
%store scatter matrix (SM) and angle tree (AT) in Int struct for output.
%Also AT can be saved to file.

for i=1:set.nr_interf                         %for every interface

    if set.out_scatter_matrix(i)              %if set to output SM

        Int(i).SM = SM(i,:);                  %append SM data to Int struct

    elseif set.out_ang_tree(i)                %if set to output AT

        ATi = AT(i,:);
        Int(i).AT = ATi;                      %append AT data to Int struct

        if set.save_ang_tree && isfield(Int,'geo_w') ... %if set to save AT
                && ~isempty(Int(i).geo_w)
            %checks if fresh ray tracing was done. (This is not the case
            %when SM or AT data were given as input.)

            %filename = ['textures\ray_data\',Lay(i).med,'_',Int(i).tex,...
            %                                      '_',Lay(i+1).med,'.mat'];
            AT_file_name = fullfile("textures","ray_data",...
                [Lay(i).med,'_',Int(i).tex,'_',Lay(i+1).med,'.mat']);

            save(AT_file_name,'ATi');            %save AT data to file
% TODO: to save the ang_tree to file, not only must save_ang_tree ==1, but
% also out_ang_tree ==1 and out_scatter_matrix == 0. Otherwise it doesn't
% work!

        end
    end
end

%---remove redundant fields that hold much data---
pw = isfield(Int,{'geo_p','geo_w'});
if pw(1), Int = rmfield(Int,'geo_p'); end
if pw(2), Int = rmfield(Int,'geo_w'); end

end
%--------------------------------------------------------------------------
function Int = AID_COATABS(QQ,Int,set)
% Append solved QQ into fluxes to Int.AID and Int.coat.abs.

% Fluxes are stored as Int.AID (angular intensity distribution) for each 
% interface. Subcoating absorptance is added to Int.coat.abs for each
% coating. This is done for all wavelengths simultaneously.
% Note: The RAT function converts Int.AID to Lay.abs.

%---INPUT---
%QQ             flux vector (AID + coat abs) for all wavelenghts
%Lay,Int,set    layer, interface and settings structs
%---OUTPUT--- (added to input struct)
%Int.AID        angular intensity distribution
%Int.coat.abs   absorbance of (sub)coating (for every wavelength)

n = set.nr_ang_interval;
I = set.nr_interf;              %nr of interfaces
x=4*I*set.nr_ang_interval;      %nr of fluxes (after this coat abs starts)

for i=1:I                       %for every interface

    %subfluxes to output struct (for AID)
    Int(i).AID(:,1,:) = QQ(block(4*(i-1)+1),:);
    Int(i).AID(:,2,:) = QQ(block(4*(i-1)+2),:);
    Int(i).AID(:,3,:) = QQ(block(4*(i-1)+3),:);
    Int(i).AID(:,4,:) = QQ(block(4*(i-1)+4),:);

    for j=1:length(set.coat_slices{i})   %for every coating

        nrz = set.coat_slices{i}(j);     %nr of depth intervals
        ix = x+(1:nrz);                  %indices of coating
        Int(i).coat(j).abs = QQ(ix,:);   %abs. depth profile in coating j
        x = ix(end);                     %increment index

    end
end

%..................................................................
    function [x]=block(b)
        %convert block coordinates to matrix coordinates.
        x = (1:n) + n * (b - 1);
        %TO DO: also used by ATM. Make shared function?
        %note: coating abs comes after final flux, so does not affect this.
    end
end
%--------------------------------------------------------------------------

function [Lay,Int] = LAYABS(Lay,Int,set)
% Calculate depth resolved Lay.abs from Int.AID (for all wavelengths)
% 
% note: First and final Lay.abs are reflectance and transmittance.
% note: (depth resolved) coating absorption was calculated in AID_coatabs.

% ---INPUT---
% Lay,Int,set   layer, interface and settings struct

% ---OUTPUT--- (added to input struct)
% Lay.abs       absorbance of layer (for every wavelength)

W = length(set.wavelength_um);            %nr of wavelengths
I = set.nr_interf;                        %nr of interfaces
aib = set.bound_ang_interval;             %angular interval boundaries
aic = (aib(1:end-1) + aib(2:end))/2;      %angular interval center [rad]

%---layer generation profile---
Lay(1).abs = sum(squeeze(Int(1).AID(:,2,:)),1);  %REFLECTANCE

for l = 2:I                                      %ABSORPTANCE (every layer)
    
    %---layer generation profile---
    % flux into LAYER l [-] for every angular interval and wavelength
    df = squeeze(Int(l-1).AID(:,4,:));    %downward from the top[-]
    uf = squeeze(Int(l).AID(:,2,:));      %upward from the bottom [-]
    
    alpha = 4*pi*imag(Lay(l).N)./set.wavelength_um;   %abs coef [1/um]
    
    % thickness vector (indicates slices)
    thi = Lay(l).thi.';                   %layer thickness (vector = slice)
    if numel(thi)==1, thi = [0,thi]; end  %#ok<AGROW> scalar to vector
    
    % OBLIQUE distances from LAYER top/bottom to interval top/bottom
    Tt = thi(1:1:end-1) * (1./cos(aic));  %LAYER TOP to interval top [um]
    Tb = thi(2:1:end) * (1./cos(aic));    %LAYER TOP to interval bot [um]

    Bt = (thi(end)-thi(1:1:end-1)) * (1./cos(aic)); %LAY BOT to i top [um]
    Bb = (thi(end)-thi(2:1:end)) * (1./cos(aic));   %LAY BOT to i bot [um]
    % note: depth profile might not be uniform, so can't simply be flipped 
    % to get distance to layer bottom.
    
    AP = zeros(length(thi)-1,W);        %initialize absorptance profile
    % contains absorptance in every spatial interval at every wavelength

    for w = 1:W                     %for every wavelength
        % applying Lambert-Beer law

        %absorptance due to downward flux
        Ad = (exp(-alpha(w) * Tt) - exp(-alpha(w) * Tb)) * df(:,w);  

        %absorptance due to upward flux
        Au = (exp(-alpha(w) * Bb) - exp(-alpha(w) * Bt)) * uf(:,w);

        %note: for the downward flux the first interface is the top 
        %interface, but for the upward flux the first interface is the
        %bottom interface.
        AP(:,w) = Ad + Au;      %total is sum of up/downward absorptances
    end

    Lay(l).abs = AP;            %append to Lay struct
end

Lay(I+1).abs = sum(squeeze(Int(end).AID(:,4,:)),1); %TRANSMITTANCE

end

%--------------------------------------------------------------------------
function table_out = TABULATE2(table_out,Lay,Int,set)
% Add absorption data to output table.
%
% ---INPUT---
% Lay,Int,set   layer, interface and settings structs

%---OUTPUT---
% table_out     table has a row for refl, layer/coating slice abs, transm
%               and a column for every wavelength.
 
%---create absorption matrix (layers/coating in oder)---
W = length(set.wavelength_um);      %nr of wavelengths
I = set.nr_interf;                  %nr of interfaces 
R = size(table_out,1);              %nr of rows in table

AA = zeros(R,W);            %initialize absorption matrix

row = 0;                            % row index

for i = 1:I                         % for every interface
    
    add_A(Lay(i));                  % layer above this interface

    if isfield(Int,'Z') && ~isempty(Int(i).Z)
        row = row(end) + 1;         % skip line for texture
    end   
    
    for j = 1:length(Int(i).coat)   % for every coating at this interface
        
        add_A(Int(i).coat(j))       % add

    end

end

add_A(Lay(I+1))

table_out{1,2} = 0;          % set inf thickness incident medium to zero
table_out{end,2} = 0;        % set inf thickness outgoing medium to zero

% TODO: is ID column needed? Yes, if table will be used as input.

%---convert to table format and combine into single table---
%TA = array2table(AA,...
%    'VariableNames',string(num2str(1000*set.wavelength_um')));

heading = strings(1,W);
for w = 1:W
    heading(w) = ['wl = ',num2str(1000*set.wavelength_um(w)),' nm'];
end

    TA = array2table(AA,'VariableNames',heading);

    table_out = [table_out,TA];

    %......................................................................
    function add_A(struc)
    %add .abs data to AA matrix and .med data to material cell.

    %--INPUT---
    %struc          Lay or Int.coat struct

    nr = max(1,length(struc.thi)-1);        %nr of layer/coating slices
    row = row(end) + (1:nr);                %row nr (range)

    AA(row,:)       = struc.abs;            %
    end

    %......................................................................

end

%--------------------------------------------------------------------------
function [tab,Lay,Int] = SPECTRUM(tab,Lay,Int,set)
% Energy conservation and calculate current by integrating over spectrum.
%
%INPUT          output table WITHOUT implied photocurrent J column
%OUTPUT         output table WITH implied photocurrent J column

wav = set.wavelength_um;                    %wavelength vector
W = length(wav);                            %nr of wavelengths

AA = table2array(tab(:,(end-W+1):end)); %extract absorptance data
%note: this are the final W columns in the table.

%---energy conservation check---
if min(min(AA))<-set.energy_conserv || max(max(AA))>1+set.energy_conserv
    %logbook{2,end+1} = 'Energy check warning: AA <0|>1';
    %disp('Energy check warning: AA <0|>1')
end
err=abs(sum(AA,1)-1);
if max(err)>set.energy_conserv
    %logbook{2,end+1} = 'R+A+T~=1';
    %disp('R+A+T~=1');
end
%---
%TODO: implement logbook

if W>1      %if multiple wavelength (otherwise integration makes no sense)

    %     %---calculate implied photocurrents---
    %     load('spectrum.mat','spec');           %load spectrum struct
    %     xx = find(strcmp({spec.name},set.spectrum),1);         %find entry
    %     if isempty(xx)
    %         logbook{2,end+1} = 'Spectrum not found.';
    %     end
    %     %abs matrix interpolated to spectum wl (extrapolation value 0!)
    %     AAs = interp1(set.wavelength_um.',AA.',spec(xx).data(:,1),'linear',0);
    %     J = AAs.'*spec(xx).data(:,2);

    % reads data from spectrum.xlsx and gives midpoint wl and interval j
    [wl_spec_um,j,status] = read_spectrum(set.spectrum,set.root);

    if status==0    % success
        %abs matrix interpolated to spectum wl (extrapolation value 0!)
        AAs = interp1(set.wavelength_um.',AA.',wl_spec_um,'linear',0);
        J = AAs.'*j;
    elseif status==1
        %logbook{2,end+1} = 'Spectrum not found.';
    else
        %logbook{2,end+1} = 'Error reading spectrum. Check units.';
    end

else
    J = zeros(size(tab,1),1);      % otherwise return 0 current
end

% add J column after Unit (column 3)
tab = addvars(tab,J,'after',3,'NewVariableNames','J [mA/cm2]');

for i = 1:set.nr_interf

    Lay(i).cur = J(Lay(i).row);        %take current from correct row
    
    for j = 1: length(Int(i).coat)

        Int(i).coat(j).cur = J(Int(i).coat(j).row);
    end

end

%table is complete, can be saved to file
if isfield(set,'results_to_file') && set.results_to_file
     % Writetable cannot format the data (e.g. make bold font, specify pre- 
     % cision). Therefore the data is written to a pre-formatted template.

    %file_template = 'results\template.xlsx';    % template has formatting
    file_template = fullfile("results","template3.xlsx");

    %file_name = ['results\',timestamp,'.xlsx'];
    file_name = fullfile("results",[timestamp,'.xlsx']);

    copyfile(file_template,file_name);          % copy template file
    
    writetable(tab,file_name,"Sheet","RAT")     % write data to copied file
end

end

%--------------------------------------------------------------------------











