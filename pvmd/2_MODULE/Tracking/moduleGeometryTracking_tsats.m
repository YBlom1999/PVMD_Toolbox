function [V_final,F_final,T_final,BF,Ncells,Acell,Amod,ML,MW] = moduleGeometryTracking_tsats(TOOLBOX_input,CELL_output,BackwardTracer,Submod_ind)
% moduleGeometryTracking_hsats builds 3D module geometry of module and frame in a tilted single-axis tracking system.
%
% This function builds the 3D module geometry based on the users input.
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Output of the CELL block
% BackwardTracer: boolean
%   Boolean indicator that indicates whether the LUX (0) simulation or the
%   backward tracer (1) is used.
% Submod_ind: double
%   Index of which submodule is simulated
%
% Returns
% -------
% V_final : double
%   xyz-coordinate of every vertex
% F_final : double
%   facet matrix (every row has 3 or 4 vertex numbers that form triangle or rectangle)
% T_final : struct
%   type of reflection and transmission (every vertex belongs to one type)
% BF : Boolean
%   bifacial module (1 = yes, 0 = no)
% Ncells : Double
%   number of cells
% Acell : Double
%   area of the cell
% Amod : Double
%   area of the module
% ML : Double
%   length of the module
% MW : Double
%   width of the module
%
% Developed by T. Martens



%===primary geometry parameters (from user input)===
CR = TOOLBOX_input.Scene.module_mounting.CellRows(Submod_ind);                        %Number of cell rows
CC = TOOLBOX_input.Scene.module_mounting.CellColumns(Submod_ind);                     %Number of cell columns
MT = TOOLBOX_input.Scene.module_mounting.ModThick(Submod_ind);                        %Module thickness [cm]
CS = TOOLBOX_input.Scene.module_mounting.CellSpacing(Submod_ind);                        %Cell spacing [cm]
ES = TOOLBOX_input.Scene.module_mounting.EdgeSpacing(Submod_ind);                        %Edge spacing [cm]
TL = TOOLBOX_input.Scene.module_mounting.ModTilt;                        %Module tilt [deg]
AZ_mod = TOOLBOX_input.Scene.module_mounting.ModAzimuth;                        %Module azimuth [deg]
HG = TOOLBOX_input.Scene.module_mounting.ModMountHeight;                  %Height to ground [cm]
CL = TOOLBOX_input.Scene.module_mounting.CellLength(Submod_ind);                      %Cell length [cm]
CW = TOOLBOX_input.Scene.module_mounting.CellWidth(Submod_ind);                       %Cell width [cm]
ML = CR*CL + (CR-1)*CS + 2*ES;                  %Module length [cm] cell + intercel spacing + edge spacing
MW = CC*CW + (CC-1)*CS + 2*ES;                  %Module width [cm]
DS = TOOLBOX_input.Scene.module_mounting.ModSideSpacing + MW;             %Module side spacing [cm]
DF = TOOLBOX_input.Scene.module_mounting.ModRowSpacing;    %Module row spacing [cm]

TA = TOOLBOX_input.Scene.module_mounting.tilt_axis; 
AZ_A = TOOLBOX_input.Scene.module_mounting.azimuth_axis;
pole_W = TOOLBOX_input.Scene.module_mounting.pole_width;


ICD = TOOLBOX_input.Scene.module_mounting.inter_column_distance;    %inter column distance
NR = TOOLBOX_input.Scene.module_mounting.number_rows;     %number of rows in front
NC = TOOLBOX_input.Scene.module_mounting.number_columns;     %number of columns in front
Nmod_side = TOOLBOX_input.Scene.module_mounting.number_modules_side;   %number of modules at the side 
cell_normal = [0 0 1];

flag_side = AZ_A - 90;
if flag_side>=360; flag_side = flag_side-360;end
if flag_side<0; flag_side = flag_side+360;end

Acell=CW*CL*1e-4;
Ncells=CR*CC;
%TODO: GEO is taken as output because it contains cell width, which is
%needed for the thermal model. (really needed?). But this means that it
%might give problems for other geometries.

%TODO: should ask user for optical properties of all surfaces
if TOOLBOX_input.Scene.module_mounting.Albedo_eff
    GA = TOOLBOX_input.Scene.module_mounting.Albedo;       %Ground albedo
else
    Ground_material = TOOLBOX_input.Scene.module_mounting.Ground_material;
    wav = CELL_output.CELL_FRONT.wav;
    specRef = readSpectralReflectivity(Ground_material,wav);
    R_ground = repelem(specRef,1,2);
    A_ground = repelem(1-specRef,1,2);
    RT_ground.RAT = cat(3,R_ground,A_ground,zeros(length(wav),2));
    RT_ground.wav = wav;
    RT_ground.aoi = [0,90];
    RT_ground.lay = {'air','ground','-'};
    GA = nan;
end

%===secondary geometry parameters (fixed or calculated from primary)===
MM = 1e-2;                       %Margin between cell and module [cm]
BF = isstruct(CELL_output.CELL_REAR); %Bifacial module (1 = yes, 0 = no)
Amod=ML*MW*1e-4;
%Note: for simplicity the module is one massive block and the encapsulated
%cells are stuck on it like 'stickers'. For bifacial modules front and rear
%stickers are used, combined with teleport for transmittance.

V_cell  = [];
F_cell = [];
T_cell = [];

V_struc = [];
F_struc = [];
T_struc = [];

N_wav = length(CELL_output.CELL_FRONT.wav);
N_ang = length(CELL_output.CELL_FRONT.aoi);


shift_x = DS.*(1:Nmod_side);
shift_x = [-shift_x, 0 , shift_x];

N_modules = Nmod_side * 2 + 1;

%---make array of cell FRONT side squares---

for c = 1:CC    %for every column
    for r = 1:CR        %for every row
        [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),MT+MM,CW,CL,0); %make a square
        [V_cell,F_cell] = combine({V_cell,V},{F_cell,F},[],BackwardTracer);    %combine new V and F into main VV and FF
    end
end

T_cell(1).Facet = [(1:2:2*Ncells);(2:2:2*Ncells)];

T_cell(1).RT = CELL_output.CELL_FRONT;          %give it the corresponding optical properties

%Add layer with ones to calculate received irradiance
T_cell(1).RT.RAT = cat(3,T_cell(1).RT.RAT(:,:,1:end-1),ones(N_wav,N_ang),T_cell(1).RT.RAT(:,:,end));
T_cell(1).RT.lay = [T_cell(1).RT.lay(1:end-1);{'Full-Abs'};T_cell(1).RT.lay(end)];

%---make horizontal array of cell REAR side squares---
if BF   %if module is bifacial,
    for c = 1:CC    %for every column
        for r = 1:CR        %for every row
            %make horizontal array of cell rectangles
            [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),-MM,CW,CL,0);  %make a square
            [V_cell,F_cell] = combine({V_cell,V},{F_cell,F},[],BackwardTracer);    %combine new V and F into main VV and FF
        end
    end
    T_cell(2).Facet = T_cell(1).Facet(end) + T_cell(1).Facet;      %TYPE 2 = solar cell rear
    T_cell(2).RT = CELL_output.CELL_REAR;     %give it the corresponding optical properties

    %Add layer with ones to calculate received irradiance
    T_cell(2).RT.RAT = cat(3,T_cell(2).RT.RAT(:,:,1:end-1),ones(N_wav,N_ang),T_cell(2).RT.RAT(:,:,end)); %To calculate complete irradiance
    T_cell(2).RT.lay = [T_cell(2).RT.lay(1:end-1);{'Full-Abs'};T_cell(2).RT.lay(end)];
end

%Tilt middle of module around the x-axis
% with tilted you can still be on 2 different sides!

V_cell = moveit(V_cell, [0,-ML/2,0]);   %move symmetry axis above x axis

if AZ_mod == flag_side                     %this is the usual
    V_cell = rotate_x(V_cell,TL);  
    cell_normal = rotate_x(cell_normal,TL);  
else                  %this means that it flipped sides
    V_cell = rotate_x(V_cell,-TL);  
    cell_normal = rotate_x(cell_normal,-TL); 
end

%rotate around y axis to get tited axis
V_cell = moveit(V_cell, [+DS*Nmod_side,0,0]);
V_cell = rotate_y(V_cell, TA);
cell_normal = rotate_y(cell_normal, TA);
V_cell = moveit(V_cell, [-DS*Nmod_side,0,0]);

% correct height
V_cell = moveit(V_cell, [0,0,HG]);

for i = 1:N_modules
    VV = [];
    FF = [];

    %---make module bulk---
    [V,F] = rectabox(0,0,0,MW,ML,0);
    
    [VV,FF] = combine({VV,V},{FF,F},[],BackwardTracer);     %combine new V and F into main VV and FF
      
    %around middle of x axis and correct x_coor
    VV = moveit(VV,[0,-ML/2,0]);
    VV = moveit(VV,[shift_x(i),0,0]);

    [V_struc,F_struc] = combine({V_struc,VV},{F_struc,FF},[],BackwardTracer);
end

T_struc(1).Facet = [(1:2:2*N_modules);(2:2:2*N_modules)];    %TYPE 3 = module bulk
if BF
    T_struc(1).RT = [0.05,0];
else
    T_struc(1).RT = [0.1,0];
end

% Rotate around x axis for correct tilt

if AZ_mod == flag_side                     %this is the usual
    V_struc = rotate_x(V_struc,TL);  
else                %this means that it flipped sides
    V_struc = rotate_x(V_struc,-TL);   
end

% Rotate around y axis for tilted axis

V_struc = moveit(V_struc, [+DS*Nmod_side,0,0]);
V_struc = rotate_y(V_struc, TA);
V_struc = moveit(V_struc, [-DS*Nmod_side,0,0]);

% Put everything at the correct height

V_struc = moveit(V_struc, [0,0,HG]);

% Poles of every structure
x_coor_poles = [-DS*Nmod_side-pole_W,-DS*Nmod_side+(2*DS*Nmod_side+MW)*cosd(TA)];
y_coor_poles = - pole_W/2;

length_x = DS*Nmod_side*2+MW;

[V,F] = rectabox(x_coor_poles(1),y_coor_poles,0,pole_W,pole_W,HG-pole_W);
[V_struc,F_struc] = combine({V_struc,V},{F_struc,F},[],BackwardTracer);

[V,F] = rectabox(x_coor_poles(2),y_coor_poles,0,pole_W,pole_W,HG-pole_W+sind(TA)*length_x);
[V_struc,F_struc] = combine({V_struc,V},{F_struc,F},[],BackwardTracer);

%optical properties
Start_index = 2*N_modules;
T_struc(2).Facet = [Start_index+(1:2:2*12);Start_index+(2:2:2*12)];    
T_struc(2).RT = [0.5,0];

% Copy structure to create one column

V_col = V_struc;
F_col = F_struc;
T_col = T_struc;

if NR~=0
    % add the bar now, otherwise interference with cells
    x_coor_bar = -DS*Nmod_side;
    
    [V,F] = rectabox(0,-pole_W/2,-pole_W/2,length_x,pole_W,pole_W);
    V = rotate_x(V,TL);
    V = rotate_y(V, TA);
    V =  moveit(V, [x_coor_bar,0,HG]);
    T.Facet = [(1:2:12);(2:2:12)]; 
    T.RT = [0.5,0];
    [V_struc,F_struc,T_struc] = combine({V_struc,V},{F_struc,F},{T_struc,T},BackwardTracer);

    y_struc_rows = 1:1:NR;
    y_struc_rows = DF.*[-y_struc_rows,y_struc_rows];
    for i = 1:numel(y_struc_rows)
        V_copy =  moveit(V_struc,[0,y_struc_rows(i),0]);
        [V_col,F_col,T_col] = combine({V_col,V_copy},{F_col,F_struc},{T_col,T_struc},BackwardTracer);
    end
else
    V_col = V_struc;
    F_col = F_struc;
    T_col = T_struc;
end

% Copy to multiple columns

V_final = V_col;
F_final = F_col;
T_final = T_col;

if NC~=0
    x_struc_col = 1:1:NC;
    x_struc_col = (ICD+length_x*cosd(TA)).*[-x_struc_col,x_struc_col];
    for i=1:numel(x_struc_col)
        V_copy = moveit(V_col,[x_struc_col(i),0,0]);
        [V_final,F_final,T_final] = combine({V_final,V_copy},{F_final,F_col},{T_final,T_col},BackwardTracer);
    end
end

% Add cells

[V_final,F_final,T_final] = combine({V_cell,V_final},{F_cell,F_final},{T_cell,T_final},BackwardTracer);

%---make unit cell---

[V,F] = rectabox(-5000,-5000,0,10000,10000,0); %Make ground plane
T.Facet = [1,2];
if TOOLBOX_input.Scene.module_mounting.Albedo_eff; T.RT = [GA 0]; else; T.RT = RT_ground; end
T.Scat = [1,1];
[V_final,F_final,T_final] = combine({V_final,V},{F_final,F},{T_final,T},BackwardTracer); %combine ground with module




%---rotate ervything to get the correct azimuth axis
V_final = rotate_z(V_final,AZ_A-90);
cell_normal = rotate_z(cell_normal,AZ_A-90);

TOOLBOX_input.Scene.module_mounting.cell_normal = cell_normal;

end
