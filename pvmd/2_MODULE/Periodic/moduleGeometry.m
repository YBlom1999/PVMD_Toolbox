function [V_final,F_final,T_final,BF,Ncells,Acell,Amod,ML,MW, TOOLBOX_input] = moduleGeometry(TOOLBOX_input,CELL_output,BackwardTracer,Submod_ind)
%moduleGeometry Builds 3D module geometry of module and frame.
%
% This function builds the 3D module geometry based on the users input.
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Output of the CELL block
%
% Returns
% -------
% VV : double
%   xyz-coordinate of every vertex
% FF : double
%   facet matrix (every row has 3 or 4 vertex numbers that form triangle or rectangle)
% TT : struct
%   type of reflection and transmission (every vertex belongs to one type)
% BF : Boolean
%   bifacial module (1 = yes, 0 = no)
% Ncells : Double
%   number of cells
% Acell : Double
%   area of the cell
% Amod : Double
%   area of the module
%
% Developed by unknown



if TOOLBOX_input.script==false
    %===ask user input about module geometry===
    Q = {'Number of cell rows:',...
        'Number of cell columns:',...
        'Module thickness [cm]',...
        'Cell spacing [cm]',...
        'Edge spacing [cm]',...
        'Module tilt [deg]',...
        'Module azimuth [deg] (0=S,90=W,180=N,270=E)',...
        'Height to ground [cm]',...
        'Module side spacing [cm]',...
        'Module row spacing [cm]'...
        'Cell length [cm]'...
        'Cell width [cm]'};
    if strcmp(CELL_output.TYPE,'T-F')
        default = {'1','170','4','0.05','3','27','0','50','100','600','120','0.5'};
    else
        default = {'12','6','0.5','0.3','1','27','0','50','100','800','15.675','15.675'};
    end
    GEO = inputdlg(Q,'MODULE GEOMETRY',1,default);
    if isempty(GEO), return, end      %if user presses 'cancel'

    %Ask for choice of Albedo
    Q ={'How do you want to consider the Albedo reflection?'};
    Choice_Albedo = listdlg('PromptString',Q,'SelectionMode','single','ListString',{'Fixed albedo','Choose material with spectral albedo'},'ListSize',[300,100]); %user choice
    if Choice_Albedo == 1
        TOOLBOX_input.Scene.module_mounting.Albedo_eff = 1;
        Q = {'Albedo [-]'};
        default = {'0.2'};
        TOOLBOX_input.Scene.module_mounting.Albedo = str2double(inputdlg(Q,'Albedo',1,default));
    elseif Choice_Albedo == 2
        TOOLBOX_input.Scene.module_mounting.Albedo_eff = 0;
        current_path = pwd;
        [~,~,data_folder] = get_folder_structure;
        cd(fullfile(data_folder, 'Material Library','SpectralReflectivityLibrary')); %go to 'Sim' folder (where simulations are stored)
        list = dir('*.mat');              %list the file names there
        Q ='Choose material:';    %ask user
        A = listdlg('PromptString',Q,'SelectionMode','single','ListString',{list.name},'ListSize',[200,100]); %user choice
        TOOLBOX_input.Scene.module_mounting.Ground_material = list(A).name;
        cd(current_path)
    else
        return
    end
    %===primary geometry parameters (from user input)===
    %Number of cell rows
    TOOLBOX_input.Scene.module_mounting.CellRows=str2double(GEO{1});

    %Number of cell columns
    TOOLBOX_input.Scene.module_mounting.CellColumns=str2double(GEO{2});

    %Module thickness [cm]
    TOOLBOX_input.Scene.module_mounting.ModThick=str2double(GEO{3});

    %Cell spacing [cm]
    TOOLBOX_input.Scene.module_mounting.CellSpacing=str2double(GEO{4});

    %Edge spacing [cm]
    TOOLBOX_input.Scene.module_mounting.EdgeSpacing=str2double(GEO{5});

    %Module tilt [deg]
    TOOLBOX_input.Scene.module_mounting.ModTilt=str2double(GEO{6});

    %Module azimuth [deg] (0=S,90=W,180=N,270=E)
    TOOLBOX_input.Scene.module_mounting.ModAzimuth=str2double(GEO{7});

    %Module height above ground [cm]
    TOOLBOX_input.Scene.module_mounting.ModMountHeight=str2double(GEO{8});

    %Module side spacing [cm]
    TOOLBOX_input.Scene.module_mounting.ModSideSpacing=str2double(GEO{9});

    %Module row spacing [cm]
    TOOLBOX_input.Scene.module_mounting.ModRowSpacing=str2double(GEO{10});

    %Cell length [cm]
    TOOLBOX_input.Scene.module_mounting.CellLength=str2double(GEO{11});

    %Cell width [cm]
    TOOLBOX_input.Scene.module_mounting.CellWidth=str2double(GEO{12});

end

%===primary geometry parameters (from user input)===
CR = TOOLBOX_input.Scene.module_mounting.CellRows(Submod_ind);                        %Number of cell rows
CC = TOOLBOX_input.Scene.module_mounting.CellColumns(Submod_ind);                     %Number of cell columns
MT = TOOLBOX_input.Scene.module_mounting.ModThick(Submod_ind);                        %Module thickness [cm]
CS = TOOLBOX_input.Scene.module_mounting.CellSpacing(Submod_ind);                        %Cell spacing [cm]
ES = TOOLBOX_input.Scene.module_mounting.EdgeSpacing(Submod_ind);                        %Edge spacing [cm]
TL = TOOLBOX_input.Scene.module_mounting.ModTilt;                        %Module tilt [deg]
AZ = TOOLBOX_input.Scene.module_mounting.ModAzimuth;                        %Module azimuth [deg]
HG = TOOLBOX_input.Scene.module_mounting.ModMountHeight;                  %Height to ground [cm]
CL = TOOLBOX_input.Scene.module_mounting.CellLength(Submod_ind);                      %Cell length [cm]
CW = TOOLBOX_input.Scene.module_mounting.CellWidth(Submod_ind);                       %Cell width [cm]
ML = CR*CL + (CR-1)*CS + 2*ES;                  %Module length [cm] cell + intercel spacing + edge spacing
MW = CC*CW + (CC-1)*CS + 2*ES;                  %Module width [cm]
DS = TOOLBOX_input.Scene.module_mounting.ModSideSpacing + MW;             %Module side spacing [cm]
DF = TOOLBOX_input.Scene.module_mounting.ModRowSpacing + ML*cosd(AZ);    %Module row spacing [cm]
Acell=CW*CL*1e-4;
Ncells=CR*CC;
%TODO: GEO is taken as output because it contains cell width, which is
%needed for the thermal model. (really needed?). But this means that it
%might give problems for other geometries.

%TODO: should ask user for optical properties of all surfaces
if TOOLBOX_input.Scene.module_mounting.Albedo_eff
    GA = TOOLBOX_input.Scene.module_mounting.Albedo;       %Ground albedo
else
    [~,~,data_folder] = get_folder_structure;
    Ground_material = TOOLBOX_input.Scene.module_mounting.Ground_material;
    lib_path = fullfile(data_folder, 'Material Library','SpectralReflectivityLibrary',Ground_material);
    load(lib_path,'lambda','specRefl');
    wav = CELL_output.CELL_FRONT.wav;
    R_ground = repelem(interp1(lambda,specRefl,wav*1e3),1,2);
    A_ground = repelem(interp1(lambda,1-specRefl,wav*1e3),1,2);
    RT_ground.RAT = cat(3,R_ground,A_ground,zeros(length(wav),2));
    RT_ground.wav = wav;
    RT_ground.aoi = [0,90];
    RT_ground.lay = {'air','ground','-'};
    GA = nan;
end
GH = 1;                          %Ground haze
GD = 1;                          %Ground diffuse exponent

%===secondary geometry parameters (fixed or calculated from primary)===
FW = 5;                          %Frame width [cm]
MM = 1e-2;                       %Margin between cell and module [cm]
BF = isstruct(CELL_output.CELL_REAR); %Bifacial module (1 = yes, 0 = no)
Amod=ML*MW*1e-4;
%Note: for simplicity the module is one massive block and the encapsulated
%cells are stuck on it like 'stickers'. For bifacial modules front and rear
%stickers are used, combined with teleport for transmittance.



V_final = [];
F_final = [];
T_final = [];

N_wav = length(CELL_output.CELL_FRONT.wav);
N_ang = length(CELL_output.CELL_FRONT.aoi);

if ~BackwardTracer %For LUX simulation
    N_modules = 1;
    shift_x = [0,0];
    shift_y = [0,0];
else %For Backward tracer simulations
    N_modules = 9;
    shift_x = [0,-DS,-DS,-DS,0,0,DS,DS,DS];
    shift_y = [0,-DF,0,DF,-DF,DF,-DF,0,DF];
end
for i = 1:N_modules
    VV = [];
    FF = [];
    %---make array of cell FRONT side squares---
    for r = 1:CR        %for every row
        for c = 1:CC    %for every column
            [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),MT+MM,CW,CL,0); %make a square
            [VV,FF] = combine({VV,V},{FF,F},[],BackwardTracer);    %combine new V and F into main VV and FF
        end
    end
    if ~BackwardTracer %For LUX simulation
        T_final(1+4*(i-1)).Facet = (Ncells+30)*(i-1)+(1:Ncells);        %TYPE 1 = solar cell front
    else %For Backward tracer simulations
        T_final(1+4*(i-1)).Facet = [(Ncells+30)*(i-1)+(1:2:2*Ncells);(Ncells+30)*(i-1)+(2:2:2*Ncells)];
    end
    T_final(1+4*(i-1)).RT = CELL_output.CELL_FRONT;          %give it the corresponding optical properties

    %Add layer with ones to calculate received irradiance
    T_final(1+4*(i-1)).RT.RAT = cat(3,T_final(1+4*(i-1)).RT.RAT(:,:,1:end-1),ones(N_wav,N_ang),T_final(1+4*(i-1)).RT.RAT(:,:,end));
    T_final(1+4*(i-1)).RT.lay = [T_final(1+4*(i-1)).RT.lay(1:end-1);{'Full-Abs'};T_final(1+4*(i-1)).RT.lay(end)];

    %---make horizontal array of cell REAR side squares---
    if BF   %if module is bifacial,
        for r = 1:CR        %for every row
            for c = 1:CC    %for every column
                %make horizontal array of cell rectangles
                [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),-MM,CW,CL,0);  %make a square
                [VV,FF] = combine({VV,V},{FF,F},[],BackwardTracer);    %combine new V and F into main VV and FF
            end
        end
        T_final(2+4*(i-1)).Facet = T_final(1+4*(i-1)).Facet(end) + T_final(1+4*(i-1)).Facet;      %TYPE 2 = solar cell rear
        T_final(2+4*(i-1)).RT = CELL_output.CELL_REAR;     %give it the corresponding optical properties

        %Add layer with ones to calculate received irradiance
        T_final(2+4*(i-1)).RT.RAT = cat(3,T_final(2+4*(i-1)).RT.RAT(:,:,1:end-1),ones(N_wav,N_ang),T_final(2+4*(i-1)).RT.RAT(:,:,end)); %To calculate complete irradiance
        T_final(2+4*(i-1)).RT.lay = [T_final(2+4*(i-1)).RT.lay(1:end-1);{'Full-Abs'};T_final(2+4*(i-1)).RT.lay(end)];
    end

    %---make module bulk---
    [V,F] = rectabox(0,0,0,MW,ML,MT);
    [VV,FF] = combine({VV,V},{FF,F},[],BackwardTracer);     %combine new V and F into main VV and FF
    if ~BackwardTracer %For LUX simulation
        T_final(3+4*(i-1)).Facet = (i+BF)*(Ncells)+30*(i-1)+(1:6);        %TYPE 3 = module bulk
    else %For Backward tracer simulations
        T_final(3+4*(i-1)).Facet = [(i+BF)*(2*Ncells)+30*(i-1)+(1:2:12);(i+BF)*(2*Ncells)+30*(i-1)+(2:2:12)];%TYPE 3 = module bulk
    end
    T_final(3+4*(i-1)).RT = [0.9,0];

    VV = rotate_x(VV,TL);   %tilt the whole module around the x-axis
    VV = moveit(VV,[shift_x(i),shift_y(i),HG]);%move the whole module in z-direction (upward)

    %---make frame (vertical bar outside each module corner)---
    [V1,F1] = rectabox(-FW+shift_x(i),        -FW+shift_y(i),0,FW,FW,HG);
    [V2,F2] = rectabox( MW+shift_x(i),        -FW+shift_y(i),0,FW,FW,HG);
    [V3,F3] = rectabox(-FW+shift_x(i),cosd(TL)*ML+shift_y(i),0,FW,FW,HG+sind(TL)*ML);
    [V4,F4] = rectabox( MW+shift_x(i),cosd(TL)*ML+shift_y(i),0,FW,FW,HG+sind(TL)*ML);

    [VV,FF] = combine({VV,V1,V2,V3,V4},{FF,F1,F2,F3,F4},[],BackwardTracer);
    if ~BackwardTracer %For LUX simulation
        T_final(4+4*(i-1)).Facet = (i+BF)*Ncells+30*(i-1)+(7:30);    %TYPE 4 = frame
    else %For Backward tracer simulations
        T_final(4+4*(i-1)).Facet = [(i+BF)*2*Ncells+30*(i-1)+(15:2:60);(i+BF)*2*Ncells+30*(i-1)+(16:2:60)];    %TYPE 4 = frame
    end
    T_final(4+4*(i-1)).RT = [0.5,0];
    [V_final,F_final] = combine({V_final,VV},{F_final,FF},[],BackwardTracer);
end

%---make unit cell---
if ~BackwardTracer %For LUX simulation
    [V,F,T] = unitcell(V_final,DS,DF,MM,GA,GH,GD); %create unit cell around it
    if ~TOOLBOX_input.Scene.module_mounting.Albedo_eff; T(1).RT = RT_ground; end %Change RT structure for ground in a spectral reflectivity is used
    [V_final,F_final,T_final] = combine({V_final,V},{F_final,F},{T_final,T},BackwardTracer); %combine unit cell with module
else %For Backward tracer simulations
    [V,F] = rectabox(-5000,-5000,0,10000,10000,0); %Make ground plane
    T.Facet = [1,2];
    if TOOLBOX_input.Scene.module_mounting.Albedo_eff; T.RT = [GA 0]; else; T.RT = RT_ground; end
    T.Scat = [1,1];
    [V_final,F_final,T_final] = combine({V_final,V},{F_final,F},{T_final,T},BackwardTracer); %combine ground with module
end

%---rotate module and unit cell by azimuth angle---
V_final = rotate_z(V_final,AZ);

if ~BackwardTracer %For LUX simulation
    Lux58(V_final,F_final,T_final); %LUX is just called to make a plot. Nr of rays is set to 0, so no rays are traced.
else %For Backward tracer simulations
    N_face = size(F_final,1);
    rgb = zeros(N_face,3);
    opa = 0.2*ones(N_face,1);
    p = patch('Faces',F_final,'Vertices',V_final);

    set(p,'FaceVertexCData',rgb,'CDataMapping','scaled',...
        'FaceColor','flat','FaceVertexAlphaData',opa,...
        'AlphaDataMapping','none','FaceAlpha','flat')
    view(30,30)
    %     compass(VV);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    xlim([-1500,1500])
    ylim([-1500,1500])
    zlim([0,500])
end

end
%==========================================================================
function [V,F] = rectabox(x,y,z,u,v,w)
%create a rectangular box with corner coordinates xyz and length
%width height uvw. It is assumed that uvw are positive. Facets:
%bot, north, east, south, west, top have normal pointing outward.
%If u==0 or v==0 or w==0 a 2D rectangle is created.

if u == 0 %make rectangle perpendicular to x-axis
    V = [x y z; x y+v z; x y+v z+w; x y z+w];
    F = [1 2 3 4]; %normal pointing in positive x-direction
elseif v == 0 %make rectangle perpendicular to y-axis
    V = [x y z; x y z+w; x+u v z+w; x+u y z];
    F = [1 2 3 4]; %normal pointing in positive y-direction
elseif w == 0 %make rectangle perpendicular to z-axis
    V = [x y z; x+u y z; x+u y+v z; x y+v z];
    F = [1 2 3 4]; %normal pointing in positive z-direction
else    %make real 3D box
    %define vertices (start at uvw)
    V = [x y z; x y+v z; x+u y+v z; x+u y z; x y z+w; x y+v z+w; x+u y+v z+w; x+u y z+w];
    %facets (bot, north, east, south, west, top) normal pointing outward
    F = [1 2 3 4; 6 7 3 2; 7 8 4 3; 8 5 1 4; 5 6 2 1; 8 7 6 5];
end
end
%--------------------------------------------------------------------------
function V = rotate_x(V0,ang)
%rotate vertex points around x-axis by angle ang
V(:,1) = V0(:,1);
V(:,2) = cosd(ang) * V0(:,2) - sind(ang) * V0(:,3);
V(:,3) = sind(ang) * V0(:,2) + cosd(ang) * V0(:,3);
end
%--------------------------------------------------------------------------
function V = rotate_z(V0,ang)
%rotate vertex points around z-axis by angle ang
%Here positive angle is defined as clockwise (0 = South)
V(:,1) =  cosd(ang) * V0(:,1) + sind(ang) * V0(:,2);
V(:,2) = -sind(ang) * V0(:,1) + cosd(ang) * V0(:,2);
V(:,3) = V0(:,3);
end
%--------------------------------------------------------------------------
function V = moveit(V,xyz)
%translate (move) vertex point in xyz-direction
V(:,1) = V(:,1) + xyz(1);
V(:,2) = V(:,2) + xyz(2);
V(:,3) = V(:,3) + xyz(3);
end
%--------------------------------------------------------------------------
function [V,F,T] = combine(VV,FF,TT,BackwardTracer)
%combine vertex (V), facet (F) and type (T) of different objects
%Vertex and facet matrices can be concatenated, but this changes
%their numbering. So the facet matrix and T.Facet has to be
%renumbered, respectively. Giving TT = [] skips generation of T.

nro = length(VV);       %nr of objects

V = VV{1};              %first object part of combined vertex matrix
v = size(VV{1},1);      %vertex shift nr 1
F = FF{1};              %first object part of combined facet matrix

if ~isempty(TT)
    f = size(FF{1},1);  %face shift nr 1
    T = TT{1};              %first object part of combined type struct
else
    T = [];
end

for o = 2:nro           %for every object
    V = [V;VV{o}];      %append vertices (TODO: check for duplicates?)
    if BackwardTracer && size(FF{o},2) == 4
        F = [F;FF{o}(:,[1,2,3])+v;FF{o}(:,[1,3,4])+v];    %facet matrix points to new vertex indices
    else
        F = [F;FF{o}+v];    %facet matrix points to new vertex indices
    end
    v = v+size(VV{o},1); %increment vertex shift
    if ~isempty(TT)
        l = length(TT{o});
        for t = 1:l         %for every type
            TT{o}(t).Facet = TT{o}(t).Facet+f;
        end
        T = concatstruct(T,TT{o}); %add Type structure
        f = f+size(FF{o},1); %increment facet shift
    end
end
%..................................................................
    function S = concatstruct(S1,S2)
        %combine 2 structs which may have different sets of fields

        ff = unique([fieldnames(S1);fieldnames(S2)]);   %get ALL fieldnames
        for n = 1:length(ff)                            %for every name
            if ~isfield(S1,ff(n)), S1(1).(ff{n})=[]; end   %if need add to S1
            if ~isfield(S2,ff(n)), S2(1).(ff{n})=[]; end   %if need add to S2
        end
        S = [S1,S2];  %standard concat requires same set of field names
    end
end
%--------------------------------------------------------------------------
function [V,F,T] = unitcell(Vm,DS,DF,MM,GA,GH,GD)
%create unit cell around module with periodic boundary condtions
%and ceiling as light source.

xmid = (max(Vm(:,1)) + min(Vm(:,1)))/2;     %mid-point x-coord
xmin = xmid - DS/2;                         %min x

ymid = (max(Vm(:,2)) + min(Vm(:,2)))/2;     %mid-point y-coord
ymin = ymid - DF/2;                         %min y

zmin = 0;                                   %ground is always at zero
zmax = max(Vm(:,3))+MM;                     %max z

[V,F] = rectabox(xmin,ymin,zmin,DS,DF,zmax);
%TODO: check whether box is larger than the module. User can make
%box arbitraly small.

%---floor---
T(1).Facet = 1;
T(1).RT = [GA 0];         %reflected albedo
T(1).Scat = [GH,GD];      %haze and diffuse exponent
%---front-back walls---
T(2).Facet = [2 4];
T(2).RT = [0 1];          %transparent
T(2).Teleport = -DF;      %with periodic boundary conditions
T(2).Invis = 1;           %make invisible in plot
%---left-right walls---
T(3).Facet = [3 5];
T(3).RT = [0 1];          %transparent
T(3).Teleport = -DS;      %with periodic boundary conditions
T(3).Invis = 1;           %make invisible in plot
%---ceiling---
T(4).Facet = 6;
T(4).RT = [0 0];          %absorbing
T(4).Emit = [0 0 0];      %light source: theta, phi and nr rays.
T(4).Invis = 1;           %make invisible in plot
end