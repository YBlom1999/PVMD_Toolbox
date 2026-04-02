function [V_final, F_final,materials,BF,cellCorners,normalSolarCell,Ncells,Acell,Amod,TOOLBOX_input,TT] = AddPanels_NonPeriodic(TOOLBOX_input,CELL_output,VV,FF,materials,Panel_i)
%add_panels adds the 3D module geometry of the module and frame.
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
% Developed by Youri Blom



if TOOLBOX_input.script==false
%===ask user input about module geometry===
    Q = {'Number of cell rows:',...
         'Number of cell columns:',...
         'Module thickness [cm]',...
         'Cell spacing [cm]',...
         'Edge spacing [cm]',...
         'Module tilt [deg]',...
         'Module azimuth [deg] (0=S,90=W,180=N,270=E)',...
         'Module side spacing [cm]',...
         'Module row spacing [cm]'...
         'Cell length [cm]'...
         'Cell width [cm]'...
         'x Coordinate [cm] (left bottom corner)'...
         'y Coordinate [cm] (left bottom corner)'...
         'z Coordinate [cm] (left bottom corner)'...
         'Albedo[0-1]'};
    if strcmp(CELL_output.TYPE,'T-F')
        default = {'1','170','4','0.05','3','27','0','100','600','120','0.5','0','0','0','0.2'};
    else
        default = {'12','6','0.5','0.3','1','27','0','100','800','15.675','15.675','0','0','0','0.2'};
    end
    GEO = inputdlg(Q,'MODULE GEOMETRY',1,default);
    if isempty(GEO), return, end      %if user presses 'cancel'
    
    %===primary geometry parameters (from user input)===
    %Number of cell rows
    TOOLBOX_input.Scene.module_mounting(Panel_i).CellRows=str2double(GEO{1});

    %Number of cell columns
    TOOLBOX_input.Scene.module_mounting(Panel_i).CellColumns=str2double(GEO{2});

    %Module thickness [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).ModThick=str2double(GEO{3});

    %Cell spacing [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).CellSpacing=str2double(GEO{4});

    %Edge spacing [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).EdgeSpacing=str2double(GEO{5});

    %Module tilt [deg]
    TOOLBOX_input.Scene.module_mounting(Panel_i).ModTilt=str2double(GEO{6});

    %Module azimuth [deg] (0=S,90=W,180=N,270=E)
    TOOLBOX_input.Scene.module_mounting(Panel_i).ModAzimuth=str2double(GEO{7});

    %Module side spacing [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).ModSideSpacing=str2double(GEO{8});

    %Module row spacing [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).ModRowSpacing=str2double(GEO{9});

    %Cell length [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).CellLength=str2double(GEO{10});

    %Cell width [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).CellWidth=str2double(GEO{11});
    
    %x Coordinate left bottom corner [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).xCoordinate=str2double(GEO{12});
    
    %y Coordinate left bottom corner [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).yCoordinate=str2double(GEO{13});
    
    %z Coordinate left bottom corner [cm]
    TOOLBOX_input.Scene.module_mounting(Panel_i).zCoordinate=str2double(GEO{14});

    %Albedo[0-1]
    TOOLBOX_input.Scene.module_mounting(Panel_i).Albedo=str2double(GEO{15});
    
end


%===primary geometry parameters (from user input)===
CR = TOOLBOX_input.Scene.module_mounting(Panel_i).CellRows;                        %Number of cell rows
CC = TOOLBOX_input.Scene.module_mounting(Panel_i).CellColumns;                     %Number of cell columns
MT = TOOLBOX_input.Scene.module_mounting(Panel_i).ModThick;                        %Module thickness [cm]
CS = TOOLBOX_input.Scene.module_mounting(Panel_i).CellSpacing;                        %Cell spacing [cm]
ES = TOOLBOX_input.Scene.module_mounting(Panel_i).EdgeSpacing;                        %Edge spacing [cm]
TL = TOOLBOX_input.Scene.module_mounting(Panel_i).ModTilt;                        %Module tilt [deg]
AZ = TOOLBOX_input.Scene.module_mounting(Panel_i).ModAzimuth;                        %Module azimuth [deg]
CL = TOOLBOX_input.Scene.module_mounting(Panel_i).CellLength;                      %Cell length [cm]
CW = TOOLBOX_input.Scene.module_mounting(Panel_i).CellWidth;                       %Cell width [cm]
ML = CR*CL + (CR-1)*CS + 2*ES;                  %Module length [cm] cell + intercel spacing + edge spacing
MW = CC*CW + (CC-1)*CS + 2*ES;  
Acell=CW*CL*1e-4;%Area cell [m2]

Ncells=CR*CC;
xLocation = TOOLBOX_input.Scene.module_mounting(Panel_i).xCoordinate;
yLocation = TOOLBOX_input.Scene.module_mounting(Panel_i).yCoordinate;
zlocation = TOOLBOX_input.Scene.module_mounting(Panel_i).zCoordinate;
%TODO: GEO is taken as output because it contains cell width, which is
%needed for the thermal model. (really needed?). But this means that it
%might give problems for other geometries.


MM = 1e-2;                       %Margin between cell and module [cm]
BF = isstruct(CELL_output.CELL_REAR); %Bifacial module (1 = yes, 0 = no)
Amod=ML*MW*1e-4;

N_environ_V = length(VV); %Number of vertices from the environment
N_environ_F = length(FF); %Number of faces from the environment
N_wav = length(CELL_output.CELL_FRONT.wav);
N_ang = length(CELL_output.CELL_FRONT.aoi);
TT = [];

%---make array of cell FRONT side squares---
for r = 1:CR        %for every row
    for c = 1:CC    %for every column
        [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),MT+MM,CW,CL,0); %make a square
        [VV,FF] = combine({VV,V},{FF,F},[]);    %combine new V and F into main VV and FF
    end
end
TT(1).Facet = 1:(CR*CC);        %TYPE 1 = solar cell front
TT(1).RT = CELL_output.CELL_FRONT;          %give it the corresponding optical properties


%Add layer with ones to calculate received irradiance
TT(1).RT.RAT = cat(3,TT(1).RT.RAT(:,:,1:end-1),ones(N_wav,N_ang),TT(1).RT.RAT(:,:,end));
TT(1).RT.lay = [TT(1).RT.lay(1:end-1);{'Full-Abs'};TT(1).RT.lay(end)];

    %---make horizontal array of cell REAR side squares---
if BF   %if module is bifacial, 
    for r = 1:CR        %for every row
        for c = 1:CC    %for every column
            %make horizontal array of cell rectangles
            [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),-MM,CW,CL,0);  %make a square
            [VV,FF] = combine({VV,V},{FF,F},[]);    %combine new V and F into main VV and FF
        end
    end
    TT(2).Facet = TT(1).Facet(end) + TT(1).Facet;      %TYPE 2 = solar cell rear
    TT(2).RT = CELL_output.CELL_REAR;     %give it the corresponding optical properties  
    
    %Add layer with ones to calculate received irradiance
    TT(2).RT.RAT = cat(3,TT(2).RT.RAT(:,:,1:end-1),ones(N_wav,N_ang),TT(2).RT.RAT(:,:,end));
    TT(2).RT.lay = [TT(2).RT.lay(1:end-1);{'Full-Abs'};TT(2).RT.lay(end)];
    
end

%---make module bulk---
[V,F] = rectabox(0,0,0,MW,ML,MT);
[VV,FF] = combine({VV,V},{FF,F},[]);     %combine new V and F into main VV and FF
TT(3).Facet = (1+BF)*(CR*CC)+(1:6);      %TYPE 3 = module bulk
TT(3).RT = [0.9,0];

VV(N_environ_V+1:end,:) = rotate_x(VV(N_environ_V+1:end,:),TL);   %tilt the whole module around the x-axis


%---rotate module and unit cell by azimuth angle---
VV(N_environ_V+1:end,:) = rotate_z(VV(N_environ_V+1:end,:),AZ);
VV(N_environ_V+1:end,:) = moveit(VV(N_environ_V+1:end,:),[xLocation,yLocation,zlocation]);

V_final = VV;
F_final = FF;
materials = [materials,repelem({'aluminum'},1,size(FF,1)-N_environ_F)];

if BF
    cellCorners_fr = reshape(VV(N_environ_V+1:N_environ_V+Ncells*4,:)',3,4,Ncells);
    cellCorners_rr = reshape(VV(N_environ_V+Ncells*4+1:N_environ_V+Ncells*8,:)',3,4,Ncells);
    cellCorners = {cellCorners_fr,cellCorners_rr};
else
    cellCorners = reshape(VV(N_environ_V+1:N_environ_V+Ncells*4,:)',3,4,Ncells);
end
normalSolarCell = repelem([sind(TL)*sind(AZ+180);sind(TL)*cosd(AZ+180);cosd(TL)],1,Ncells);


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
    function [V,F,T] = combine(VV,FF,TT)
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
            F = [F;FF{o}(:,[1,2,3])+v;FF{o}(:,[1,3,4])+v];    %facet matrix points to new vertex indices
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
    
end