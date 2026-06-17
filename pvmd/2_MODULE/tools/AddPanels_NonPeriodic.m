function [V_final, F_final,BF,cellCorners,normalSolarCell,Ncells,Acell,Amod,TT] = AddPanels_NonPeriodic(TOOLBOX_input,CELL_output,VV,FF,Panel_i,Submod_i)
% AddPanels_NonPeriodic adds the PV modules to the environmental geometry
% of the PV system in a non periodic simulation
%
% This function extends the simulated geometry and adds the panels to the
% environment
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   Simulation parameters
% CELL_output : struct
%   Output of the CELL block
% VV:   double
%   xyz-coordinate of every vertex at the input
% FF : double
%   facet matrix (every row has 3 vertex numbers that form triangle)
% Panel_i: double
%   index of which panel is being added
% Submod_i: double
%   index of which submodule is being added
%
% Returns
% -------
% V_final : double
%   xyz-coordinate of every vertex
% F_final : double
%   facet matrix (every row has 3 or 4 vertex numbers that form triangle or rectangle)
% BF: boolean
%   Indicators whether the cell is a bifacial device
% cell_corners : double
%   the xyz coordinates of the cell corners
% normalSolarCell: double
%   the normal vector of all solar cells
% Ncells : Double
%   number of cells
% Acell : Double
%   area of the cell
% Amod : Double
%   area of the module
% TT: struct
%   type of reflection and transmission (every vertex belongs to one type)
%
% Developed by Youri Blom



%===primary geometry parameters (from user input)===
CR = TOOLBOX_input.Scene.module_mounting(Panel_i).CellRows(Submod_i);                        %Number of cell rows
CC = TOOLBOX_input.Scene.module_mounting(Panel_i).CellColumns(Submod_i);                     %Number of cell columns
MT = TOOLBOX_input.Scene.module_mounting(Panel_i).ModThick(Submod_i);                        %Module thickness [cm]
CS = TOOLBOX_input.Scene.module_mounting(Panel_i).CellSpacing(Submod_i);                        %Cell spacing [cm]
ES = TOOLBOX_input.Scene.module_mounting(Panel_i).EdgeSpacing(Submod_i);                        %Edge spacing [cm]
TL = TOOLBOX_input.Scene.module_mounting(Panel_i).ModTilt(1);                        %Module tilt [deg]
AZ = TOOLBOX_input.Scene.module_mounting(Panel_i).ModAzimuth(1);                        %Module azimuth [deg]
CL = TOOLBOX_input.Scene.module_mounting(Panel_i).CellLength(Submod_i);                      %Cell length [cm]
CW = TOOLBOX_input.Scene.module_mounting(Panel_i).CellWidth(Submod_i);                       %Cell width [cm]
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
N_wav = length(CELL_output.CELL_FRONT.wav);
N_ang = length(CELL_output.CELL_FRONT.aoi);
TT = [];

%---make array of cell FRONT side squares---
for r = 1:CR        %for every row
    for c = 1:CC    %for every column
        [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),MT+MM,CW,CL,0); %make a square
        [VV,FF] = combine({VV,V},{FF,F},[],1);    %combine new V and F into main VV and FF
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
            [VV,FF] = combine({VV,V},{FF,F},[],1);    %combine new V and F into main VV and FF
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
[VV,FF] = combine({VV,V},{FF,F},[],1);     %combine new V and F into main VV and FF
TT(3).Facet = (1+BF)*(CR*CC)+(1:6);      %TYPE 3 = module bulk
TT(3).RT = [0.9,0];

VV(N_environ_V+1:end,:) = rotate_x(VV(N_environ_V+1:end,:),TL);   %tilt the whole module around the x-axis


%---rotate module and unit cell by azimuth angle---
VV(N_environ_V+1:end,:) = rotate_z(VV(N_environ_V+1:end,:),AZ);
VV(N_environ_V+1:end,:) = moveit(VV(N_environ_V+1:end,:),[xLocation,yLocation,zlocation]);

V_final = VV;
F_final = FF;

cellCorners = reshape(VV(N_environ_V+1:N_environ_V+Ncells*4,:)',3,4,Ncells);
if BF
    cellCornersRear = reshape(VV(N_environ_V+Ncells*4+1:N_environ_V+Ncells*8,:)',3,4,Ncells);
    cellCorners = {cellCorners,cellCornersRear};
end
normalSolarCell = repelem([sind(TL)*sind(AZ+180);sind(TL)*cosd(AZ+180);cosd(TL)],1,Ncells);

end