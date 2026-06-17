function [V_final,F_final,T_final,BF,Ncells,Acell,Amod,ML,MW] = moduleGeometry(TOOLBOX_input,CELL_output,BackwardTracer,Submod_ind)
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
% Developed by unknown


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
DF = TOOLBOX_input.Scene.module_mounting.ModRowSpacing + ML*cosd(TL);    %Module row spacing [cm]
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
    R_ground = repelem(specRef',1,2);
    A_ground = repelem(1-specRef',1,2);
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
    for c = 1:CC    %for every column
        for r = 1:CR        %for every row
            [V,F] = rectabox(ES+(c-1)*(CW+CS),ES+(r-1)*(CL+CS),MT+MM,CW,CL,0); %make a square
            [VV,FF] = combine({VV,V},{FF,F},[],BackwardTracer);    %combine new V and F into main VV and FF
        end
    end
    if ~BackwardTracer %For LUX simulation
        T_final(1).Facet = 1:Ncells;        %TYPE 1 = solar cell front
    else %For Backward tracer simulations
        Start_index = 2*(Ncells*(1+BF)+30)*(i-1);
        T_final(1+4*(i-1)).Facet = [Start_index+(1:2:2*Ncells);Start_index+(2:2:2*Ncells)];
    end
    T_final(1+4*(i-1)).RT = CELL_output.CELL_FRONT;          %give it the corresponding optical properties

    %Add layer with ones to calculate received irradiance
    T_final(1+4*(i-1)).RT.RAT = cat(3,T_final(1+4*(i-1)).RT.RAT(:,:,1:end-1),ones(N_wav,N_ang),T_final(1+4*(i-1)).RT.RAT(:,:,end));
    T_final(1+4*(i-1)).RT.lay = [T_final(1+4*(i-1)).RT.lay(1:end-1);{'Full-Abs'};T_final(1+4*(i-1)).RT.lay(end)];

    %---make horizontal array of cell REAR side squares---
    if BF   %if module is bifacial,
        for c = 1:CC    %for every column
            for r = 1:CR        %for every row
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
        T_final(3).Facet = (1+BF)*(Ncells)+(1:6);        %TYPE 3 = module bulk
    else %For Backward tracer simulations
        Start_index = 2*Ncells*(1+BF)*i+2*30*(i-1);
        T_final(3+4*(i-1)).Facet = [Start_index+(1:2:12);Start_index+(2:2:12)];    %TYPE 3 = module bulk
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
        T_final(4).Facet = (1+BF)*Ncells+(7:30);    %TYPE 4 = frame
    else %For Backward tracer simulations
        Start_index = 2*Ncells*(1+BF)*i+2*30*(i-1);
        T_final(4+4*(i-1)).Facet = [Start_index+(13:2:60);Start_index+(14:2:60)];    %TYPE 4 = frame
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

if TOOLBOX_input.Scene.plotScene
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
end