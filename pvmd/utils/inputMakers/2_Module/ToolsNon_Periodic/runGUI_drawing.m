%code for opening the GUI to make drawning for non-periodic environments
%This opens the GUI "GUI_drawning.mlapp".
%Author: Youri Blom

%Define the ground
size_x = 50;
size_y = 50;

V = [-size_x,-size_y,0;
    size_x,-size_y,0;
    size_x,size_y,0;
    -size_x,size_y,0];
F = [1,2,3;
    1,3,4];

materials = {'ground','ground'};
rgb = [0,1,0;
    0,1,0];
opa = [0.2;0.2];

[~,~,data_folder] = get_folder_structure;
libpath = fullfile(data_folder,'Material Library','SpectralReflectivityLibrary');
[reflectivity,~] = getMaterialReflectivity({'ground'},libpath);
Albedo = [reflectivity;reflectivity];
Scattering = [1;1];


%The initial condition of the ground is shown
fig = figure ;
hold on
p = patch('Faces',F,'Vertices',V(:,1:3));
set(p,'FaceVertexCData',rgb,'CDataMapping','scaled',...
            'FaceColor','flat','FaceVertexAlphaData',opa,...
            'AlphaDataMapping','none','FaceAlpha','flat')
view(30,30)
xlabel('X')
ylabel('Y')
zlabel('Z')

%The GUI is loaded
GUI_drawing(V,F,materials,fig,rgb,opa,Albedo,Scattering);