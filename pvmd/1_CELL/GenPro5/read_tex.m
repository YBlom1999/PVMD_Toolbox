function [xy,Z,units,model,substrate] = read_tex(tex_file,root)
%Imports texture data from Excel texture file.

%INPUT: 
%tex_file       name of texture file (without folder or .xlsx)

%OUTPUT:
%xy             lateral xy size of height matrix [um] 1x2 matrix
%Z              texture height map (square matrix)
%model          model name (for calculating scatter matrices)

TexFilePath = fullfile(root,'textures',[tex_file,'.xlsx']);

Z    = readmatrix(TexFilePath,'Sheet','Z');       %read Z sheet data as matrix

info = readcell(TexFilePath,'Sheet','info');      %read info sheet data as cell

units     = info{find(strcmp(info(:,1),'units'),1),2};      %find unit (um,nm)
model     = info{find(strcmp(info(:,1),'model'),1),2};      %find model name
xy        = info{find(strcmp(info(:,1),'xy'),1),2};         %find xy value
substrate = info{find(strcmp(info(:,1),'substrate'),1),2};  %find substrate value

if numel(xy)==1, xy = xy*[1,1];end  %if single xy, assume square x = y
%keep the option of non-square (rectangular) height maps in future version.  

%TODO: what to do if texture file is not found
%TODO: what if info is incomplete
%TODO: could automatically decide on model based on unit (um = ray)

end
