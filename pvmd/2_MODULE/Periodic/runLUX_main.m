function [MODULE_output, TOOLBOX_input] = runLUX_main(TOOLBOX_input, CELL_output)
%runLUX_main Runs the LUX software for the rac-tracing
%
% This function calculates the sensitivity map of the PV system
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
% MODULE_output : struct
%   Output of this module block
% TOOLBOX_input : struct
%   Simulation parameters
%
% Developed by unknown (Rudi?)


%===construct geometry vertex, facet and type based on the user input===
[V,F,T,BF,Ncells,Acell,Amod,ML,MW,TOOLBOX_input] = moduleGeometry(TOOLBOX_input,CELL_output,0,1);

if TOOLBOX_input.script==false
    Q = {'Number of Rays'};
    Recommended=25000; %Rays
    default = {num2str(Recommended)};
    TOOLBOX_input.Scene.module_mounting.NRays=str2double(inputdlg(Q,'Values',1,default));
    TOOLBOX_input.Scene.module_mounting.NRays=min(max(TOOLBOX_input.Scene.module_mounting.NRays,100),75000);% lower and upper bounds

    %TODO: is it possible to pass on CELL_output, so that it does not need to
    %be loaded for a second time inside the module builder?

    %===grouping schemes to obtain individual/average sensitivity map===
    Q ='Individual or average sensitivity map?';
    GS = listdlg('PromptString',Q,'SelectionMode','single','ListString',...
    {'Individual cell sensitivity','Average sensitivity'});
    %TODO: make it possible for users to group the cells in different ways,
    %like cells in one row or one column

    if GS == 2 
        TOOLBOX_input.Scene.module_mounting.avgSensitivity=true;
    else    
        TOOLBOX_input.Scene.module_mounting.avgSensitivity=false;
    end
    
end



disp('Ray-tracing geometry. This may take a few minutes...')

%===ray-tracing to obtain sensitivity map===
N_refinement = 3;
[Vs,Fs,As,~,zenith,azimuth] = icohemisphere(N_refinement,1);     %calculate light source angles of incidence THIS DETERMINES NR OF ANGLES!!!
A = length(Fs);                                    %nr of angles of incidence for ray-tracing
C = length(T(1).Facet);                            %nr of cells (Type1 is always cell)
sz = size(T(1).RT.RAT);                               
if numel(sz)==3, L = sz(3)-1; W = sz(1); elseif numel(sz)==2, L = sz(2)-1; W = 1; else, L = 1; W = 1; end


et = find(~cellfun(@isempty,{T.Emit}),1); %find light source in geometry
T(et).Emit(1) = TOOLBOX_input.Scene.module_mounting.NRays;                      %SET NR OF RAYS HERE!!!!
%initialize sensitivity map (-0.1 means not calculated yet)
SM_f = -0.1*ones(A,C,L,W);                      %dim1: angles, dim2: facets, dim3: layers
h_f = 99;                                 %figure handle
if BF, SM_r = -0.1*ones(A,C,L,W); h_r = 98; end %if bifacial, also initialize rear side sensitivity map

for a = 1:A                               %for every angle of incidence
    T(et).Emit(2:3) = [zenith(a),azimuth(a)];   %set azimuth and zenith angle
    SS = Lux58(V,F,T);                    %calculate sensitivity of all facets and layers for that angle
    SM_f(a,:,:,:) = SS(1:C,:,:);              %add to 4D sensitivity array (dim1: angle, dim2: facet (cell only), dim3: layer, dim4: wav)
    %plot progress sensitivity map (total absorptance averaged over all cells and all wavelengths)
    h_f = flatplot3(Vs,Fs,mean(mean(SM_f(:,:,1,:),4),2),h_f); %plot total absorptance, wavelength average
    if BF
        SM_r(a,:,:,:) = SS((C+1):2*C,:,:);
        h_r = flatplot3(Vs,Fs,mean(SM_r(:,:,1,30),2),h_r); %plot progress sensitivity map (cell average value)
    end
end

if TOOLBOX_input.Scene.module_mounting.avgSensitivity 
    SM_f = ones(size(SM_f)).*mean(SM_f,2); %use average sensitivity and keep the size of output constant
    if BF  
        SM_r = ones(size(SM_r)).*mean(SM_r,2); 
    end
end

disp('Sensitivity map completed')

%scan through all angles to make sensitivity map for every cell front and
%rear, also for every wavelength.
%This requires matrix (layer x wavelength) accumulator to be implemented in
%LUX

MODULE_output.skydome.AZA = [azimuth,zenith,As];    %skydome zenith, azimuth, area for every triangle center
MODULE_output.skydome.Vs = Vs;                      %skydome vertices (for plotting SKYmap)
MODULE_output.skydome.Fs = Fs;                      %skydome facets (for plotting SKYmap)
MODULE_output.SM_f{1} = SM_f;
MODULE_output.wav = CELL_output.CELL_FRONT.wav;     %pass on wavelength information
MODULE_output.N=Ncells;
MODULE_output.A=Acell;
MODULE_output.Amod=Amod;
% Z (L-24 as well)
MODULE_output.ML=ML*1e-2;
MODULE_output.MW=MW*1e-2;
% -
MODULE_output.ModTilt = TOOLBOX_input.Scene.module_mounting.ModTilt;
if BF, MODULE_output.SM_r{1} = SM_r; end
end
    



