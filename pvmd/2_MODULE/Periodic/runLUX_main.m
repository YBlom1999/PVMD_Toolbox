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
% Developed by unknown (Rudi?), Adjusted by Y. Blom

N_submodules = max(CELL_output.SUBMOD_IND);

SM_f = cell(1,N_submodules);
SM_r = cell(1,N_submodules);
Ncells = zeros(1,N_submodules);
Acell = zeros(1,N_submodules);
Amod = zeros(1,N_submodules);

%If there is more than 1 submodule, the CELL_output needs to be adjusted
for Submod_i = 1:N_submodules

    CELL_adjusted = CELL_output;
    Submod_ind = find([1,CELL_output.SUBMOD_IND == Submod_i,1]);
    CELL_adjusted.CELL_FRONT.RAT = CELL_output.CELL_FRONT.RAT(:,:,Submod_ind);
    CELL_adjusted.CELL_FRONT.lay = CELL_output.CELL_FRONT.lay(Submod_ind);
    if isstruct(CELL_output.CELL_REAR)
        CELL_adjusted.CELL_REAR.RAT = CELL_output.CELL_REAR.RAT(:,:,Submod_ind);
        CELL_adjusted.CELL_REAR.lay = CELL_output.CELL_REAR.lay(Submod_ind);
    end

    %===construct geometry vertex, facet and type based on the user input===
    [V,F,T,BF,Ncells(Submod_i),Acell(Submod_i),Amod(Submod_i),ML,MW] = moduleGeometry(TOOLBOX_input,CELL_adjusted,0,Submod_i);


    disp('Ray-tracing geometry. This may take a few minutes...')

    %===ray-tracing to obtain sensitivity map===
    [Vs,Fs,As,~,zenith,azimuth] = icohemisphere(2,1);     %calculate light source angles of incidence THIS DETERMINES NR OF ANGLES!!!
    A = length(Fs);                                    %nr of angles of incidence for ray-tracing
    C = length(T(1).Facet);                            %nr of cells (Type1 is always cell)
    sz = size(T(1).RT.RAT);
    if numel(sz)==3, L = sz(3)-1; W = sz(1); elseif numel(sz)==2, L = sz(2)-1; W = 1; else, L = 1; W = 1; end


    et = find(~cellfun(@isempty,{T.Emit}),1); %find light source in geometry
    T(et).Emit(1) = TOOLBOX_input.Scene.module_mounting.NRays;                      %SET NR OF RAYS HERE!!!!
    %initialize sensitivity map (-0.1 means not calculated yet)
    SM_f_i = -0.1*ones(A,C,L,W);                      %dim1: angles, dim2: facets, dim3: layers
    h_f = 99;                                 %figure handle
    if BF, SM_r_i = -0.1*ones(A,C,L,W); h_r = 98; end %if bifacial, also initialize rear side sensitivity map

    for a = 1:A                               %for every angle of incidence
        T(et).Emit(2:3) = [zenith(a),azimuth(a)];   %set azimuth and zenith angle
        SS = Lux58(V,F,T);                    %calculate sensitivity of all facets and layers for that angle
        SM_f_i(a,:,:,:) = SS(1:C,:,:);              %add to 4D sensitivity array (dim1: angle, dim2: facet (cell only), dim3: layer, dim4: wav)
        %plot progress sensitivity map (total absorptance averaged over all cells and all wavelengths)
        if TOOLBOX_input.Scene.plotScene
            h_f = flatplot3(Vs,Fs,mean(mean(SM_f_i(:,:,1,:),4),2),[-0.1,1],parula(512),'Sensitivity [-]',h_f); %plot total absorptance, wavelength average
        end
        if BF
            SM_r_i(a,:,:,:) = SS((C+1):2*C,:,:);
            if TOOLBOX_input.Scene.plotScene
                h_r = flatplot3(Vs,Fs,mean(SM_r_i(:,:,1,30),2),[-0.1,1],parula(512),'Sensitivity [-]',h_r); %plot progress sensitivity map (cell average value)
            end
        end
    end

    if TOOLBOX_input.Scene.module_mounting.avgSensitivity
        SM_f_i = ones(size(SM_f_i)).*mean(SM_f_i,2); %use average sensitivity and keep the size of output constant
        if BF
            SM_r_i = ones(size(SM_r_i)).*mean(SM_r_i,2);
        end
    end
    SM_f{Submod_i} = SM_f_i;
    if BF; SM_r{Submod_i} = SM_r_i; end
    disp('Sensitivity map completed')
end

%scan through all angles to make sensitivity map for every cell front and
%rear, also for every wavelength.
%This requires matrix (layer x wavelength) accumulator to be implemented in
%LUX

MODULE_output.skydome.AZA = [azimuth,zenith,As];    %skydome zenith, azimuth, area for every triangle center
MODULE_output.skydome.Vs = Vs;                      %skydome vertices (for plotting SKYmap)
MODULE_output.skydome.Fs = Fs;                      %skydome facets (for plotting SKYmap)
MODULE_output.SM_f = SM_f;
MODULE_output.wav = CELL_output.CELL_FRONT.wav;     %pass on wavelength information
MODULE_output.N=Ncells;
MODULE_output.A=Acell;
MODULE_output.Amod=Amod;
% Z (L-24 as well)
MODULE_output.ML=ML*1e-2;
MODULE_output.MW=MW*1e-2;
% -
MODULE_output.ModTilt = TOOLBOX_input.Scene.module_mounting.ModTilt;
if BF, MODULE_output.SM_r = SM_r; end
end




