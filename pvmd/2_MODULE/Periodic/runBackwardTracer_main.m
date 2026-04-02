function [MODULE_output, TOOLBOX_input] = runBackwardTracer_main(TOOLBOX_input, CELL_output)
%runBackwardTracer_main Runs the backward raytracing.
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
% Developed by Youri Blom

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

    [V,F,T,BF,Ncells(Submod_i),Acell(Submod_i),Amod(Submod_i),ML,MW,TOOLBOX_input] = moduleGeometry(TOOLBOX_input,CELL_adjusted,1,Submod_i);
    Ncells_i = Ncells(Submod_i);
    %===obtain parameters needed for the backward ray tracer
    wav = CELL_output.CELL_FRONT.wav;

    % [SM_f, SM_r, skydome] = Calculate_SM(V{1},F{1},T{1},BF,Ncells(1), wav,TOOLBOX_input);


    cellCorners = reshape(V(1:Ncells_i*4,:)',3,4,Ncells_i);
    TL = TOOLBOX_input.Scene.module_mounting.ModTilt;
    AZ = TOOLBOX_input.Scene.module_mounting.ModAzimuth;
    normalSolarCell = repelem([sind(TL)*sind(AZ+180);sind(TL)*cosd(AZ+180);cosd(TL)],1,Ncells_i);

    %===create an array for the albedo and scattering of every vertex
    Albedo = zeros(length(F),length(wav));
    Scattering = zeros(length(F),1);
    for F_i = 1:length(F)
        for T_i = 1:length(T)
            ind = find(T(T_i).Facet == F_i,1);
            if ~isempty(ind)
                if isstruct(T(T_i).RT)
                    Albedo(F_i,:) = mean(T(T_i).RT.RAT(:,:,1),2);
                else
                    Albedo(F_i,:) = T(T_i).RT(1);
                end
                if ~isempty(T(T_i).Scat)
                    Scattering(F_i) = 1;
                end
                break
            end
        end
    end
    %=== Run backward raytracer
    if TOOLBOX_input.Scene.module_mounting.avgSensitivity
        cellCorners = mean(cellCorners,3);
        normalSolarCell = mean(normalSolarCell,2);
        [SM_mean,Vs,Fs,azimuth,zenith,As] = BackwardTracer(V,F,1,cellCorners,normalSolarCell,Albedo,Scattering,T(1).RT);
        SM_f{Submod_i} = ones(size(SM_mean,1),Ncells_i,size(SM_mean,3),size(SM_mean,4)).*SM_mean;
    else
        [SM_f{Submod_i},Vs,Fs,azimuth,zenith,As] = BackwardTracer(V,F,Ncells_i,cellCorners,normalSolarCell,Albedo,Scattering,T(1).RT);
    end
    flatplot3(Vs,Fs,mean(mean(SM_f{Submod_i}(:,:,1,:),4),2),1);
    if BF
        cellCorners = reshape(V((Ncells_i*4+1):(Ncells_i*8),:)',3,4,Ncells_i);
        if TOOLBOX_input.Scene.module_mounting.avgSensitivity
            cellCorners = mean(cellCorners,3);
            [SM_mean,Vs,Fs,azimuth,zenith,As] = BackwardTracer(V,F,1,cellCorners,-1*normalSolarCell,Albedo,Scattering,T(2).RT);
            SM_r{Submod_i} = ones(size(SM_mean,1),Ncells_i,size(SM_mean,3),size(SM_mean,4)).*SM_mean;
        else
            [SM_r{Submod_i},Vs,Fs,azimuth,zenith,As] = BackwardTracer(V,F,Ncells,cellCorners,-1*normalSolarCell,Albedo,Scattering,T(2).RT);
            
        end
        flatplot3(Vs,Fs,mean(mean(SM_r{Submod_i}(:,:,1,:),4),2),2);
    end

end

MODULE_output.skydome.AZA = [azimuth,zenith,As];    %skydome zenith, azimuth, area for every triangle center
MODULE_output.skydome.Vs = Vs;                      %skydome vertices (for plotting SKYmap)
MODULE_output.skydome.Fs = Fs;                      %skydome facets (for plotting SKYmap)
MODULE_output.SM_f = SM_f;
MODULE_output.wav = wav;     %pass on wavelength information
MODULE_output.N=Ncells;
MODULE_output.A=Acell;
MODULE_output.Amod=Amod;
MODULE_output.ML=ML*1e-2;
MODULE_output.MW=MW*1e-2;
MODULE_output.ModTilt = TOOLBOX_input.Scene.module_mounting.ModTilt;
if BF, MODULE_output.SM_r = SM_r; end
end

