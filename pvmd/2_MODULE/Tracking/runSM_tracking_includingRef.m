function  SM_stored=runSM_tracking_includingRef(CELL_output,wav_WEATHER, az_mod,tilt_mod,SM_stored,TOOLBOX_input)
%SM_new_single Calculates the sensitivity map of a single orientation for a
%given vertex and zenith
%
% This function calculates the sensitivity map of the PV module
%
% Parameters
% ----------
% CELL_output : struct
%   Output of the CELL block
% azim_vertex : double
%   the azimuth of all vertices
% zenith_vertex : double
%   the zenith of all vertices
% wav_WEATHER : double
%   the wavelengths that are considered by the WEATHER module
% albedo : double
%   the spectral albedo of the ground
% az_mod : double
%   the azimuth of the module
% tilt_mod : double
%   the tilt of the module
%
% Returns
% -------
% SM_new : double
%   The calculated sensitivity map of the PV module
%
% Developed by Orestis Chatzilampos, integrated by Youri Blom

% here wav= moddule_wav=46 and wav_Weather = SMARTS= 187

%SOS change the layers when running single junction cell
% This functions returns the SM_new of my method values for 1 orientation and
% an 160x3x46 --> extented SMARTS --> 160x3x187


%% Get the correct geometry to calculate SM

TOOLBOX_input.Scene.module_mounting.ModTilt = tilt_mod;
TOOLBOX_input.Scene.module_mounting.ModAzimuth = az_mod;

if strcmp(TOOLBOX_input.Scene.trackingType,'HSATS')
    [vertices,faces,T,Bf,Ncells,Acell,Amod,~,~] = moduleGeometryTracking_hsats(TOOLBOX_input,CELL_output,1,1);  %Stems from runBackwardTracer_main.melseif TOOLBOX_input.Scene.TSATS
elseif strcmp(TOOLBOX_input.Scene.trackingType,'TSATS')
    [vertices,faces,T,Bf,Ncells,Acell,Amod,~,~] = moduleGeometryTracking_tsats(TOOLBOX_input,CELL_output,1,1);  %Stems from runBackwardTracer_main.m
elseif strcmp(TOOLBOX_input.Scene.trackingType,'DATS')
    [vertices,faces,T,Bf,Ncells,Acell,Amod,~,~] = moduleGeometryTracking_dats(TOOLBOX_input,CELL_output,1,1);  %Stems from runBackwardTracer_main.m
end

%% Define solar cell
cellCorners = reshape(vertices(1:Ncells*4,:)',3,4,Ncells);

if strcmp(TOOLBOX_input.Scene.trackingType,'TSATS')
    cellNormal = TOOLBOX_input.Scene.module_mounting.cell_normal';
else
    cellNormal = [sind(tilt_mod)*sind(az_mod+180);sind(tilt_mod)*cosd(az_mod+180);cosd(tilt_mod)];
end
cellNormal = repelem(cellNormal,1,Ncells);
 
%% To extend values

EQE_wav=CELL_output.CELL_FRONT.wav;

%% extra info needed for reflection

 %===create an array for the albedo and scattering of every vertex
reflectivity = zeros(length(faces),length(EQE_wav));
Scattering = zeros(length(faces),1);
for F_i = 1:length(faces)
    for T_i = 1:length(T)
        ind = find(T(T_i).Facet == F_i,1);
        if ~isempty(ind)
            if isstruct(T(T_i).RT)
                reflectivity(F_i,:) = mean(T(T_i).RT.RAT(:,:,1),2);
            else
                reflectivity(F_i,:) = T(T_i).RT(1);
            end
            if ~isempty(T(T_i).Scat)
                Scattering(F_i) = 1;
            end
            break
        end
    end
end

settings.N_refinement_normal = TOOLBOX_input.Scene.N_refinement_normal;
settings.N_refinement_reduced = TOOLBOX_input.Scene.N_refinement_reduced;
settings.reducedRays = 1;
SM = BackwardTracer(vertices, faces, Ncells, cellCorners, cellNormal, reflectivity, Scattering, T(1).RT,settings);

if Bf==1
    iCellContainingPlane_tot = TOOLBOX_input.Scene.module_mounting.cell_inplane_rear;
    cellCorners = reshape(vertices((Ncells*4+1):(Ncells*8),:)',3,4,Ncells);
    centroidSolarCell = squeeze(mean(cellCorners,2));

    SM_r = calculateSM_tracking(vertices, faces, Ncells, centroidSolarCell,-1*cellNormal, reflectivity, Scattering, CELL_output.CELL_REAR,iCellContainingPlane_tot);
    SM = SM +SM_r;
end

%use less memory, error in order of 10^-8
SM_new = single(SM); 

if strcmp(TOOLBOX_input.Scene.trackingType,'HSATS') || strcmp(TOOLBOX_input.Scene.trackingType,'TSATS')
    N_rep_before = find(wav_WEATHER == EQE_wav(1))-1;
    N_rep_after = length(wav_WEATHER) - find(wav_WEATHER == EQE_wav(end));

    SM_before = repmat(SM_new(:,:,:,1), 1,1,1,N_rep_before);
    SM_after  = repmat(SM_new(:,:,:,end), 1,1,1,N_rep_after);
    SM_new = cat(4, SM_before, SM_new, SM_after); % made it from 160x72x3x46 to 160x72x3x187 extend along the third=3 dimension 46-- 
end

SM_stored{az_mod+1, tilt_mod+1} = SM_new;
end