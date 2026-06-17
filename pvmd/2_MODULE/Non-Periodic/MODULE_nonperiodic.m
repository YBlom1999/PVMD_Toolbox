function [MODULE_output] = MODULE_nonperiodic(TOOLBOX_input, CELL_output)
%WEATHER_periodic Calculates the sensitivity map for non-periodic simulations
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


if TOOLBOX_input.Scene.loadSimulation
    load(CELL_output.TYPE,TOOLBOX_input.Scene.SimulationFile,'MODULE_output')
    TOOLBOX_input_mount = load(CELL_output.TYPE,TOOLBOX_input.Scene.SimulationFile,'TOOLBOX_input');
    if exist('TOOLBOX_input_mount','var')
        TOOLBOX_input.Scene.module_mounting = TOOLBOX_input_mount.TOOLBOX_input.Scene.module_mounting;
        TOOLBOX_input.Scene.N_panels = TOOLBOX_input_mount.TOOLBOX_input.Scene.N_panels;
       
    end
    TOOLBOX_input.Scene.N_panels = length(MODULE_output.SM_f);
else
    N_panels = TOOLBOX_input.Scene.N_panels;
    N_submodules = max(CELL_output.SUBMOD_IND);
    %Calculate SM map
    SM_f = cell(N_panels,N_submodules);
    SM_r = cell(N_panels,N_submodules);

    %Settings for the backward ray tracer
    settings.N_refinement_normal = TOOLBOX_input.Scene.N_refinement_normal;
    settings.N_refinement_reduced = TOOLBOX_input.Scene.N_refinement_reduced;
    settings.reducedRays = 0;

    %Load the environment
    load(fullfile(TOOLBOX_input.Scene.NonPeriodic_Environment),'V','F','opa','rgb','Albedo','Scattering')
    opa_env = opa;
    rgb_env = rgb;
    N_env = size(F,1);
    Albedo_env = Albedo;
    Scattering_env = Scattering;
    V = 100*V;

    %===create an array for the albedo and scattering of every vertex
    wav = CELL_output.CELL_FRONT.wav;
    Albedo = zeros(length(F),length(wav));
    Scattering = zeros(1,length(F));
    Albedo(1:N_env,:) = Albedo_env.*ones(N_env,length(wav));
    Scattering(1:N_env) = Scattering_env;

    %=== Create structure for all the panels
    Panels.Ncells = zeros(1,N_submodules);
    Panels.Acell = zeros(1,N_submodules);
    Panels.Amod = zeros(1,N_submodules);
    Panels = repelem(Panels,N_panels,1);
    ModTilt = zeros(N_panels,1);
    for Panel_i = 1:N_panels
        for Submod_i = 1:N_submodules

            CELL_adjusted = CELL_output;
            Submod_ind = find([1,CELL_output.SUBMOD_IND == Submod_i,1]);
            CELL_adjusted.CELL_FRONT.RAT = CELL_output.CELL_FRONT.RAT(:,:,Submod_ind);
            CELL_adjusted.CELL_FRONT.lay = CELL_output.CELL_FRONT.lay(Submod_ind);
            if isstruct(CELL_output.CELL_REAR)
                CELL_adjusted.CELL_REAR.RAT = CELL_output.CELL_REAR.RAT(:,:,Submod_ind);
                CELL_adjusted.CELL_REAR.lay = CELL_output.CELL_REAR.lay(Submod_ind);
            end

            [V, F,BF,cellCorners,normalSolarCell,Ncells_i,Acell_i,Amod_i,T] = AddPanels_NonPeriodic(TOOLBOX_input,CELL_adjusted,V,F,Panel_i,Submod_i);
            Panels(Panel_i,1).Ncells(Submod_i) = Ncells_i;
            Panels(Panel_i,1).Acell(Submod_i) = Acell_i;
            Panels(Panel_i,1).Amod(Submod_i) = Amod_i;


            %Update Albedo and scattering for panels
            for t = 1:length(T)
                if ~isempty(T(t).Facet)
                    index_start = 2*T(t).Facet(1)-1;
                    index_end = 2*(T(t).Facet(end)-T(t).Facet(1))+index_start+1;
                    index_rgb = (index_start:index_end)+N_env+2*T(end).Facet(end)*(Panel_i-1);
                    if isstruct(T(t).RT)
                        Albedo(index_rgb,:) = mean(T(t).RT.RAT(:,:,1),2)'.*ones(length(index_rgb),length(wav));
                    else
                        Albedo(index_rgb,:) = T(t).RT(1);
                    end
                    if isfield(T(t),'Scat') && ~isempty(T(t).Scat)
                        Scattering(index_rgb) = 1;
                    end
                end
            end

            [SM_f{Panel_i,Submod_i},Vs,Fs,azimuth,zenith,As] = BackwardTracer(V,F,Ncells_i,cellCorners,normalSolarCell,Albedo,Scattering,T(1).RT,settings);
            if TOOLBOX_input.Scene.plotScene
                flatplot3(Vs,Fs,mean(mean(SM_f{Panel_i,Submod_i}(:,:,1,:),4),2),[-0.1,1],parula(512),'Sensitivity [-]',Panel_i+1);
            end
            if BF
                [SM_r{Panel_i,Submod_i},~,~,~,~,~] = BackwardTracer(V,F,Ncells_i,cellCorners,-normalSolarCell,Albedo,Scattering,T(1).RT);
                plot_i = TOOLBOX_input.Scene.N_panels+Panel_i;
                if TOOLBOX_input.Scene.plotScene
                    flatplot3(Vs,Fs,mean(mean(SM_r{Panel_i,Submod_i}(:,:,1,:),4),2),[-0.1,1],parula(512),'Sensitivity [-]',plot_i);
                end
            end
        end
        ModTilt(Panel_i) = TOOLBOX_input.Scene.module_mounting(Panel_i).ModTilt;
    end
    if TOOLBOX_input.Scene.plotScene
        plot_Environment(F,V,N_env,rgb_env,opa_env,TOOLBOX_input,Panels);
    end

    MODULE_output.skydome.AZA = [azimuth,zenith,As];    %skydome zenith, azimuth, area for every triangle center
    MODULE_output.skydome.Vs = Vs;                      %skydome vertices (for plotting SKYmap)
    MODULE_output.skydome.Fs = Fs;                      %skydome facets (for plotting SKYmap)
    MODULE_output.SM_f = SM_f;
    MODULE_output.wav = CELL_output.CELL_FRONT.wav;     %pass on wavelength information
    MODULE_output.Panels = Panels;
    MODULE_output.ModTilt = ModTilt;
    if BF, MODULE_output.SM_r = SM_r; end
end
end


function plot_Environment(F,V,N_env,rgb_env,opa_env,TOOLBOX_input,Panels)
%plot_Environment plots the environment of the simulation
%
% This function makes a figure of the environment and the PV modules
%
% Parameters
% ----------
% F : double
%   The faces of the environment (including the modules)
% V : double
%   The vertices in the environment (including the modules)
% N_env : double
%   The number of faces belonging to the environment (excluding the modules)
% rgb_env : double
%   The color of faces belonging to the environment (excluding the modules)
% opa_env : double
%   The opacity of faces belonging to the environment (excluding the modules)
% TOOLBOX_input : struct
%   Simulation parameters
% Panels : struct
%   Properties of all the PV panels
%
% Returns
% -------
%
% Developed by Youri Blom
nrf = size(F,1);

rgb = zeros(nrf,3);                      %initialize
opa = zeros(nrf,1);                      %initialize
rgb(1:N_env,:) = rgb_env;
opa(1:N_env) = opa_env;

for Panel_i = 1:TOOLBOX_input.Scene.N_panels
    T = Panels(Panel_i).T;
    for t = 1:length(T)
        if ~isempty(T(t).Facet)
            index_start = 2*T(t).Facet(1)-1;
            index_end = 2*(T(t).Facet(end)-T(t).Facet(1))+index_start+1;
            index_rgb = (index_start:index_end)+N_env+2*T(end).Facet(end)*(Panel_i-1);
            if isstruct(T(t).RT)
                rgb(index_rgb,3) = 1;      %give cell blue color
                opa(index_rgb) = 1-mean(mean(T(t).RT.RAT(:,:,end)));
            else
                rgb(index_rgb,:) = T(t).RT(1)*ones(2*numel(T(t).Facet),3);   %others gray-scale
                opa(index_rgb) = 1-T(t).RT(2);
            end
        end
    end
end
p = patch('Faces',F,'Vertices',V);       %plot facets
%set rgb and opacity
set(p,'FaceVertexCData',rgb,'CDataMapping','scaled',...
    'FaceColor','flat','FaceVertexAlphaData',opa,...
    'AlphaDataMapping','none','FaceAlpha','flat')

axis equal off tight                    %xyz aspect ratio equal
view(30,30)                            %standard view from side
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-1000,1000]);
ylim([-1000,1000]);
zlim([0,1000]);

end