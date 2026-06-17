function [opt_az,opt_tilt,A,J,Irr,SM_stored]=CrossMethod(CELL_output,wav,Bi,Bf,opt_az_last_hour,opt_tilt_last_hour,SM_stored,TOOLBOX_input,Ncells)
% cross_method_fun_per_hour Finds the optimal orientation by using the cross
% method
%
% This function optimizes the tilt and azimuth of the PV module by using
% the cross method. This method works by using the previous optimum as
% initial guess and iteratively checking the neighbouring orientations.
%
% Parameters
% ----------
% CELL_output : struct
%   Output of the CELL block
% wav : double
%   the wavelengths that are considered by the WEATHER module
% azim_vertex : double
%   the azimuth of all vertices
% zenith_vertex : double
%   the zenith of all vertices
% albedo : double
%   the spectral albedo of the ground
% Bi : double
%   the irradiance from all vertices
% Bf : double
%   the photon flux from all vertices
% opt_az_last_hour : double
%   the optimal azimuth of the last hour
% opt_tilt_last_hour : double
%   the optimal tilt of the last hour
% SM_stored : cell
%   all stored sensitivity maps
% TOOLBOX_input: struct
%   A structure with all inputs
% N_cells: double
%   Number of cells
%
% Returns
% -------
% opt_az : double
%   The optimal azimuth found
% opt_tilt : double
%   The optimal tilt found
% A: double
%   The absorbed irradiance for the optimal orientation
% J: double
%   The absorbed photon flux in each cell for the optimal orientation
% Irr: double
%   The received irradiance for the optimal orientation
%
% Developed by Orestis Chatzilampos, integrated by Youri Blom

lay = length(CELL_output.CELL_FRONT.lay)-2;
EQE_wav=CELL_output.CELL_FRONT.wav;

mod_tilt = opt_tilt_last_hour;   %This is the tilt of the module!!
mod_az = opt_az_last_hour;

if strcmp(TOOLBOX_input.Scene.trackingType,'DATS')
    mod_az = floor(mod_az/20) * 20;     %this jumps to nearest multiple of 20
    mod_tilt = floor(mod_tilt/20) * 20;     %will significantly increase computational time

    n_max_iter_fine = 2;    %maximum iterations with finer search
else
    n_max_iter_fine = 1;
end

if isempty(SM_stored{mod_az+1, mod_tilt+1})
    SM_stored=runSM_tracking_includingRef(CELL_output,wav,mod_az,mod_tilt,SM_stored,TOOLBOX_input); % no need to put azim0+1 because you do not save and call from matrix them you calculate them
end

SM = SM_stored{mod_az+1,mod_tilt+1};

if strcmp(TOOLBOX_input.Scene.trackingType,'DATS')                % in DATS only a smaller portion of SM is stored to reduce required memory
    N_rep_before = find(wav == EQE_wav(1))-1;
    N_rep_after = length(wav) - find(wav == EQE_wav(end));
    SM_before = repmat(SM(:,:,:,1), 1,1,1,N_rep_before);
    SM_after  = repmat(SM(:,:,:,end), 1,1,1,N_rep_after);
    SM = cat(4, SM_before, SM, SM_after); % made it from 160x72x3x46 to 160x72x3x187 extend along the third=3 dimension 46--
end

current_Irr = zeros(1,Ncells);

for iCell = 1:Ncells
    current_Irr(1,iCell) = trapz(wav,sum(Bi.*squeeze(SM(:,iCell,end,:)))');
end

harm_mean = numel(current_Irr) / sum(1 ./ current_Irr);
current_max = harm_mean;

fine_search = false;    %to increase search speed for dats
n_iter_fine = 0;

while n_iter_fine <= n_max_iter_fine    % stop generating neigbour orientations

    % Neighbor points

    if strcmp(TOOLBOX_input.Scene.trackingType,'HSATS') || strcmp(TOOLBOX_input.Scene.trackingType,'TSATS')

        if fine_search
            neighbor_tilt = [mod_tilt+6,mod_tilt+5,mod_tilt+4,mod_tilt+3,mod_tilt+2,mod_tilt+1, mod_tilt-1,mod_tilt-2, mod_tilt-3,mod_tilt-4,mod_tilt-5, mod_tilt-6];
        else
            step = 10;
            neighbor_tilt = [mod_tilt+3*step,mod_tilt+2*step,mod_tilt+step,mod_tilt-step,mod_tilt-2*step,mod_tilt-3*step];
        end

        neighbor_azim = [opt_az_last_hour];
        neighbor_azim = repelem(neighbor_azim,1,numel(neighbor_tilt));

    elseif strcmp(TOOLBOX_input.Scene.trackingType,'DATS')

        if fine_search
            step = floor(step/2);
            neighbor_tilt = [mod_tilt+step,mod_tilt+step,mod_tilt-step,mod_tilt-step];
            neighbor_azim = [mod_az+step, mod_az-step,mod_az+step, mod_az-step];
        else
            step = 20;
            neighbor_tilt = [mod_tilt+step,mod_tilt+step,mod_tilt-step,mod_tilt-step,mod_tilt,mod_tilt, mod_tilt+step,mod_tilt-step];
            neighbor_azim = [mod_az+step, mod_az-step,mod_az+step, mod_az-step, mod_az+step, mod_az-step, mod_az, mod_az];
        end

    end

    neighbor_tilt = neighbor_tilt(neighbor_tilt<=60);
    neighbor_azim = neighbor_azim(neighbor_tilt<=60);

    neighbor_azim(neighbor_tilt<0) = neighbor_azim(neighbor_tilt<0)+180;
    neighbor_tilt(neighbor_tilt<0) =  neighbor_tilt(neighbor_tilt<0).*(-1);

    neighbor_azim(neighbor_azim>359) = neighbor_azim(neighbor_azim>359)-360;
    neighbor_azim(neighbor_azim<0) = neighbor_azim(neighbor_azim<0)+360;


    if strcmp(TOOLBOX_input.Scene.trackingType,'DATS')
        valid_idx = neighbor_azim <= 120 | neighbor_azim >= 240;       %% Remove neighbor azimuths that fall inside the forbidden interval (120° to 240°)
        neighbor_azim = neighbor_azim(valid_idx);
        neighbor_tilt = neighbor_tilt(valid_idx);
    end

    % This is what you try to max, can be irr, absorption, J
    neighbor_max = zeros(1,numel(neighbor_tilt));


    for i = 1:numel(neighbor_tilt)
        if isempty(SM_stored{neighbor_azim(i)+1, neighbor_tilt(i)+1})

            SM_stored= runSM_tracking_includingRef(CELL_output,wav,neighbor_azim(i),neighbor_tilt(i),SM_stored,TOOLBOX_input); % no need to put azim0+1 because you do not save and call from matrix them you calculate them

        end

        SM = SM_stored{neighbor_azim(i)+1,neighbor_tilt(i)+1};

        if strcmp(TOOLBOX_input.Scene.trackingType,'DATS')                 % in DATS only a smaller portion of SM is stored to reduce required memory
            SM_before = repmat(SM(:,:,:,1), 1,1,1,N_rep_before);
            SM_after  = repmat(SM(:,:,:,end), 1,1,1,N_rep_after);
            SM = cat(4, SM_before, SM, SM_after); % made it from 160x72x3x46 to 160x72x3x187 extend along the third=3 dimension 46--
        end

        Irr = zeros(1,Ncells);
        for iCell = 1:Ncells
            Irr(1,iCell) = trapz(wav,sum(Bi.*squeeze(SM(:,iCell,end,:)))');
        end
        harm_mean = numel(Irr) / sum(1 ./ Irr);
        neighbor_max(i) = harm_mean;


    end

    % Find neighbor point with maximum absorption
    [max_neighbor_Irr, max_neighbor_idx] = max(neighbor_max);


    if max_neighbor_Irr <= current_max
        fine_search = true;     %now we will go look for optimal orientation with finer resolution

        n_iter_fine = n_iter_fine+1;
    else
        % Move to neighbor point with maximum absorption
        mod_tilt = neighbor_tilt(max_neighbor_idx);
        mod_az = neighbor_azim(max_neighbor_idx);
        current_max = max_neighbor_Irr;   % and the current_abs is updated with the max

        if fine_search
            n_iter_fine = n_iter_fine+1;
        end
    end
end




opt_az=mod_az;   % azim0 beacause the central point is the optimum before it exits
opt_tilt=mod_tilt;
Irr = zeros(1, TOOLBOX_input.Scene.module_mounting.Ncells);
A = zeros(1, TOOLBOX_input.Scene.module_mounting.Ncells);
J =zeros(TOOLBOX_input.Scene.module_mounting.Ncells,lay);


SM_J = SM_stored{opt_az+1,opt_tilt+1};

if strcmp(TOOLBOX_input.Scene.trackingType,'DATS')                 % in DATS only a smaller portion of SM is stored to reduce required memory
    SM_before = repmat(SM_J(:,:,:,1), 1,1,1,N_rep_before);
    SM_after  = repmat(SM_J(:,:,:,end), 1,1,1,N_rep_after);
    SM_J = cat(4, SM_before, SM_J, SM_after); % made it from 160x72x3x46 to 160x72x3x187 extend along the third=3 dimension 46--
end

for lay_i = 1:lay
    for iCell = 1:Ncells
        Irr(1,iCell) = trapz(wav,sum(Bi.*squeeze(SM_J(:,iCell,end,:)))');
        A(1,iCell) = trapz(wav,sum(Bi.*squeeze(SM_J(:,iCell,1,:)))');
        J(iCell,lay_i)=trapz(wav,sum(Bf.*squeeze(SM_J(:,iCell,lay_i+1,:)))');
    end

end
end