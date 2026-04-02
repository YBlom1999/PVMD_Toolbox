function [opt_az,opt_tilt,A,J,Irr]=cross_method_fun_per_hour(CELL_output,azim_vertex, zenith_vertex, ...
   wav_WEATHER,albedo,Bi,Bf,opt_az_last_hour,opt_tilt_last_hour)
%cross_method_fun_per_hour Finds the optimal orientation by using the cross
%method
%
% This function optimizes the tilt and azimuth of the PV module by using
% the cross method. This method works by using the previous optimum as
% initial guess and iteratively checking the neighbouring orientations.
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
% Bi : double
%   the irradiance from all vertices
% Bf : double
%   the photon flux from all vertices
% opt_az_last_hour : double
%   the optimal azimuth of the last hour
% opt_tilt_last_hour : double
%   the optimal tilt of the last hour
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

%In case want random initial guess uncomment below and in main script 


% function [opt_az_per_hour,opt_tilt_per_hour,opt_abs_per_hour,J_per,J_si]=cross_method_fun_per_hour(CELL_output,azim_vertex, zenith_vertex, ...
%     wav_WEATHER,albedo,layers,Bi,Bf)

%counter_iter=0 %count iterations per_hour to converge
wav=wav_WEATHER; % to be used in the integration, SM function uses as input the wav_Weather imported when calling cross method
lay = length(CELL_output.CELL_FRONT.lay)-2;
azim0=opt_az_last_hour; % make initail guess optimum of last hour
tilt0=opt_tilt_last_hour;

SM=calculate_SM_new_single_orientation(CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
    albedo,azim0,tilt0);

% Calculate initial absorption
current_Irr= trapz(wav,sum(Bi.*squeeze(SM(:,end,:)))'); % 1=total cell layer for Abs, no need to put azim0+1, see below

while true

    % Neighbor points
    neighbor_azim = [azim0+1, azim0-1, azim0, azim0];
    neighbor_tilt = [tilt0, tilt0, tilt0+1, tilt0-1];

    %  Shift azim  at the edge cases
    neighbor_azim(neighbor_azim==360) = 0;   % to make azim goes around
    neighbor_azim(neighbor_azim==-1) = 359;
    neighbor_tilt(neighbor_tilt==-1) = 0;  % to prevent 91 and -1 indices
    neighbor_tilt(neighbor_tilt==91) = 90;

    % Calculate absorption in each neighbor point
    neighbor_Irr = zeros(1,4);

    for i = 1:4
        if (neighbor_tilt(i) >= 0) && (neighbor_tilt(i) <= 90)

            SM=calculate_SM_new_single_orientation(CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
                albedo,neighbor_azim(i),neighbor_tilt(i)); % no need to put azim0+1 because you do not save and call from matrix them you calculate them

            neighbor_Irr(i)=trapz(wav,sum(Bi.*squeeze(SM(:,end,:)))');
        end
    end


    % Find neighbor point with maximum absorption
    [max_neighbor_Irr, max_neighbor_idx] = max(neighbor_Irr);

    % If maximum absorption is found in current point, break loop
    if max_neighbor_Irr <= current_Irr % it exits if none the neighbor points is greater than the current point

        opt_az=azim0;   % azim0 beacause the central point is the optimum before it exits
        opt_tilt=tilt0;
        Irr=current_Irr;

        % Jph current in optimum place
        SM_J=calculate_SM_new_single_orientation(CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
                albedo,opt_az,opt_tilt);

        J = zeros(1,lay);
        A = trapz(wav,sum(Bi.*squeeze(SM_J(:,1,:)))');
        for lay_i = 1:lay
            J(lay_i)=trapz(wav,sum(Bf.*squeeze(SM_J(:,lay_i+1,:)))');
        end

        break
    end
    % Move to neighbor point with maximum absorption
    azim0 = neighbor_azim(max_neighbor_idx);  % every time the initial azim0 is updated
    tilt0 = neighbor_tilt(max_neighbor_idx);
    current_Irr = max_neighbor_Irr;   % and the current_abs is updated with the max
    %counter_iter=counter+1;   %count iterations per_hour to converge
end

end