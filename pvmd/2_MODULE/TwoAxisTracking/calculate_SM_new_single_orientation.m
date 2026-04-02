function  SM_new_single=calculate_SM_new_single_orientation(CELL_output,azim_vertex, zenith_vertex, wav_WEATHER, ...
    albedo,az_mod,tilt_mod)
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

%Redefine input
wav_MODULE = CELL_output.CELL_FRONT.wav;
az_mod_new=az_mod;
tilt_mod_new=tilt_mod;
alt_mod_new=90-tilt_mod_new;
az_vertex=azim_vertex;
az_vertex(az_vertex >= -180 & az_vertex < 0) = az_vertex (az_vertex >= -180 & az_vertex  < 0) + 360;
zen_vertex=zenith_vertex;
alt_vertex=90-zen_vertex;
angles_cell=CELL_output.CELL_FRONT.aoi;

layers = size(CELL_output.CELL_FRONT.RAT,3);
EQE = CELL_output.CELL_FRONT.RAT;
EQE(:,:,1) = 1-sum(CELL_output.CELL_FRONT.RAT(:,:,[1,end]),3);
EQE(:,:,end) = 1;
wav=CELL_output.CELL_FRONT.wav;

cos_AOI_f=zeros(160,1);
cos_AOI_r=zeros(160,1);
AOI_f=zeros(160,1);
AOI_r=zeros(160,1);

SM_new_single=zeros(160,layers,numel(wav));
for vertex=1:160  % think of changing the am =90-θm so that is θm only dependent

    cos_AOI_f(vertex)=max(cosd(alt_mod_new)*cosd(alt_vertex(vertex))*cosd(az_mod_new-az_vertex(vertex))+sind(alt_mod_new)*sind(alt_vertex(vertex)),0); % 160x1 ignore negative cos_AOI values which means that angle is 90-180o set it to zero deg when negative
    cos_AOI_r(vertex)=max(cosd(alt_mod_new)*cosd(-alt_vertex(vertex))*cosd(az_mod_new-az_vertex(vertex))+sind(alt_mod_new)*sind(-alt_vertex(vertex)),0); % 160x1 rear altitude vertex=-altitude_vertex

    AOI_f(vertex)=acosd(cos_AOI_f(vertex));   % 160x1  returns the angle in the interval [0 180]
    AOI_r(vertex)=acosd(cos_AOI_r(vertex));  % 160x1

    [~,index_f]=min(abs(angles_cell-AOI_f(vertex)));  % Find the closest angle with AOI in cell_block
    [~,index_r]=min(abs(angles_cell-AOI_f(vertex)));

    SM_new_single(vertex,:,:)=squeeze(EQE(:,index_f,:)*cos_AOI_f(vertex)+albedo.*EQE(:,index_r,:)*cos_AOI_r(vertex))';

end


N_rep_before = find(wav_WEATHER == wav_MODULE(1))-1;
N_rep_after = length(wav_WEATHER) - find(wav_WEATHER == wav_MODULE(end));


SM_new_single = cat(3,repelem(SM_new_single(:,:,1),1,1,N_rep_before),SM_new_single,repelem(SM_new_single(:,:,end),1,1,N_rep_after)); % made it from 160x3x46 to 160x3x187 extend along the third=3 dimension 46-- 