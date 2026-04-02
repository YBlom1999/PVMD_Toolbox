function Vopt = Vopt_calculator(wav,photon_spec,E_g,T_cell,Angle_emit)
%Vopt_calculator Calculates the optimal voltage for the solar cell
%
% This function calculates the optimal voltage for a solar cell. The
% optimal voltage is the voltage that maximises the power of an ideal soalr
% cell
%
% Parameters
% ----------
% wav : double
%   The wavelenghts
% photon_spec : double
%   The spectral photon flux.
% E_g : double
%   The bandgap energy
% T_cell : double
%   The temperature of the cell
% Angle_emit : double
%   The solid angle of emission
%
% Returns
% -------
% Vopt : double
%   The optimal voltage for this incoming irradiance and temperature.
%
% Developed by Y. Blom
h = 6.62607004e-34;
q = 1.60217662e-19;
c = 299792458;
k = 1.380649e-23;
lambda_g = h*c/(E_g*q);
N_lambda = find(wav > h*c/(q*E_g),1)-1;
E = h*c./wav;
N = 100;
Vopt = linspace(0,E_g,N);
J_abs = zeros(size(photon_spec,1),N);
for i = 1:N
    photon_BB_cell = (2*Angle_emit*c./(wav.^4)*1./(exp((E-q*Vopt(i))*(1./(k*T_cell))')-1))';
    J_in = (photon_spec(:,1:N_lambda)-photon_BB_cell(:,1:N_lambda))*q;
    J_abs(:,i) = trapz(wav(1:N_lambda),J_in')';
end
P = J_abs.*Vopt;
[~,ind_m