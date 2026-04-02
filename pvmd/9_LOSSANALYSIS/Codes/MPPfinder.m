function [Impp, Vmpp] = MPPfinder(Iph,Rs,Rsh,n,I0,T)
%MPPfinder Finds the maximum power point current and voltage for a solar
%cell
%
% This function finds the maximum power point voltage and current based on
% the parameters of the one diode model.
%
% Parameters
% ----------
% Iph : double
%   The photo generated current
% Rs : double
%   The series resistance
% Rsh : double
%   The shunt resistance current
% n : double
%   The ideality factor
% I0 : double
%   The saturation current
% T : double
%   The cell temperature
%
% Returns
% -------
% Impp : double
%   The maximum power point current
% Vmpp : double
%   The maximum power point voltage
%
% Developed by Y. Blom
    V=[-15:-1 0:5e-3:1.5];
    k = 1.380649e-23;
    q = 1.60217662e-19;
    Vth = k*T/q;
    z=(Rs.*I0./(n.*Vth.*(1+Rs./Rsh))).*exp((Rs.*(Iph+I0)+V)./(n.*Vth.*(1+Rs./Rsh)));
    I=(Iph+I0-ones(size(z)).*V./(Rsh))./(1+Rs./Rsh)-lambertw(z).*(n.*Vth)./Rs;
    P = I.*V;
    [~,ind] = max(P,[],2);
    
    %For STC simulations
    if size(I,1) == 1
        Impp = I(ind);
    %For operating conditions
    else
        %Linear indexing is needed
        Index_help = (1 : size(I, 1)) .';
        new_index = sub2ind(size(I), Index_help, ind);
        Impp 