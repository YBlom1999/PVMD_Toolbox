function [V_module,I_module,ind_Diode] = connectCells2ModuleSeries(Vcell,Nby,numCells,I,day,night)
% connectCellsSeries combines the cell IV curves into a module IV curve for a standard module.
%
% The IV curves of the cells are interconnected in a series connection,
% with the inclusion of bypass diodes
%
% Parameters
% ----------
% Vcell : double
%   The voltage of all cell IV curves
% Nby : double
%   The number of bypass diodes
% numCells : double
%   The number of cells
% I : double
%   The current of all cell IV curves
% day : double
%   The indices of which times are day time
% night : double
%   The indices of which times are night time
%
% Returns
% -------
% V__ : double
%   The voltage of the module IV curve
% I_module : double
%   The current of the module IV curve
% ind_Diode : boolean
%   An indicator of when the bypass diodes are active
%
% Written by M. R. Vogt
% Adjusted by Y. Blom

Vcell=reshape(Vcell,numCells/Nby,length(day)*Nby,length(I));

%Sum accross all cells conneted to one bypass diode
V_sum=sum(Vcell);

%Reshape so that timeslotts and bypass diodes have their own dimension
V_sum=reshape(V_sum,Nby,length(day),length(I));

%module I-V
%Eliminate negative values and sum up all strings
V__=sum(max(V_sum,0),1);

%indicate when the bypass diodes are active
ind_Diode = zeros(Nby,length(day)+length(night),length(I));
ind_Diode(:,day,:) = V_sum < 0;

I_module = I;


%Add 0 curves
V_module=-999*ones(length(day)+length(night),length(I_module));
V_module(night,:)=zeros(length(night),length(I_module));
V_module(day,:)=V__;
V_module=max(V_module,0);

end