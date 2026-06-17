function [V_module,I_module,ind_Diode] = connectCells2ModuleButterfly(TOOLBOX_input,V_,Nby,numCells,I,day,night)
% ConnectCellsStandard combines the cell IV curves into a module IV curve for a butterfly module.
%
% The IV curves of the cells are interconnected in a butterfly topology,
% with the inclusion of bypass diodes
%
% Parameters
% ----------
% TOOLBOX_input : struct
%   All inputs for the simulation
% V_ : double
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
% V_module : double
%   The voltage of the module IV curve
% I_module : double
%   The current of the module IV curve
% ind_Diode : boolean
%   An indicator of when the bypass diodes are active
%
% Written by M. R. Vogt
% Adjusted by Y. Blom
% Reshape for each branch (per bypass diodes and of two sides)

%Make new index for butterfly topology
N_rows = TOOLBOX_input.Scene.module_mounting.CellRows;
N_colums = TOOLBOX_input.Scene.module_mounting.CellColumns;
ind_new = 1:numCells*length(day);
for t = 1:length(day)
    for i = 1:N_colums/2
        Start_ind = 2*(i-1)*N_rows + numCells*(t-1); 
        places1 = Start_ind + (N_rows/2+1:N_rows);
        places2 = Start_ind + (N_rows+1:1.5*N_rows);
        ind_new(places1) = places2;
        ind_new(places2) = places1;
    end
end


V_=reshape(V_(ind_new,:),numCells/Nby/2,length(day)*Nby*2,length(I));

%Sum accross all cells conneted to one string
V_sum=sum(V_);

%Reshape of the two blocks and diodes for each timeslot
V_sum=reshape(V_sum,Nby*2,length(day),length(I));

%Make a step arrays for the voltages and currents
V_steps = linspace(-10,max(max(max(V_sum))),1000);
I_strings = zeros(Nby*2,length(day),length(V_steps));
I_total = zeros(Nby,length(day),length(V_steps));
V_strings = zeros(Nby,length(day),length(I));

%The calculation is done for each hour of the day seperately
for i = 1:length(day)
    for k = 1:Nby*2
        if max(V_sum(k,i,:)) > 0
            %Switch from current steps to voltage steps (from V(I) to I(V)) such that currents can be added
            I_strings(k,i,:) = interp1(squeeze(V_sum(k,i,:))',I,V_steps,'linear','extrap');
        end
    end
    % The current for bypass diode is the sum of the two strings
    I_intermediate = reshape(I_strings(:,i,:),2,Nby,length(V_steps));
    I_total(:,i,:) = sum(I_intermediate);

    % A new current step array is created
    I_module = 2*I;
    for k = 1:Nby
        if max(I_total(k,i,:)) > 0
            %Switch from voltage steps to current steps (from I(V) to V(I)) such that voltages can be added
            V_strings(k,i,:) = interp1(squeeze(I_total(k,i,:))',V_steps,I_module,'linear','extrap');
        end
    end
end


%module I-V
%Eliminate negative values and sum up all strings
%V__=zeros(length(day),length(I));
V__=sum(max(V_strings,0),1);

%indicate when the bypass diodes are active
ind_Diode = zeros(Nby,length(day)+length(night),length(I_module));
ind_Diode(:,day,:) = V_strings < 0;

%Add 0 curves
V_module=-999*ones(length(day)+length(night),length(I_module));
V_module(night,:)=zeros(length(night),length(I_module));
V_module(day,:)=V__;
V_module=max(V_module,0);

end