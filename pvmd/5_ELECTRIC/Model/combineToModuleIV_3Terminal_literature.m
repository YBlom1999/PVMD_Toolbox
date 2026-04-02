function [V_module,I_module] = combineToModuleIV_3Terminal_literature(V_1,V_2,Nby,numCells,I,day,night,R,configuration,VM_ratio_m,VM_ratio_n,fig)
%combineToModuleIV_3Terminal_literature Calculates the IV curve of the
%module with 3T cells
%
% This function calculates the module IV curve of 3T tandem cells. It is
% based on the work of McMahon
% [https://doi.org/10.1109/JPHOTOV.2021.3068325]
%
% Parameters
% ----------
% V_1 : double
%   Voltages of the IV curve of top cell
% V_2 : double
%   Voltages of the IV curve of bottom cell cell
% N_by : double
%   Number of bypass diodes (should still be included in the model)
% numCells : double
%   number of cells in the module
% I : double
%   The currents of the IV curve
% day : double
%   Indices of the day
% night : double
%   Indices of the night
% R : double
%   Resistance of the interconnection
% configuration : string
%   Indictator of the cell is series or reverse connected.
% VM_ratio m : double
%   index of how many top cells are connected in parallel
% V_1 : double
%   index of how many bottom cells are connected in parallel.
%
% Returns
% -------
% V_module : double
%   The voltage of the IV curve of the module
% I : double
%   The current of the IV curve of the module
%
% Developed by Y. Blom


%-- The number of active cells are calculated
VM_ratio = VM_ratio_n/VM_ratio_m;
Sum_VM = VM_ratio_m+VM_ratio_n;
if strcmp(configuration,'Series')
    N_active = numCells+1-Sum_VM;
elseif strcmp(configuration,'Reverse')
    N_active = numCells+1-max(VM_ratio_m,VM_ratio_n);
end


%-- The IV curves with nan are replaced with negative numbers
for i=1:length(V_1(:,1))
    k=0;
    for j=1:length(V_1(i,:))
        if isnan(V_1(i,j))
            V_1(i,j)=-20-k*0.001;
            k=k+1;
        end
    end
end
for i=1:length(V_2(:,1))
    k=0;
    for j=1:length(V_2(i,:))
        if isnan(V_2(i,j))
            V_2(i,j)=-20-k*0.001;
            k=k+1;
        end
    end
end

%-- THe voltages are reshaped into 3 dimensional matrices
V_1=reshape(V_1-R(1)*I,numCells,length(day),length(I));
V_2=reshape(V_2-R(2)*I,numCells,length(day),length(I));

I_module = I*Sum_VM;

%-- the IV curve is calculated
if length(day) == 1
    V = 0:0.01:2;
    vstc_top_day = V_1(1,:)/VM_ratio_m;
    vstc_bot_day = V_2(1,:)/VM_ratio_n;
    Istc_top3t = interp1(vstc_top_day,I,V,"linear","extrap");
    Istc_bot3t = interp1(vstc_bot_day,I,V,"linear","extrap");
    Istc_string = (VM_ratio_m*Istc_top3t)+(VM_ratio_n*Istc_bot3t);
    V_module = interp1(Istc_string,N_active*V,I_module);
    
else
        
    V = 0:0.01:2;
    V_module = zeros(length(day)+length(night),length(I));
    for i = 1:length(day)
        vstc_top_day = squeeze(V_1(1,i,:))/VM_ratio_m;
        vstc_bot_day = squeeze(V_2(1,i,:))/VM_ratio_n;
        Istc_top3t = interp1(vstc_top_day,I,V,"linear","extrap");
        Istc_bot3t = interp1(vstc_bot_day,I,V,"linear","extrap");
        Istc_string = (VM_ratio_m*Istc_top3t)+(VM_ratio_n*Istc_bot3t);
        V_module(day(i),:) = interp1(Istc_string,N_active*V,I_module,"linear","extrap");
        
             
        
    end
end
end