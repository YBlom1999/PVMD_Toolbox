function [V_module,I_module,ind_Diode] = combineToModuleIV_StandardMod(V_,Nby,numCells,I,day,night,R,ButterFly,TOOLBOX_input)
% combineToModuleIV_StandardMod combines the cell IV curves into a module IV curve for a standard module.
%
% The IV curves of the cells, calclated in the previous step, are
% integrated into a single module IV curve
%
% Parameters
% ----------
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
% R : double
%   The interconnection resistance
% Butterfly : boolean
%   Indicator if the module is has the butterfly topology
% TOOLBOX_input : struct
%   All inputs for the simulation
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

% Account for the series resistance
V_=V_ - I*R;
V_diode = TOOLBOX_input.electric.forwardVoltageDiode;

if ~ButterFly %Non-butterfly topology
    [V__,I_module,ind_Diode] = ConnectCellsStandard(V_,Nby,numCells,I,day,night,V_diode,TOOLBOX_input);
else %Butterfly topology
    [V__,I_module,ind_Diode] = ConnectCellsButterfly(V_,Nby,numCells,I,day,night,V_diode,TOOLBOX_input);
end



%Add 0 curves
V_module=-999*ones(length(day)+length(night),length(I_module));
V_module(night,:)=zeros(length(night),length(I_module));
V_module(day,:)=V__;
V_module=max(V_module,0);

end


function [V__,I_module,ind_Diode] = ConnectCellsStandard(V_,Nby,numCells,I,day,night,V_diode,TOOLBOX_input)
% ConnectCellsStandard combines the cell IV curves into a module IV curve for a standard module.
%
% The IV curves of the cells are interconnected in a series connection,
% with the inclusion of bypass diodes
%
% Parameters
% ----------
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
% V_diode : double
%   The forward bypass voltage of the bypass diodes
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
N_rows = TOOLBOX_input.Scene.module_mounting.CellRows;
N_columns = TOOLBOX_input.Scene.module_mounting.CellColumns;

if N_rows>N_columns
    ind_new = zeros(numCells*length(day),1);
    for t = 1:length(day)
        Start_ind = numCells*(t-1);
        index = [];
        for row_i = 1:N_rows
            for col_i = 1:N_columns
                index_add = row_i+(col_i-1)*N_rows;
                index = [index, index_add];
            end
        end
        ind_new(Start_ind+1:numCells) = Start_ind+index;
    end
    V_(ind_new,:) = V_;

end

V_=reshape(V_,numCells/Nby,length(day)*Nby,length(I));


%Sum accross all cells conneted to one bypass diode
V_sum=sum(V_);

%Reshape so that timeslotts and bypass diodes have their own dimension
V_sum=reshape(V_sum,Nby,length(day),length(I));

%module I-V
%Eliminate negative values and sum up all strings
V__=sum(max(V_sum,-1*V_diode),1);

%indicate when the bypass diodes are active
ind_Diode = zeros(Nby,length(day)+length(night),length(I));
ind_Diode(:,day,:) = V_sum < -1*V_diode;

I_module = I;

end

function [V__,I_module,ind_Diode] = ConnectCellsButterfly(V_,Nby,numCells,I,day,night,V_diode,TOOLBOX_input)
% ConnectCellsStandard combines the cell IV curves into a module IV curve for a butterfly module.
%
% The IV curves of the cells are interconnected in a butterfly topology,
% with the inclusion of bypass diodes
%
% Parameters
% ----------
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
% V_diode : double
%   The forward bypass voltage of the bypass diodes
% TOOLBOX_input : struct
%   All inputs for the simulation
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
N_columns = TOOLBOX_input.Scene.module_mounting.CellColumns;

if N_columns > N_rows
    ind_new = 1:numCells*length(day);
    for t = 1:length(day)
        for i = 1:N_rows/2
            Start_ind = 2*(i-1)*N_colums + numCells*(t-1);
            places1 = Start_ind + (N_colums/2+1:N_colums);
            places2 = Start_ind + (N_colums+1:1.5*N_colums);
            ind_new(places1) = places2;
            ind_new(places2) = places1;
        end
    end
else
    ind_new = zeros(numCells*length(day),1);
    for t = 1:length(day)
        Start_ind = numCells*(t-1);
        index = [];
        for i =1:2 % for the top and bottom strings
            for row_i = 1:N_rows/2
                
                for col_i = 1:N_columns/2
                    index_add1 = (col_i-1)*2*N_rows+(i-1)*N_rows+row_i;
                    index_add2 = index_add1+N_rows/2;
                    index = [index, index_add1,index_add2];
                end
            end
        end
        ind_new(Start_ind+(1:numCells)) = Start_ind+index;
    end
end
V_(ind_new,:) = V_;

V_=reshape(V_,numCells/Nby/2,length(day)*Nby*2,length(I));

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
            [~,ind_un] = unique(squeeze(I_total(k,i,:))');
            if isfield(TOOLBOX_input.electric,'TaylorParam')
                V_intp = V_steps(ind_un);
                [V_intp,I_intp] = CurrentTaylor(squeeze(I_total(k,i,ind_un))',V_intp,TOOLBOX_input.electric.TaylorParam);
                [I_intp,ind_un] = unique(I_intp);
                ind_un(end) = length(V_intp);
                V_intp = V_intp(ind_un);
            else
                V_intp = V_steps(ind_un);
                I_intp = squeeze(I_total(k,i,ind_un))';
            end
            V_strings(k,i,:) = interp1(I_intp,V_intp,I_module,'linear','extrap');
        end
    end
end


%module I-V
%Eliminate negative values and sum up all strings
%V__=zeros(length(day),length(I));
V__=sum(max(V_strings,-1*V_diode),1);

%indicate when the bypass diodes are active
ind_Diode = zeros(Nby,length(day)+length(night),length(I_module));
ind_Diode(:,day,:) = V_strings < -1*V_diode;

end

function [V_out,I_out] = CurrentTaylor(I_in,V,Param)
% CurrentTaylor calculates the IV curve of the substring in a taylor modules
% after the buck converter
%
% The IV curve of the substrings are fed into the buck converter and the
% output IV curve is calculated.
%
% Parameters
% ----------
% I_in: double
%   The current of the string on the input
% V : double
%   The voltage of the string
% Param: struct
%   The parameters of the converter
%
% Returns
% -------
% I_out: double
%   The current of the string on the output
%
% Written by  Y. Blom

R_loss = Param.Rloss;
Iout_max = Param.I0max;
Dmax = Param.Dmax;

V_out = Dmax*V;
P_new = I_in.*V-R_loss*I_in.^2;
[Pmpp,ind_mpp] = max(P_new);
Vmpp = V(ind_mpp);
Impp = I_in(ind_mpp);

P_new(V<Vmpp) = Pmpp - R_loss.*(Pmpp./V_out(V<Vmpp));
I_out = P_new./V_out;
Imax_ind = find(((I_out>Iout_max)+(V_out<=0)),1);
I_out(Imax_ind:end) = Iout_max;



end