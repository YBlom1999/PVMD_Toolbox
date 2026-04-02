function [V_module] = combineToModuleIV_2Terminal(V_1,V_2,Nby,numCells,I,day,night,R,eff_LC)
%%Combines all cell VI curves into one module IV curve
%Input
%V_ -Cell VI curves
%Nby - number of bypass diodes
%numCells - number of cells
%I - current values ocrresponding to V_
%day - time instances with non-zero VI-cures
%night - time instances with  VI-cures=0
%R-Resistance
%Ouput: Module VI curve
%Author: Malte Vogt

%Account for LC
[~,Isc_ind] = min(abs(V_1),[],2);
Isc1 = I(Isc_ind);
I_rec1 = max(Isc1'-I,0);
for i = 1:size(V_2,1)
    V_2(i,:) = interp1(I,V_2(i,:),(I-eff_LC*I_rec1(i,:))',"linear","extrap");
end

%merge vectors with alternating rows
V_=zeros(2*length(V_1(:,1)),length(I));

%Odd numbers are top cell
top=1:2:2*length(V_1(:,1));

%Even numbers are bottom cell
bottom=top+1;

%merge
V_(top,:)=V_1-R*I;
V_(bottom,:)=V_2;

% short circuit of each branch with bypass diodes
V_=reshape(V_,2*numCells/Nby,length(day)*Nby,length(I));


%Sum accross all cells conneted to one bypass diode
%V_sum=zeros(length(day)*Nby,length(I));
V_sum=sum(V_);

%Reshape so that timeslotts and bypass diodes have their own dimension
V_sum=reshape(V_sum,Nby,length(day),length(I));

%module I-V
%Eliminate negative values and sum up all strings
%V__=zeros(length(day),length(I));
V__=sum(max(V_sum,0),1);

%Add 0 curves
V_module=-999*ones(length(day)+length(night),length(I));
V_module(night,:)=zeros(length(night),length(I));
V_module(day,:)=V__;
V_module=max(V_module,0);
end