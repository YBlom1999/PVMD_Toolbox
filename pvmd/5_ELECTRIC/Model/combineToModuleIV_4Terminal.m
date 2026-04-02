function [V_module] = combineToModuleIV_4Terminal(V_1,V_2,Nby,numCells,I,day,night,R)
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
 %top cell
    %Vector to matrix
    V_1=reshape(V_1-R(1)*I,numCells(1)/Nby,length(day)*Nby,length(I));

    %Sum accross all cells conneted to one bypass diode 
    %V_sum1=zeros(length(day)*Nby,length(I));
    V_sum1=sum(V_1,1);

    %Reshape so that timeslotts and bypass diodes have their own dimensioninterp1
    V_sum1=reshape(V_sum1,Nby,length(day),length(I));

    %module I-V

    %Eliminate negative values and sum up all strings
    %V__1=zeros(length(day),length(I));
    V__1=sum(max(V_sum1,0),1);
    
    %Add 0 curves
    V_module1=-999*ones(length(day)+length(night),length(I));
    V_module1(night,:)=zeros(length(night),length(I));
    V_module1(day,:)=V__1;

    %bottom cell
    V_2=reshape(V_2-R(2)*I,numCells(2)/Nby,length(day)*Nby,length(I));
    %Sum accross all cells conneted to one bypass diode 
    %V_sum2=zeros(length(day)*Nby,length(I));
    V_sum2=sum(V_2,1);
    
    %Reshape so that timeslotts and bypass diodes have their own dimension
    V_sum2=reshape(V_sum2,Nby,length(day),length(I));

    %module I-V

    %Eliminate negative values and sum up all strings
    %V__2=zeros(length(day),length(I));
    V__2=sum(max(V_sum2,0),1);
    %Add 0 curves
    V_module2=-999*ones(length(day)+length(night),length(I));
    V_module2(night,:)=zeros(length(night),length(I));
    V_module2(day,:)=V__2;
    
    V_module=[V_module1;V_module2];
end