function [V_module,I] = combineToModuleIV_3Terminal(V_1,V_2,Nby,numCells,I,day,night,R)
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

    %odd and even cell numbers
    odd_cells=1:2:length(V_2(:,1));
    even_cells=2:2:length(V_2(:,1));
    
    %add V(i) curves of bottom cells in pairs
    V_botPairs = V_2(odd_cells,:) + V_2(even_cells,:);
    
    %Split up top cells in odd and even
    V_topOdd=V_1(odd_cells,:);
    V_topEven=V_1(even_cells,:);
    
   
    %replace nan with V<=-20 (nan results in a error with interpolate)
    for i=1:length(V_topOdd(:,1))
        k=0;
        for j=1:length(V_topOdd(i,:))
            if isnan(V_topOdd(i,j))
                V_topOdd(i,j)=-20-k*0.001;
                k=k+1;
            end
        end
    end    
    for i=1:length(V_topEven(:,1))
        k=0;
        for j=1:length(V_topEven(i,:))
            if isnan(V_topEven(i,j))
                V_topEven(i,j)=-20-k*0.001;
                k=k+1;
            end
        end
    end    
    for i=1:length(V_botPairs(:,1))
        k=0;
        for j=1:length(V_botPairs(i,:))
            if isnan(V_botPairs(i,j))
                V_botPairs(i,j)=-20-k*0.001;
                k=k+1;
            end
        end
    end    

    
    %Prepare inversion to I(V_steps) for paralllel connection
    V_max=max([max(V_topOdd),max(V_topEven),max(V_botPairs)])*1.05;
    V_steps=0:0.001:round(V_max,3);
    I_topOdd=zeros(length(V_topOdd(:,1)),length(V_steps));
    I_topEven=zeros(length(V_topEven(:,1)),length(V_steps));
    I_botPairs=zeros(length(V_botPairs(:,1)),length(V_steps));
    
    %Inversion, Interpolate only woks with vectors
    for i=1:length(V_topOdd(:,1)) 
        if V_topOdd(i,1)>0 &&  V_topEven(i,1)>0 && V_botPairs(i,1)>0%zeros cause problems with interpolate
        I_topOdd(i,:)=interp1(V_topOdd(i,:),I,V_steps);
        I_topEven(i,:)=interp1(V_topEven(i,:),I,V_steps);
        I_botPairs(i,:)=interp1(V_botPairs(i,:),I,V_steps);
        end
   end
    
    %Sum up parrallel unit element (three parallel strings )
    I_unitCell=I_topOdd+I_topEven+I_botPairs;
    
    
    %Invert back to V(I)
    %Resolution of the IV curve in A
    Imax=max(max(I_unitCell))*1.05;
    I=0:5e-3:Imax;    
    %replace nan with I<=-20 (nan results in a error with interpolate)
    for i=1:length(I_unitCell(:,1))
        k=0;
        for j=1:length(I_unitCell(i,:))
            if isnan(I_unitCell(i,j))
                I_unitCell(i,j)=-20-k*0.001;
                k=k+1;
            end
        end
    end    
    
    V_unitCell=zeros(length(I_unitCell(:,1)),length(I));
    for i=1:length(I_unitCell(:,1)) 
        if I_unitCell(i,1)>0 %zeros cause problems with interpolate
            V_unitCell(i,:)=interp1(I_unitCell(i,:),V_steps,I);
        end
    end
    
    % short circuit of each branch with bypass diodes
    V_=reshape(V_unitCell,0.5*numCells/Nby,length(day)*Nby,length(I));


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
    V_module=V_module-R*I/3;
    V_module=max(V_module,0);
end