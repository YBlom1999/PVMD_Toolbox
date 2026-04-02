function [V_module] = combineToModuleIV_TripleTandem(V_1,V_2,V_3,Nby,numCells,I,day,night,R,eff_LC,Type,SubMod_ind)
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

if strcmp(Type,'A')
    [V_module1] = CalculateIVmodSingle(V_1,I,R(1),numCells(1),Nby,day,night);

    %Account for LC
    [~,Isc_ind] = min(abs(V_2),[],2);
    Isc2 = I(Isc_ind);
    I_rec2 = max(Isc2'-I,0);
    for i = 1:size(V_2,1)
        V_3(i,:) = interp1(I,V_3(i,:),(I-eff_LC(3)*I_rec2(i,:))',"linear","extrap");
    end


    [V_module2] = CalculateIVmodDouble(V_2,V_3,I,R(2),numCells(2),Nby,day,night);

    V_module=[V_module1;V_module2];

elseif strcmp(Type,'B')
    %Account for LC
    [~,Isc_ind] = min(abs(V_1),[],2);
    Isc1 = I(Isc_ind);
    I_rec1 = max(Isc1'-I,0);
    for i = 1:size(V_2,1)
        V_2(i,:) = interp1(I,V_2(i,:),(I-eff_LC(1)*I_rec1(i,:))',"linear","extrap");
    end


    %Due to the large sizes of the matrix, this steps needs to be defined
    %into different steps.
    if length(day) > 1000
        V_module1 = zeros(length(day)+length(night),length(I));

        N_steps = 10;
        StepSize = (length(day)+length(night))/N_steps;
        N_begin = (0:N_steps-1)*StepSize+1;
        N_end = (1:N_steps)*StepSize;
        for step_i = 1:N_steps
            day_ind = logical((day>=N_begin(step_i)).*(day<=N_end(step_i)));
            day_short = day(day_ind)-N_begin(step_i)+1;
            night_ind = logical((night>=N_begin(step_i)).*(night<=N_end(step_i)));
            night_short = night(night_ind);
            day_pos = find(day_ind);
            V_ind1 = (day_pos(1)-1)*numCells(1)+1;
            V_ind2 = day_pos(end)*numCells(1);
            V_1_short = V_1(V_ind1:V_ind2,:);
            V_2_short = V_2(V_ind1:V_ind2,:);
            V_module1(N_begin(step_i):N_end(step_i),:) = CalculateIVmodDouble(V_1_short,V_2_short,I,R(1),numCells(1),Nby,day_short,night_short);
        end
    else
        V_module1 = CalculateIVmodDouble(V_1,V_2,I,R(1),numCells(1),Nby,day,night);
    end
    [V_module2] = CalculateIVmodSingle(V_3,I,R(2),numCells(2),Nby,day,night);





    V_module=[V_module1;V_module2];

elseif strcmp(Type,'C')

    %Account for LC
    [~,Isc_ind] = min(abs(V_1),[],2);
    Isc1 = I(Isc_ind);
    I_rec1 = max(Isc1'-I,0);
    for i = 1:size(V_2,1)
        V_2(i,:) = interp1(I,V_2(i,:),(I-eff_LC(1)*I_rec1(i,:))',"linear","extrap");
    end
    [~,Isc_ind] = min(abs(V_2),[],2);
    Isc2 = I(Isc_ind);
    I_rec2 = max(Isc2'-I,0);

    for i = 1:size(V_3,1)
        V_3(i,:) = interp1(I,V_3(i,:),(I-eff_LC(2)*I_rec1(i,:)-eff_LC(3)*I_rec2(i,:))',"linear","extrap");


    end
    [V_module] = CalculateIVmodTriple(V_1,V_2,V_3,I,R,numCells,Nby,day,night);
end


end

function [V_mod] = CalculateIVmodSingle(V,I,R,numCells,Nby,day,night)
V=reshape(V-R*I,numCells/Nby,length(day)*Nby,length(I));

%Sum accross all cells conneted to one bypass diode
%V_sum1=zeros(length(day)*Nby,length(I));
V_sum=sum(V,1);

%Reshape so that timeslotts and bypass diodes have their own dimensioninterp1
V_sum=reshape(V_sum,Nby,length(day),length(I));

%module I-V

%Eliminate negative values and sum up all strings
%V__1=zeros(length(day),length(I));
V_tot=sum(max(V_sum,0),1);

%Add 0 curves
V_mod=zeros(length(day)+length(night),length(I));
V_mod(day,:)=V_tot;

end

function [V_mod] = CalculateIVmodDouble(V_1,V_2,I,R,numCells,Nby,day,night)
%Calculate the IV of the module consisting of the middle and bottom cell
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
V_tot=sum(max(V_sum,0),1);

%Add 0 curves
V_mod=zeros(length(day)+length(night),length(I));
V_mod(day,:)=V_tot;
V_mod=max(V_mod,0);

end

function [V_mod] = CalculateIVmodTriple(V_1,V_2,V_3,I,R,numCells,Nby,day,night)
%Calculate the IV of the module consisting of the middle and bottom cell
V_=zeros(3*length(V_1(:,1)),length(I));

%Odd numbers are top cell
top=1:3:3*length(V_1(:,1));

%Even numbers are bottom cell
middle=top+1;
bottom = top+2;

%merge
V_(top,:)=V_1-R*I;
V_(middle,:)=V_2;
V_(bottom,:)=V_3;

% short circuit of each branch with bypass diodes
V_=reshape(V_,3*numCells/Nby,length(day)*Nby,length(I));


%Sum accross all cells conneted to one bypass diode
%V_sum=zeros(length(day)*Nby,length(I));
V_sum=sum(V_);

%Reshape so that timeslotts and bypass diodes have their own dimension
V_sum=reshape(V_sum,Nby,length(day),length(I));

%module I-V
%Eliminate negative values and sum up all strings
%V__=zeros(length(day),length(I));
V_tot=sum(max(V_sum,0),1);

%Add 0 curves
V_mod=zeros(length(day)+length(night),length(I));
V_mod(day,:)=V_tot;
V_mod=max(V_mod,0);

end