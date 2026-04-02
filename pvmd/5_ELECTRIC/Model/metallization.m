function [Op,R] = metallization(CellArea,numBusbars,numFingers,BusbarWdith,FingerWdith)
%%Function to calculate the shading and resistance of the cell metall grid
%%Model taken from: G. Papakonstantinou, “Investigation and Optimization 
%of the Front Metal Contact of Silicon Heterojunction Solar Cells,” no. July, 2014.
%The code was first implement by Abdallah Nour El Din and later reworked and comented by Malte Vogt
%%Input
%CellArea used to calculate the unit cell width 
%numBusbars number of busbars per cell
%numFingers number of fingers per cell
%BusbarWdith; busbar width (m)
%FingerWdith; %finger width (m)
%%Output
%Op - percentatge of the cell, which is shaded
%R - resistance of the cellmetalization
%%TO-Do only tested for square shaped cells not recttanuglar cells

%%Caluclate optical losses


%define variables
Length_unitCell=sqrt(CellArea)/numBusbars;%length of a unit cell convert cm-->m
Fingers_unitCell=numFingers/numBusbars;%numbder fingers per unit cell
Length_finger_unitCell= (Length_unitCell-BusbarWdith)/2; %finger length
half_FingerPitch=Length_unitCell./(2*Fingers_unitCell); %half finger pitch
%optical loss
Op=100*(FingerWdith.*Length_finger_unitCell+half_FingerPitch.*BusbarWdith)./(2*Fingers_unitCell.*half_FingerPitch.^2); 

%Resistance


%define variables
cur_extractPoint=2*(Length_unitCell*100)^2; %current extraction points
S=half_FingerPitch-FingerWdith/2; %half finger spacing
H=Length_unitCell./(2*cur_extractPoint);
RITOsq= 60;
Rfsq= 0.015;
Rbsq= 0.015;
rc= 2e-6;

%Calculate resistance
Lt=sqrt(rc/RITOsq);
RITO=RITOsq*S./(3*Length_finger_unitCell);
Rf=2*Rfsq*Length_finger_unitCell./(3*FingerWdith);
Rc=(RITOsq*Lt./Length_finger_unitCell).*coth(FingerWdith./(2*Lt));
Rb=Rbsq*H./(3*BusbarWdith);
%Resistance
R=(1/numBusbars^2)*(RITO+Rf+Rc)./(4*Fingers_unitCell)+Rb/(2*cur_extractPoint);
