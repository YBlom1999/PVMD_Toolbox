function [Shading,R] = metallization(CellArea,numBusbars,numFingers,BusbarWidth,FingerWidth)
% metallization calculates the shading and resistance of the cell metal
% grid. The model is taken from: G. Papakonstantinou, “Investigation and Optimization 
% of the Front Metal Contact of Silicon Heterojunction Solar Cells,” no. July, 2014.
%
% Parameters
% ----------
% CellArea : double
%   Area of the cell
% numBusbars: double
%   The number of busbars
% numFingers: double
%   The number of fingers
% BusbarWidth: double
%   The width of the busbar
% FingerWidth: double
%   The width of the finger
%
% Returns
% -------
% Shading: double
%   The percentage of shading (non active area due to the metalization)
% R: double
%   resistance of the cellmetalization
%
% Developed by Abdallah Nour El Din and later reworked and comented by Malte Vogt

%define variables
Length_unitCell=sqrt(CellArea)/numBusbars;%length of a unit cell convert cm-->m
Fingers_unitCell=numFingers/numBusbars;%numbder fingers per unit cell
Length_finger_unitCell= (Length_unitCell-BusbarWidth)/2; %finger length
half_FingerPitch=Length_unitCell./(2*Fingers_unitCell); %half finger pitch
%optical loss
Shading=100*(FingerWidth.*Length_finger_unitCell+half_FingerPitch.*BusbarWidth)./(2*Fingers_unitCell.*half_FingerPitch.^2); 

%Resistance


%define variables
cur_extractPoint=2*(Length_unitCell*100)^2; %current extraction points
S=half_FingerPitch-FingerWidth/2; %half finger spacing
H=Length_unitCell./(2*cur_extractPoint);
RITOsq= 60;
Rfsq= 0.015;
Rbsq= 0.015;
rc= 2e-6;

%Calculate resistance
Lt=sqrt(rc/RITOsq);
RITO=RITOsq*S./(3*Length_finger_unitCell);
Rf=2*Rfsq*Length_finger_unitCell./(3*FingerWidth);
Rc=(RITOsq*Lt./Length_finger_unitCell).*coth(FingerWidth./(2*Lt));
Rb=Rbsq*H./(3*BusbarWidth);
%Resistance
R=(1/numBusbars^2)*(RITO+Rf+Rc)./(4*Fingers_unitCell)+Rb/(2*cur_extractPoint);
