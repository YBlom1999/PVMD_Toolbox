function [V_,day,night] = InterpolateCell(Acell,pvModIrr,pvModTemp,TOOLBOX_input,numCells,I_,JphSTC)
% Generate IV curves by solving the double diode model equatinon
% [viDataset,gtDataset,cellParam] = generateSCIVcurvesDDM()
% [viDataset,gtDataset,cellParam] = generateSCIVcurvesDDM(...)
% Order of optional parameters: cellArea (cm2) ; Rs (ohm) ; Rp (ohm) ; Is1 (A) ;
% n1 ;Is2 (A) ; n2 ; Iph0 (A); Eg (eV); TIPH1(1/degC) ; TRP1 ; TRS1; TXIS1;
% TXIS2; minIrr (W/m2); stepIrr (W/m2); maxIrr (W/m2); minTemp (degC); 
% stepTemp (degC); maxTemp (degC); filename
% Check parameters in https://uk.mathworks.com/help/physmod/sps/ref/solarcell.html


day=find(sum(pvModIrr,2)>0);
night=find(sum(pvModIrr,2)==0); 

Iph0 = JphSTC*Acell*(1-TOOLBOX_input.electric.shading/100);
Rs0 = TOOLBOX_input.electric.IVparameters(1);
Rp0 = TOOLBOX_input.electric.IVparameters(2);
Is10 = TOOLBOX_input.electric.IVparameters(3);
n1 = TOOLBOX_input.electric.IVparameters(4);
Is20 = TOOLBOX_input.electric.IVparameters(5);
n2 = TOOLBOX_input.electric.IVparameters(6);
Eg = TOOLBOX_input.electric.IVparameters(7);
T_Iph = TOOLBOX_input.electric.IVparameters(8);
Tr_p1 = TOOLBOX_input.electric.IVparameters(9);
Tr_s1 = TOOLBOX_input.electric.IVparameters(10);
TX_Is1 = TOOLBOX_input.electric.IVparameters(11);
TX_Is2 = TOOLBOX_input.electric.IVparameters(12);

minIrr = TOOLBOX_input.electric.IVconditions(1);
stepIrr = TOOLBOX_input.electric.IVconditions(2);
maxIrr = TOOLBOX_input.electric.IVconditions(3);
minTemp = TOOLBOX_input.electric.IVconditions(4);
stepTemp = TOOLBOX_input.electric.IVconditions(5);
maxTemp = TOOLBOX_input.electric.IVconditions(6);

irrs = minIrr:stepIrr:maxIrr;
temps = minTemp:stepTemp:maxTemp;

GTs = [repelem(irrs,length(temps))' repmat(temps,1,length(irrs))'];
V = [0:0.025:0.3,0.325:0.005:0.42,0.421:0.001:0.8];%Adjust if needed (only if you understand what you are doing)
G0 = 1000;
T0 = 25+273.15;
%
nIVs = size(GTs,1);
viDataset = repmat(single(V'),1,2,nIVs);
vDataset = zeros(nIVs,length(I_));
viDataset(:,2,:) = NaN; %current values are made NaN before trying to solve the 2-diode model
for ix=1:nIVs
    G = GTs(ix,1);
    T = GTs(ix,2)+273.15;
    Vt = 1.38064852e-23/1.60217662e-19*T;
    Iph = Iph0*G/G0*(1+T_Iph*(T-T0));
    Rs = Rs0*(T/T0)^Tr_s1;
    Rp = Rp0*(T/T0)^Tr_p1;
    Is1 = Is10*(T/T0)^(TX_Is1)*exp(Eg*(T/T0-1)/(n1*Vt));
    Is2 = Is20*(T/T0)^(TX_Is2)*exp(Eg*(T/T0-1)/(n2*Vt));
    fun = @(I) Iph-Is1*(exp((V+I*Rs)/(n1*Vt))-1)-Is2*(exp((V+I*Rs)/(n2*Vt))-1)-(V+I*Rs)/Rp-I;
    current = fsolve(fun,zeros(1,length(V)),optimset('Display','none'));
    viDataset(:,2,ix) = single(current);
    vDataset(ix,:) = interp1(current,V,I_,'linear','extrap');
end

%% Saving
GTs(:,2) = GTs(:,2)+273.15;
gtDataset = int16(GTs);

% The irradiance and temperature of the solar cells must be integer values.
% If the input is not integer, then the irradiance and temperatures are
% converted to integers.
if ~isa(pvModIrr,'int16')
    aux = int16(pvModIrr);
    aux1 = pvModIrr-single(aux);
    aux1 = max(abs(aux1(:)));
    pvModIrr = aux;
    clear aux1 aux;
end

if ~isa(pvModTemp,'int16')
    aux = int16(pvModTemp);
    aux1 = pvModTemp-single(aux);
    aux1 = max(abs(aux1(:)));
    pvModTemp = aux;
    clear aux1 aux;
end



irr = int16(unique(gtDataset(:,1))); %irradiance values in cell's dataset
nIrr = length(irr);
temp = int16(unique(gtDataset(:,2))); %temperature values in cell's dataset
nTemp = length(temp);

%Check if gtDataset is correctly equispaced
if sum(gtDataset(:,1) ~= repelem(irr,nTemp))
    error('gtDataset: Some irradiance values might be missing.');
elseif sum(gtDataset(:,2) ~= repmat(temp,[nIrr,1]))
    error('gtDataset: Some temperature values might be missing.');
end
if sum(diff(diff(irr))) || sum(diff(diff(temp)))
    error('gtDataset: Temperature and irradiance vectors spacing must be equispaced.');
end

%For each module irradiance and tempearture value, we look for index to
%the closest values in the temperature and irradiance vectors of the cell's
%I-V dataset. To avoid loops and since the module irradiance and
%temperature can be a 3D matrix, we use 4D matrices to find the indices.
irr4dim(1,1,1,:) = irr;
temp4dim(1,1,1,:) = temp;
pvModIrr = reshape(pvModIrr(day,:)',numCells*length(day),1);
pvModTemp = reshape(pvModTemp(day,:)',numCells*length(day),1);
[~,iIrr] = min(abs(pvModIrr-irr4dim),[],4);%indices to find closest value in the irradiance vector (gs)
iIrr = uint16(iIrr);
[~,iTemp] = min(abs(pvModTemp-temp4dim),[],4);%indices to find closest value in the temperature vector (vs)
iTemp = uint16(iTemp);

%Calculate the indeces for the cell's IV dataset combining the indices for
%the irradiance and temperature vectors. The data type is uint32 or uint16
%depending on the size of the I-V dataset.
if size(gtDataset,1) <= 2^16-1
    iCellDataset = (iIrr(:)-1)*nTemp+iTemp(:);
elseif size(gtDataset,1) <= 2^32-1
    iCellDataset = (uint32(iIrr(:))-1)*nTemp+uint32(iTemp(:));
else
    error('The I-V dataset is too large.');
end

V_=zeros(length(day)*numCells,length(I_));
for i=1:nIVs
    index=find(iCellDataset == i);
    V_(index,:)=repmat(vDataset(i,:),length(index),1);
end
end



