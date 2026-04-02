function [moduleTemp,thermalNet,nodeTypeMat] = cellTemp_incropera2D(pvModule,weather,varargin)
% Sintax:
% [moduleTemp,thermalNet,nodeTypeMat] = incropera2D(pvModule,weather)
% [moduleTemp,thermalNet,nodeTypeMat] = incropera2D(pvModule,weather,key1,val1,....)
% pvModule is a 1D array of PvLayer objects. The first object in the array
% is the top layer of the PV module (the illuminated side).
% weather is a struct with the following fields: localTime (1D datetime
% array); tAmb (1D double array, Celcius); tGnd (1D double array, Celsius);
% tSky (1D double array, Celsius); wind (1D double array, m/s); qGen (2D
% double array, #rows = length(pvModule), #columns = length(weather.localTime), W/m3).
%
% For more details on this thermal model: https://doi.org/10.1002/pip.3504
%
% Valid optional keys and values:
% key: 'rear convection' >> value: double (= 1 by default). It indicates the
% ratio between the wind speed in the rear and front sides. For example 0.5
% means that the wind speed on the rear side is 50% of that on the front
% side. If the value is 0 it is assumed that there is no forced nor free
% convection on the rear side. If only free convection on the back side is
% desired, then value should be ~1e-5.
% key: 'rear temperature' >> value: double (= NaN by default). It fixes the 
% tempearture on the back side of the module. The temperature must be
% specified in Celsius.

%% Parsing optional inputs
rearConvCoef = 1;
rearAirTemp = NaN;
numargin = length(varargin);
if ~rem(numargin,2)%there must be a value for each key
    for k=1:2:numargin
        val = varargin{k+1};
        switch lower(varargin{k})
            case 'rear convection' %change if necessary
                if isnumeric(val)
                    rearConvCoef = val;
                else
                    error('Rear convection value must be double');
                end
            case 'rear temperature'
                if isnumeric(val)
                    rearAirTemp = val+273.15;%(K)
                else
                    error('Rear temperature value must be double');
                end
            otherwise
                try
                    error([varargin{k},' is an invalid input key']);
                catch %if you can't concatenate then the key is not a char
                    error('One of the keys is not a char');
                end
        end
    end
else
    error('Invalid input argument format');
end

%% Turning layers in PV module into structure to accelerate calcualtions
pvModule.layers = obj2struct(pvModule.layers,{'name','thickness','k','rho','cp','dx','dy','emissivity','nIntNodeY'});

%%
nNodesX = 1 + round(pvModule.width/pvModule.layers(1).dx);
nLayers = length(pvModule.layers);
% Count nodes in Y direction
nNodesY = (1+nLayers) + sum([pvModule.layers.nIntNodeY]);%#interfacial nodes + #internal nodes

%% Creating thermal network
thermalNet = ThermalNode.empty(0);%1D vector with thermal nodes
iTopLayer = nan(nNodesY,1);%row indices for top layer
iBotLayer = nan(nNodesY,1);%row indices for botom layer
i1=2;
i2=1;
for iLay = 1:nLayers
    iTopLayer(i1:i1+pvModule.layers(iLay).nIntNodeY) = iLay;
    iBotLayer(i2:i2+pvModule.layers(iLay).nIntNodeY) = iLay;
    i1 = i1 + pvModule.layers(iLay).nIntNodeY+1;
    i2 = i2 + pvModule.layers(iLay).nIntNodeY+1;
end
isInterfaceNode = iTopLayer~=iBotLayer;% Because at interface nodes the top and bottom layers are different

nodeTypeMat = nan(nNodesY,nNodesX); % Only for debugging purposes

% n is the index for the y direction (rows)
% m is the index for the x direction (columns)
% m and n start counting from the top left corner
for iNode = 1:nNodesX*nNodesY
    
    m = ceil(iNode/nNodesY);%Faster than [n,m] = ind2sub([nNodesY,nNodesX],iNode);
    n = 1+rem(iNode-1,nNodesY);%Faster than [n,m] = ind2sub([nNodesY,nNodesX],iNode);
    
    thermalNet(iNode) = ThermalNode();
    thermalNet(iNode).position = [n,m];
    
    % Assign layer info. For the top-most nodes there is no topLayer
    % data and for the bottom-most nodes there is no botLayer data
    if ~isnan(iTopLayer(n))
        thermalNet(iNode).topLayer = pvModule.layers(iTopLayer(n));
    else
        thermalNet(iNode).topLayer = PvLayer('name','None');%Some nodes don't have a top layer
    end
    if ~isnan(iBotLayer(n))
        thermalNet(iNode).botLayer = pvModule.layers(iBotLayer(n));
    else
        thermalNet(iNode).botLayer = PvLayer('name','None');%Some nodes don't have a bottom layer
    end
    
    %Assign node Type
    thermalNet(iNode).type = assignNodeType(n,m,nNodesY,nNodesX,isInterfaceNode);
    nodeTypeMat(iNode) = thermalNet(iNode).type; % Only for debugging purposes
    
end
thermalNet = obj2struct(thermalNet,{'position','type','topLayer','botLayer'});

%% Accessing arrays is faster than accesing struct members. The following
% arrays help to speed up the exetuion time
nodeType = [thermalNet.type];

layerNames = {'None',pvModule.layers.name};
layerIdxOrder = (1:length(layerNames))-1;

topLayer_iq =  nan(nNodesX*nNodesY,1,'double');%layer indices for Qgen
topLayer_thickness = nan(nNodesX*nNodesY,1,'double');
topLayer_k = nan(nNodesX*nNodesY,1,'double');
topLayer_rho = nan(nNodesX*nNodesY,1,'double');
topLayer_cp = nan(nNodesX*nNodesY,1,'double');
topLayer_dx = nan(nNodesX*nNodesY,1,'double');
topLayer_dy = nan(nNodesX*nNodesY,1,'double');
topLayer_emissivity = nan(nNodesX*nNodesY,1,'double');

botLayer_iq =  nan(nNodesX*nNodesY,1,'double');%layer indices for Qgen
botLayer_thickness = nan(nNodesX*nNodesY,1,'double');
botLayer_k = nan(nNodesX*nNodesY,1,'double');
botLayer_rho = nan(nNodesX*nNodesY,1,'double');
botLayer_cp = nan(nNodesX*nNodesY,1,'double');
botLayer_dx = nan(nNodesX*nNodesY,1,'double');
botLayer_dy = nan(nNodesX*nNodesY,1,'double');
botLayer_emissivity = nan(nNodesX*nNodesY,1,'double');

for iNode = 1:nNodesX*nNodesY
    topLayer_iq(iNode) = layerIdxOrder(strcmp(thermalNet(iNode).topLayer.name,layerNames));
    topLayer_thickness(iNode) = thermalNet(iNode).topLayer.thickness;
    topLayer_k(iNode) = thermalNet(iNode).topLayer.k;
    topLayer_rho(iNode) = thermalNet(iNode).topLayer.rho;
    topLayer_cp(iNode) = thermalNet(iNode).topLayer.cp;
    topLayer_dx(iNode) = thermalNet(iNode).topLayer.dx;
    topLayer_dy(iNode) = thermalNet(iNode).topLayer.dy;
    topLayer_emissivity(iNode) = thermalNet(iNode).topLayer.emissivity;
    botLayer_iq(iNode) = layerIdxOrder(strcmp(thermalNet(iNode).botLayer.name,layerNames));
    botLayer_thickness(iNode) = thermalNet(iNode).botLayer.thickness;
    botLayer_k(iNode) = thermalNet(iNode).botLayer.k;
    botLayer_rho(iNode) = thermalNet(iNode).botLayer.rho;
    botLayer_cp(iNode) = thermalNet(iNode).botLayer.cp;
    botLayer_dx(iNode) = thermalNet(iNode).botLayer.dx;
    botLayer_dy(iNode) = thermalNet(iNode).botLayer.dy;
    botLayer_emissivity(iNode) = thermalNet(iNode).botLayer.emissivity;
end
%%

dts = ones(1,length(weather.tAmb))*3600;
nInstants = length(dts);

qGen = weather.qGen;
tAmb = weather.tAmb+273.15;%(K)
tGnd = weather.tGnd+273.15;%(K)
tSky = weather.tSky+273.15;%(K)
wind = weather.windSpeed;%1*ones(1,nInstants);
moduleTemp = nan(nNodesY,nNodesX,nInstants,'double');%(K)
moduleTemp(:,:,1) = tAmb(1);%Initial module temperature guess
sigmaB = 5.670374419e-8;%Stefan-Boltzmann (W/m2/K4)
gravity = 9.8; %[m/s^2] Earth's gravity
epsFront = pvModule.layers(1).emissivity;
epsBack = pvModule.layers(end).emissivity;


dh = 2*pvModule.length*pvModule.width/(pvModule.length+pvModule.width);%hydraulic diameter
svf = (1+cosd(pvModule.tilt))/2;%sky view factor for front surface
gvf = (1-cosd(180-pvModule.tilt))/2;%ground view factor for back surface

for t = 1:nInstants-1
    dt = dts(t);
    Tfront = moduleTemp(1,:,t);
    Trear = moduleTemp(end,:,t);
    
    % Radiative heat transfer
    hRadFrontSky = sigmaB*(Tfront.^2+tSky(t)^2).*(Tfront+tSky(t))/((1-epsFront)/epsFront+1/svf);
    hRadRearGnd = sigmaB*(Trear.^2+tGnd(t)^2).*(Trear+tGnd(t))/((1-epsBack)/epsFront+1/svf);
    
    %% Convective heat transfer at Front side
    TFilmAirF = (Tfront+tAmb(t))/2;
    
    [rhoAirF,muAirF,prandtlAirF,kAirF,cpAirF] = calcAirPropK(TFilmAirF);
    nuAirF = muAirF./rhoAirF; % kinematic viscosity
    reynoldsF = wind(t)*dh./nuAirF;
    if reynoldsF<=5e5 %laminar regime
        c1F = 0.664; 
    else
        c1F = 0.86;
    end
    nusseltForcedF = c1F*reynoldsF.^(1/2).*prandtlAirF.^(1/3);
    grashofF = gravity./TFilmAirF.*abs(Tfront-TFilmAirF)*dh^3./nuAirF.^2;%Must be double-checked: I am using abs(deltaT)
    rayleighF = prandtlAirF.*grashofF;
    
    nusseltFreeF = nan(size(rayleighF));
    isHotPlate = Tfront>TFilmAirF;

    %https://doc.simulationx.com/4.2/1033/Content/Libraries/HeatTransfer/Accessories/BasicHeatTransfer/HeatTransferModels/FlatPlateHeatTransfer.htm?TocPath=Libraries%7CHeat%20Transfer%7CAccessories%7CHeat%20Transfer%7CHeat%20Transfer%20Models%7C_____8
    if pvModule.tilt>=0 && pvModule.tilt<= 5
        phiH = (1+(0.322./prandtlAirF).^(11/20)).^(-20/11);
        phiC = (1+(0.492./prandtlAirF).^(9/16)).^(-16/9);
        caseA = rayleighF.*phiH>7e4 & isHotPlate;
        caseB = rayleighF.*phiH<=7e4 & isHotPlate;
        caseC = ~isHotPlate;%cold Plate
        aux = 0.15*(rayleighF.*phiH).^(1/3);
        nusseltFreeF(caseA) = aux(caseA);
        aux = 0.766*(rayleighF.*phiH).^(1/5);
        nusseltFreeF(caseB) = aux(caseB);
        aux = 0.6*(rayleighF.*phiC).^(1/5);
        nusseltFreeF(caseC) = aux(caseC);
    elseif pvModule.tilt<= 85
        nusseltFreeF = (0.825+(0.387*(rayleighF*sind(pvModule.tilt)).^(1/6))./(1+(0.492./prandtlAirF).^(9/16)).^(8/27)).^2;
    else %vertical
        nusseltFreeF = 0.616*rayleighF.^(1/5).*(prandtlAirF./(0.8+prandtlAirF)).^(1/5);
    end

    nusseltFreeF(Tfront == TFilmAirF) = 0; % Zero convection heat transfer when mediums are at the same temperature
    if sum(isnan(nusseltFreeF))
        error('Invalid Nusselt number at front.');
    end
    hConvFreeF = nusseltFreeF.*kAirF/dh;
    hConvForcedF = nusseltForcedF.*kAirF/dh;
    hConvFront = (hConvFreeF.^3+hConvForcedF.^3).^(1/3);
    
    
    %% Convective heat transfer at Rear side
    if isnan(rearAirTemp)
        TFilmAirR = (Trear+tAmb(t))/2;
    else %if the air temperature on the back side is fixed
        TFilmAirR = (Trear+rearAirTemp)/2;
    end
    [rhoAirR,muAirR,prandtlAirR,kAirR,cpAirR] = calcAirPropK(TFilmAirR);
    nuAirR = muAirR./rhoAirR;% kinematic viscosity
    reynoldsR = (rearConvCoef*wind(t))*dh./nuAirR;
    if reynoldsR<=5e5 %laminar regime
        c1R = 0.664;
    else
        c1R = 0.86;%doi.org/10.1115/1.3450946
    end
    nusseltForcedR = c1R*reynoldsR.^(1/2).*prandtlAirR.^(1/3);
    grashofR = gravity./TFilmAirR.*abs(Trear-TFilmAirR).*dh^3./nuAirR.^2;%Must be double-checked: I am using abs(deltaT)
    rayleighR = prandtlAirR.*grashofR;
    nusseltFreeR = nan(size(rayleighR));
    isHotPlate = Trear>TFilmAirR;  
    %https://doc.simulationx.com/4.2/1033/Content/Libraries/HeatTransfer/Accessories/BasicHeatTransfer/HeatTransferModels/FlatPlateHeatTransfer.htm?TocPath=Libraries%7CHeat%20Transfer%7CAccessories%7CHeat%20Transfer%7CHeat%20Transfer%20Models%7C_____8
    if pvModule.tilt>=0 && pvModule.tilt<= 5
        phiH = (1+(0.492./prandtlAirR).^(9/16)).^(-16/9);
        phiC = (1+(0.322./prandtlAirR).^(11/20)).^(-20/11);
        caseA = isHotPlate;%hot Plate
        caseB = rayleighR.*phiC>7e4 & ~isHotPlate;
        caseC = rayleighR.*phiC<=7e4 & ~isHotPlate;
        
        aux = 0.6*(rayleighR.*phiH).^(1/5);
        nusseltFreeR(caseA) = aux(caseA);
        
        aux = 0.15*(rayleighR.*phiC).^(1/3);
        nusseltFreeR(caseB) = aux(caseB);
        
        aux = 0.766*(rayleighR.*phiC).^(1/5);
        nusseltFreeR(caseC) = aux(caseC);
    elseif pvModule.tilt<= 85
        nusseltFreeR = (0.825+(0.387*(rayleighR*sind(pvModule.tilt)).^(1/6))./(1+(0.492./prandtlAirR).^(9/16)).^(8/27)).^2;
    else %vertical
        nusseltFreeR = 0.616*rayleighR.^(1/5).*(prandtlAirR./(0.8+prandtlAirR)).^(1/5);
    end
    
    
    nusseltFreeR(Trear == TFilmAirR) = 0; % Zero convection heat transfer when mediums are at the same temperature
    if sum(isnan(nusseltFreeR))
        error('Invalid Nusselt number at rear.');
    end
    
    if rearConvCoef ~= 0 %if the rear side is completely insulated
        hConvFreeR = nusseltFreeR.*kAirR/dh;
    else
        hConvFreeR = zeros(size(nusseltFreeR));
    end
    hConvForcedR = nusseltForcedR.*kAirR/dh;
    hConvRear = (hConvFreeR.^3+hConvForcedR.^3).^(1/3);
    
    
    %% Generating matrix A and vector C to calculate module temp. A x T = C
    A = zeros(nNodesX*nNodesY,nNodesX*nNodesY);
    C = zeros(nNodesX*nNodesY,1);
    for iNode = 1:nNodesX*nNodesY
        m = ceil(iNode/nNodesY);% Faster than [n,m] = ind2sub([nNodesY,nNodesX],iNode);
        n = 1+rem(iNode-1,nNodesY);% Faster than [n,m] = ind2sub([nNodesY,nNodesX],iNode);
        switch nodeType(iNode)
            case 1 %top left node
                rho = botLayer_rho(iNode);
                cp = botLayer_cp(iNode);
                k = botLayer_k(iNode);
                dx  = botLayer_dx(iNode);
                dy  = botLayer_dy(iNode);
                
                BiRad = hRadFrontSky(m)*dy/k;
                BiY = hConvFront(m)*dy/k;
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = botLayer_iq(iNode);%layer index for Qgen
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + 2*tAmb(t)*BiY*FoY + 2*tSky(t)*BiRad*FoY + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX+2*BiY*FoY+2*BiRad*FoY;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -2*FoY;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -2*FoX;
            case 2
                rho = botLayer_rho(iNode);
                cp = botLayer_cp(iNode);
                k = botLayer_k(iNode);
                dx  = botLayer_dx(iNode);
                dy  = botLayer_dy(iNode);
                
                BiRad = hRadFrontSky(m)*dy/k;
                BiY = hConvFront(m)*dy/k;
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = botLayer_iq(iNode);%layer index for Qgen
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + 2*tAmb(t)*BiY*FoY + 2*tSky(t)*BiRad*FoY + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX+2*BiY*FoY+2*BiRad*FoY;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -2*FoY;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -FoX;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -FoX;
            case 3 %top right node
                rho = botLayer_rho(iNode);
                cp = botLayer_cp(iNode);
                k = botLayer_k(iNode);
                dx  = botLayer_dx(iNode);
                dy  = botLayer_dy(iNode);
                
                BiRad = hRadFrontSky(m)*dy/k;
                BiY = hConvFront(m)*dy/k;
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = botLayer_iq(iNode);%layer index for Qgen
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + 2*tAmb(t)*BiY*FoY + 2*tSky(t)*BiRad*FoY + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX+2*BiY*FoY+2*BiRad*FoY;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -2*FoY;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -2*FoX;
            case 4
                rho = botLayer_rho(iNode);
                cp = botLayer_cp(iNode);
                k = botLayer_k(iNode);
                dx  = botLayer_dx(iNode);
                dy  = botLayer_dy(iNode);
                
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = botLayer_iq(iNode);%layer index for Qgen

                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + moduleTemp(n,m,t);
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -FoY;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol) = -FoY;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -2*FoX;
            case 5 %all internal nodes
                rho = botLayer_rho(iNode);
                cp = botLayer_cp(iNode);
                k = botLayer_k(iNode);
                dx  = botLayer_dx(iNode);
                dy  = botLayer_dy(iNode);
                
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = botLayer_iq(iNode);%layer index for Qgen
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol)  = -FoY;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -FoY;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -FoX;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -FoX;
            case 6
                rho = botLayer_rho(iNode);
                cp = botLayer_cp(iNode);
                k = botLayer_k(iNode);
                dx  = botLayer_dx(iNode);
                dy  = botLayer_dy(iNode);
                
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = botLayer_iq(iNode);%layer index for Qgen
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol)  = -FoY;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -FoY;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -2*FoX;
            case 7
                rho1 = topLayer_rho(iNode);
                cp1 = topLayer_cp(iNode);
                k1 = topLayer_k(iNode);
                dx  = topLayer_dx(iNode);
                dy1  = topLayer_dy(iNode);
                
                rho2 = botLayer_rho(iNode);
                cp2 = botLayer_cp(iNode);
                k2 = botLayer_k(iNode);
                %dx  = botLayer_dx(iNode);
                dy2  = botLayer_dy(iNode);
                
                R1 = 2*k1*dt / (rho1*cp1*dy1^2 + rho2*cp2*dy1*dy2);
                R2 = 2*k2*dt / (rho1*cp1*dy1*dy2 + rho2*cp2*dy2^2);
                R3 = 2*k1*dy1*dt / (rho1*cp1*dx^2*dy1 + rho2*cp2*dx^2*dy2);
                R4 = 2*k2*dy2*dt / (rho1*cp1*dx^2*dy1 + rho2*cp2*dx^2*dy2);
                
                iqTop = topLayer_iq(iNode);%layer index for Qgen1
                iqBot = botLayer_iq(iNode);%layer index for Qgen1
                
                topQ = qGen(iqTop,t)*dy1*dt / (rho1*cp1*dy1 + rho2*cp2*dy2);
                botQ = qGen(iqBot,t)*dy2*dt / (rho1*cp1*dy1 + rho2*cp2*dy2);
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = topQ + botQ + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+R1+R2+R3+R4;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol)  = -R1;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -R2;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -(R3+R4);
            case 8
                rho1 = topLayer_rho(iNode);
                cp1 = topLayer_cp(iNode);
                k1 = topLayer_k(iNode);
                dx  = topLayer_dx(iNode);
                dy1  = topLayer_dy(iNode);
                
                rho2 = botLayer_rho(iNode);
                cp2 = botLayer_cp(iNode);
                k2 = botLayer_k(iNode);
                %dx  = botLayer_dx(iNode);
                dy2  = botLayer_dy(iNode);
                
                R1 = 2*k1*dt / (rho1*cp1*dy1^2 + rho2*cp2*dy1*dy2);
                R2 = 2*k2*dt / (rho1*cp1*dy1*dy2 + rho2*cp2*dy2^2);
                R3 = 2*dt / (rho1*cp1*dx^2*dy1 + rho2*cp2*dx^2*dy2) * (k1*dy1/2 + k2*dy2/2);
                
                iqTop = topLayer_iq(iNode);%layer index for Qgen1
                iqBot = botLayer_iq(iNode);%layer index for Qgen1
                
                topQ = qGen(iqTop,t)*dy1*dt / (rho1*cp1*dy1 + rho2*cp2*dy2);
                botQ = qGen(iqBot,t)*dy2*dt / (rho1*cp1*dy1 + rho2*cp2*dy2);
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = topQ + botQ + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+R1+R2+2*R3;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol)  = -R1;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -R2;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -R3;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -R3;
            case 9
                rho1 = topLayer_rho(iNode);
                cp1 = topLayer_cp(iNode);
                k1 = topLayer_k(iNode);
                dx  = topLayer_dx(iNode);
                dy1  = topLayer_dy(iNode);
                
                rho2 = botLayer_rho(iNode);
                cp2 = botLayer_cp(iNode);
                k2 = botLayer_k(iNode);
                %dx  = botLayer_dx(iNode);
                dy2  = botLayer_dy(iNode);
                
                R1 = 2*k1*dt / (rho1*cp1*dy1^2 + rho2*cp2*dy1*dy2);
                R2 = 2*k2*dt / (rho1*cp1*dy1*dy2 + rho2*cp2*dy2^2);
                R3 = 2*k1*dy1*dt / (rho1*cp1*dx^2*dy1 + rho2*cp2*dx^2*dy2);
                R4 = 2*k2*dy2*dt / (rho1*cp1*dx^2*dy1 + rho2*cp2*dx^2*dy2);
                
                iqTop = topLayer_iq(iNode);%layer index for Qgen1
                iqBot = botLayer_iq(iNode);%layer index for Qgen1
                
                topQ = qGen(iqTop,t)*dy1*dt / (rho1*cp1*dy1 + rho2*cp2*dy2);
                botQ = qGen(iqBot,t)*dy2*dt / (rho1*cp1*dy1 + rho2*cp2*dy2);
                
                iRow = (m-1)*nNodesY + n;
                C(iRow) = topQ + botQ + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+R1+R2+R3+R4;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol)  = -R1;
                % coefficient for [m,n+1]
                iCol = (m-1)*nNodesY + (n+1);
                A(iRow,iCol) = -R2;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -(R3+R4);
            case 10 %bottom left node
                rho = topLayer_rho(iNode);
                cp = topLayer_cp(iNode);
                k = topLayer_k(iNode);
                dx  = topLayer_dx(iNode);
                dy  = topLayer_dy(iNode);
                
                BiRad = hRadRearGnd(m)*dy/k;
                BiY = hConvRear(m)*dy/k;
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = topLayer_iq(iNode);%layer index for Qgen
                
                
                if isnan(rearAirTemp)
                    tRearSurf = tAmb(t);
                else
                    tRearSurf = rearAirTemp;
                end
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + 2*tRearSurf*BiY*FoY + 2*tGnd(t)*BiRad*FoY + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX+2*BiY*FoY+2*BiRad*FoY;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol) = -2*FoY;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -2*FoX;
            case 11
                rho = topLayer_rho(iNode);
                cp = topLayer_cp(iNode);
                k = topLayer_k(iNode);
                dx  = topLayer_dx(iNode);
                dy  = topLayer_dy(iNode);
                
                BiRad = hRadRearGnd(m)*dy/k;
                BiY = hConvRear(m)*dy/k;
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = topLayer_iq(iNode);%layer index for Qgen
                
                
                if isnan(rearAirTemp)
                    tRearSurf = tAmb(t);
                else
                    tRearSurf = rearAirTemp;
                end
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + 2*tRearSurf*BiY*FoY + 2*tGnd(t)*BiRad*FoY + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX+2*BiY*FoY+2*BiRad*FoY;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol) = -2*FoY;
                % coefficient for [m+1,n]
                iCol = ((m+1)-1)*nNodesY + n;
                A(iRow,iCol) = -FoX;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -FoX;
            case 12 %bottom right node
                rho = topLayer_rho(iNode);
                cp = topLayer_cp(iNode);
                k = topLayer_k(iNode);
                dx  = topLayer_dx(iNode);
                dy  = topLayer_dy(iNode);
                
                BiRad = hRadRearGnd(m)*dy/k;
                BiY = hConvRear(m)*dy/k;
                FoX = k/rho/cp*dt/dx^2;
                FoY = k/rho/cp*dt/dy^2;
                
                iq = topLayer_iq(iNode);%layer index for Qgen
                
                
                if isnan(rearAirTemp)
                    tRearSurf = tAmb(t);
                else
                    tRearSurf = rearAirTemp;
                end
                iRow = (m-1)*nNodesY + n;
                C(iRow) = qGen(iq,t)*dt/rho/cp + 2*tRearSurf*BiY*FoY + 2*tGnd(t)*BiRad*FoY + moduleTemp(n,m,t);
                
                % main diagonal (coefficient [m,n])
                iCol = (m-1)*nNodesY + n;
                A(iRow,iCol) = 1+2*FoY+2*FoX+2*BiY*FoY+2*BiRad*FoY;
                % coefficient for [m,n-1]
                iCol = (m-1)*nNodesY + (n-1);
                A(iRow,iCol) = -2*FoY;
                % coefficient for [m-1,n]
                iCol = ((m-1)-1)*nNodesY + n;
                A(iRow,iCol) = -2*FoX;
            otherwise
                error('Invalid node type');
        end
    end
    %Calculate Tmodule in t+1
    moduleTemp(:,:,t+1) = reshape(A\C,nNodesY,nNodesX);

end
moduleTemp = moduleTemp-273.15;%from K to C
thermalNet = reshape(thermalNet,nNodesY,nNodesX);
end


%% Temperature-dependent air properties
function [rhoAir,muAir,prAir,kAir,cpAir] = calcAirPropK(airTempK)
cpCoeff = [1.999999999999753e-04 -0.069999999999983 1.007999999999997e+03];
rhoCoeff = [2.20978564347549e-16,-6.91544160945619e-13,8.95101570636355e-10,-6.19895332646621e-07,0.000247648193236670,-0.0570020732669734,6.99490734259050];
muCoeff = [-3.683200000000047e-11 7.118832160000028e-08 6.639794524799587e-07];
prCoeff = [4.797099999999812e-07 -4.166255729999887e-04 0.789109669632473];
kCoeff = [-4.090500000000105e-08 9.903240150000067e-05 4.012593151373989e-04];

cpAir = cpCoeff(1)*airTempK.^2 + cpCoeff(2)*airTempK + cpCoeff(3);%Heat capacity (J/kg/K).
% tempC = airTempK - 273.15;
% RH = 80/100;%air relative humidity (-) This is the average relative humidity in the netherlands
% psat = 0.61121e3*exp((18.678-tempC/234.5).*(tempC./(257.14+tempC)));%Buck's equation (Pa)
% pwater = RH*psat;
% rhoAir = ((101325-pwater)*0.0289652+pwater*0.018016)./(8.31446*airTempK);
rhoAir = rhoCoeff(1)*airTempK.^6 + rhoCoeff(2)*airTempK.^5 + rhoCoeff(3)*airTempK.^4 + rhoCoeff(4)*airTempK.^3 + rhoCoeff(5)*airTempK.^2 + rhoCoeff(6)*airTempK + rhoCoeff(7);%Density (kg/m3)
muAir = muCoeff(1)*airTempK.^2 + muCoeff(2)*airTempK + muCoeff(3);%Dynamic viscosity (Pa*s)
prAir = prCoeff(1)*airTempK.^2 + prCoeff(2)*airTempK + prCoeff(3);%Prandtl number
kAir = kCoeff(1)*airTempK.^2 + kCoeff(2)*airTempK + kCoeff(3);%Thermal conductivity (W/m/K)
end

%% Node numbering
function nodeType = assignNodeType(n,m,nNodesY,nNodesX,isInterfaceNode)
if m==1 && n==1 %top left
    nodeType = 1;
elseif m<nNodesX && n==1 %top middle
    nodeType = 2;
elseif m==nNodesX && n==1 %top right
    nodeType = 3;
elseif m==1 && ~isInterfaceNode(n)
    nodeType = 4;
elseif m<nNodesX && ~isInterfaceNode(n)
    nodeType = 5;
elseif m==nNodesX && ~isInterfaceNode(n)
    nodeType = 6;
elseif m==1 && n<nNodesY && isInterfaceNode(n)
    nodeType = 7;
elseif m<nNodesX && n<nNodesY && isInterfaceNode(n)
    nodeType = 8;
elseif m==nNodesX && n<nNodesY && isInterfaceNode(n)
    nodeType = 9;
elseif m==1 && n==nNodesY %bottom left
    nodeType = 10;
elseif m<nNodesX && n==nNodesY %bottom middle
    nodeType = 11;
elseif m==nNodesX && n==nNodesY %bottom right
    nodeType = 12;
else
    error('Invalid node type');
end
end

function mystruct = obj2struct(objArray,fieldNames)
M=numel(objArray);
N=numel(fieldNames);

for m=M:-1:1
    for n=1:N
        field = fieldNames{n};
        mystruct(m).(field)=objArray(m).(field);
    end
end
end