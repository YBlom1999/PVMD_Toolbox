function [viModule, blockIsc] = simulateModIntepIVcurves(viDataset,gtDataset,numBypass,numCells,pvModType,pvModIrr,pvModTemp,varargin)
% Sintax:
% viModule = simulateModIntepIVcurves(viDataset,gtDataset,pvModLayoutMask,pvModConnMat,pvModType,pvModIrr,pvModTemp)
% [viModule,blockIsc] = simulateModIntepIVcurves(viDataset,gtDataset,pvModLayoutMask,pvModConnMat,pvModType,pvModIrr,pvModTemp,key1,value1,...)
%
% Description:
% It calculates the I-V curve of a PV module with a specific topology when
% the solar cells are at the specified irradiance and temperature conditions.
% The module's I-V curve is calculated by interpolating and adding together
% the I-V curves of the solar cells and bypass diodes. The function returns
% the I-V curve of the P-V module and the short-circuit current of each
% block of cells.
%
% Input arguments:
%
% viDataset: 3D array with the I-V curves of the solar cells at different
% temperatures and irradiances. The number of points per I-V curve is equal
% to the number of rows. Voltage and current values are in the first and
% second columns, respectively. In Each page we have an I-V curve for a
% different combination of temperature and irradiance.
%
% gtDataset: 2D array. The number of temperature and irradiance
% combinations (i.e, pages)in viDataset must be equal to the number of rows
% in gtDataset. Irradiance and temperature values are in the first and
% second columns, respectively. Both irradiance and temperature vectors must
% be linearly spaced (independently) and the values must be integers.
% If needed the inputs could be scaled to obtain resolutions better than
% 1 K and 1 W/m2. However, that's almost never necessary because the
% uncertainty in irradiance and temerpature values is always larger than
% 1 K and 1 W/m2.
%
% pvModType = 'sp' (series-parallel). In the future will be other options
% but for now you can only simulate series-parllel topologies with or without
% bypass diodes.
%
% pvModLayoutMask: 2D array to define the module layout, i.e. how many
% blocks of cells there are in the module. Blocks of cells in a SP
% (series-parallel) module are the strings of series-connected cells under
% the same bypass diode. If the module has no bypass diodes, then the blocks
% are each of the parallel-connected strings of cells.
% In TCT (total-cross-tied) modules, blocks are the subgroups of
% parallel-connected cells. The dimensions of pvModLayoutMask are determined
% by the number of cells in the module. Cells that belong to the same block are
% assigned the same integer value.
%
% pvModConnMat: 2D array to defines how the blocks are connected with each
% other. The dimensions of pvModConnMat are the same of pvModLayoutMask.
% When simulating SP modules, blocks with the same value are connected in
% series, then series-connected blocks are connected in parallel. When
% simulating TCT modules, blocks with the same value are connected in
% parallel, then parallel-connected blocks are connected in series.
%
% pvModIrr and pvModTemp: 3D arrays used to specify the operating
% conditions of each of the solar cells in the module. Each page in these
% arrays is a differnt module simulation. The number of rows and columns in
% these arrays must be the same as in  pvModLayoutMask and pvModConnMat.
% The temperature and irradiance value in each cell is used to select the
% cell's I-V curve from the viDataset by looking for the closest match in
% gtDataset.
%
% Optional keys and values:
%
% key: 'bypass' >> values: 'none' (default value) or 'ideal' or 2D double
% array or the forward I-V curve of a bypass element.
% To simulate a BPD with a specific I-V curve, use a 2D array with two
% columns with the forward I-V curve of the diode. The first column must
% contain the diode voltage (V) and the second the diode current (A). The
% diode voltage must start V = 0  and the diode maximum current value
% must be at least as high as the highest Isc in the solar cell's dataset.
%
% key: 'parallel mode' >> values: false (default), true (to start the parallel pool)
%
% key: 'verbose' >> values: false, true (default) (to display some information during the simulations)
%
% key: 'debugging' >> values: false (default), true (to plot the results of some of the intermediate steps)
%
% key: 'parallel points' >> number of points used for resampling I-V curves in parallel connections (default value = 2000)
%
% key: 'series points' >> number of points used for resampling I-V curves in series connections (default value = 2000)
%
% Output arguments:
%
% viModule: I-V curves of the module for each of the simulated instants.
% Each page corresponds to a different simulation. The first column are the
% module's voltage values and the second are module's current values.
%
% blockIsc: This is the short circuit current of each of the blocks in the
% module.
%
% Examples of use:
% All the following examples require the creation of an I-V dataset for the
% solar cell. A quick way to do this is the following:
%
% cellArea = 155;%(cm^2)
% R_tab_thick = 0.0122/16; %(ohm)
% R_rib_eff = 0.0028; %(ohm)
% Rs = 0.00273 + R_tab_thick + R_rib_eff; % Series resistance [ohms] << Includes the tabbing resistance
% Rp = 5.0067; % Shunt resistance[ohms]
% Is1=5.305e-11; % Saturation current [A]
% n1=1; %Ideality factor [-]
% Is2= 9.7e-6; % Saturation current [A]
% n2=2; %Ideality factor [-]
% Iph0 = 5.7818; % STC short circuit current [A]
% Eg = 1.12; %[eV]
% TIPH1= 0.00024; % Isc temp coeff [1/degC]
% TRP1=0;
% TRS1=0;
% TXIS1=3;
% TXIS2=3;
%
% minIrr  = 10;
% stepIrr = 10;
% maxIrr = 1200;
% minTemp = -10;
% stepTemp = 5;
% maxTemp = 40;
%
% [viDataset,gtDataset] = generateSCIVcurvesDDM(cellArea,Rs,Rp,Is1,n1,Is2,n2,Iph0,Eg, ...
%     TIPH1,TRP1,TRS1,TXIS1,TXIS2, ...
%     minIrr,stepIrr,maxIrr, ...
%     minTemp,stepTemp,maxTemp);
%
% 1) To simulate a conventional 72-cell PV module (all-series) with and 3
% bypass diodes at STC.
% pvModLayoutMask = [ones(2,12);2*ones(2,12);3*ones(2,12)];
% pvModConnMat = ones(6,12);
% irr = 1000*ones(6,12);
% temp  = 25*ones(6,12);
% pvModType = 'sp';
% viModule =  simulateModIntepIVcurves(viDataset,gtDataset,pvModLayoutMask,...
% pvModConnMat,pvModType,irr,temp,'bypass','ideal');
% plotIrrOnModule(irr);
% plotModuleIVcurves(viModule);
%
% 2) To simulate a partially shaded series-parallel module without bypass
% diodes and with the same layout as a conventional module:
% pvModLayoutMask = [ones(2,12);2*ones(2,12);3*ones(2,12)];
% pvModConnMat = [ones(2,12);2*ones(2,12);3*ones(2,12)];
% irr = 1000*ones(6,12);
% irr(1:3) = 100;
% temp  = 25*ones(6,12);
% pvModType = 'sp';
% viModule = simulateModIntepIVcurves(viDataset,gtDataset,pvModLayoutMask,...
% pvModConnMat,pvModType,irr,temp);
% plotIrrOnModule(irr);
% plotModuleIVcurves(viModule);
%
% 3) To simulate a partially shaded PV modules with 72 cells connected in
% series and one bypass diode per cell:
% pvModLayoutMask = reshape(1:72,6,12);
% pvModConnMat = ones(6,12);
% irr = 1000*ones(6,12);
% irr(1,1) = 100;
% irr(5,2:4) = 500;
% temp  = 25*ones(6,12);
% pvModType = 'sp';
% viModule = simulateModIntepIVcurves(viDataset,gtDataset,pvModLayoutMask,...
% pvModConnMat,pvModType,irr,temp,'bypass','ideal');
% plotIrrOnModule(irr);
% plotModuleIVcurves(viModule);


% NOTE FOR DEVELOPERS: This function can be further developed to simulate
% 'tct' and 'butterfly' topologies. The code for these 2 topologies must be
% rewritten following the logic of the 'sp' section. Doing this is quite
% simple and it works well because I have done it before and it worked very
% well. After making some improvements to handle I-V datasets of solar cells
% reverse characteristics, I did not update the code for the TCT and BUTTERFLY
% (HALF-CELL) topologies (I only did it for the SP topology). Hence, the
% code for these topologies became obsolete.
% Additional input parameter for developers:
% key: 'debugging' >> values: false (default), true (to plot internal I-V curves)

if numCells == 72 && numBypass == 3
    pvModLayoutMask = [ones(2,12);2*ones(2,12);3*ones(2,12)];
    pvModConnMat = ones(6,12);
else
    pvModLayoutMask = ones(6,12);
    pvModConnMat = ones(6,12);
end
if numCells == 72
    pvModIrr = reshape(pvModIrr',6,12,size(pvModIrr,1));
    pvModTemp = reshape(pvModTemp',6,12,size(pvModTemp,1));

end

[debugMode,parallelMode,bpdMode,bpdType,verbose,parPoints,serPoints] = parseOptionalArg(varargin);

if parallelMode
    parforArg = Inf;%To use as many threads as possible (machine-dependent)
else
    parforArg = 0;
end

isRevUndefined = true; %Check whether the cell's i-v curves include reverse characteristics
if sum(viDataset(1,1,:)<0) == size(viDataset,3)
    isRevUndefined = false;
end

if verbose
    fprintf('\n\n');
    disp('************************************************************');
    disp('******** Starting PV module I-V curve interpolation ********');
    disp('************************************************************');
    fprintf('\n');
end

% The irradiance and temperature of the solar cells must be integer values.
% If the input is not integer, then the irradiance and temperatures are
% converted to integers.
if ~isa(pvModIrr,'int16')
    aux = int16(pvModIrr);
    aux1 = pvModIrr-single(aux);
    aux1 = max(abs(aux1(:)));
    if verbose
        disp('<strong>Irradiance matrix converted to type int16.</strong>');
        fprintf('The largest cell irradiance rounding error is %.2f W/m2.\n\n',aux1);
    end
    pvModIrr = aux;
    clear aux1 aux;
end

if ~isa(pvModTemp,'int16')
    aux = int16(pvModTemp);
    aux1 = pvModTemp-single(aux);
    aux1 = max(abs(aux1(:)));
    if verbose
        disp('<strong>Temperature matrix converted to type int16.</strong>');
        fprintf('The largest cell temperature rounding error is %.2f K.\n\n',aux1);
    end
    pvModTemp = aux;
    clear aux1 aux;
end

%Check if the dimensions of pvModLayoutMask, irr and temp match
irrSize = size(pvModIrr);
tempSize = size(pvModTemp);
connMatSize = size(pvModConnMat);
if isequal(irrSize, tempSize)
    %Check if it is necesarry to transpose the module mask
    if ~isequal(irrSize(1:2),connMatSize)
        if isequal(flip(irrSize(1:2)),connMatSize)
            pvModConnMat = pvModConnMat';
            if verbose
                disp('<strong>Module mask has been transposed</strong>');
            end
        else
            error('Dimensions of irradiance matrix and module mask do not match');
        end
    end
else
    error('Dimensions of irradiance and temperature matrices do not match');
end

nCellGroups = max(pvModConnMat(:));%number of groups of cells in the module
nBlocks = max(pvModLayoutMask(:)); %total number of blocks of cells in the module - used to define where the bypass diodes are
nInstants = size(pvModIrr,3); %number of module I-V curves to be simulated

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

%Clipping input irradiance and temperature values to the limits defined by
%the gtDataset matrix
maxPVIrr =  max(pvModIrr(:));
maxPVTemp =  max(pvModTemp(:));
minPVIrr = min(pvModIrr(:));
minPVTemp = min(pvModTemp(:));

if verbose
    if maxPVIrr>max(irr)
        warning('The maximum module irradiance is %d W/m2 but the maximum irradiance in the cell dataset is %d W/m2.',maxPVIrr,max(irr));
    end
    if minPVIrr<min(irr)
        warning('The minimum module irradiance is %d W/m2 but the minimum irradiance in the cell dataset is %d W/m2.',minPVIrr,min(irr));
    end
    if maxPVTemp>max(temp)
        warning('The maximum module temperature is %d K but the maximum temperature in the cell dataset is %d K.',maxPVTemp,max(temp));
    end
    if minPVTemp<min(temp)
        warning('The minimum module temperature is %d K but the minimum temperature in the cell dataset is %d K.',minPVTemp,min(temp));
    end
end

%For each module irradiance and tempearture value, we look for index to
%the closest values in the temperature and irradiance vectors of the cell's
%I-V dataset. To avoid loops and since the module irradiance and
%temperature can be a 3D matrix, we use 4D matrices to find the indices.
irr4dim(1,1,1,:) = irr;
temp4dim(1,1,1,:) = temp;
[~,iIrr] = min(abs(pvModIrr-irr4dim),[],4);%indices to find closest value in the irradiance vector (gs)
iIrr = uint16(iIrr);
[~,iTemp] = min(abs(pvModTemp-temp4dim),[],4);%indices to find closest value in the temperature vector (vs)
iTemp = uint16(iTemp);

%Clip the module irradiance and temperature to the values irradiance and
%temperature values in the cell's I-V dataset.
% clippedIrr = gs(i4gs);
% clippedTemp = ts(i4ts);
% flatIrr = clippedIrr(:);
% flatTemp = clippedTemp(:);

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
iCellDataset = reshape(iCellDataset,size(pvModIrr,1),size(pvModIrr,2),size(pvModIrr,3));%rearraging the values make a 3D matrix again
clear irr temp clippedIrr clippedTemp i4gs i4ts % Don't need these variables anymore

blockIsc = single(nan(nInstants,nBlocks));

%%
% The debuggin mode is automatically disabled if the number of time
% instants to be simulated are higher than 20 to prevent the creation of
% thousands of plots when then user commits a mistake.
if nInstants>20
    debugMode = false;
end
%%

switch lower(pvModType)
    case {'series-parallel','sp','s-p'}
        %Resample the I-V curves to be able to "connect" them in series
        [i4s,v4s] = iv4series('vi',viDataset,'points',serPoints);%i4s is a 1-D array, v4s is a 2-D array of "voltage for series" connections
        viModule = single(zeros(parPoints,2,nInstants));



        %                 for iTime = 1:nInstants %Replace the parforline below with this line when debugging with breakpoints
        for iTime = 1:nInstants %Comment this line when debugging with breakpoints
            iCurve = iCellDataset(:,:,iTime);%Cells's I-V curve indexes
            %             iCurve = iCellDataset(:,:,1);%Cells's I-V curve indexes
            groupVolt = single(zeros(size(v4s,2),nCellGroups));%Voltage of each group of blocks
            if debugMode, f1=figure(); ax1=axes(f1,'NextPlot','add'); end
            blockIscAux = nan(1,nBlocks);
            for iGroup = 1:nCellGroups
                blocksInThisGroup = unique(pvModLayoutMask((pvModConnMat==iGroup))); % blocks in the current group
                nBlocksInGroup = length(blocksInThisGroup);%number of blocks in this group
                blockVolt = single(zeros(serPoints,nBlocksInGroup));%voltage of each block
                for iBlock = 1:nBlocksInGroup
                    iCurve2 =  iCurve(pvModLayoutMask==blocksInThisGroup(iBlock));%indices of I-V curves of the cells in the current block
                    blockVolt(:,iBlock) = sum(v4s(iCurve2,:),1);%connect in series all the cells
                    %If isRevUndefined blockVolt(:,iBlock) is NaN for all
                    %voltages that correspond to current values higher than Isc.
                    [~,iMinVolt] = min(abs(blockVolt(:,iBlock)));
                    blockVolt(iMinVolt,iBlock) = 0;%to extend the I-V curve to the y-axis
                    thisBlockIsc = i4s(iMinVolt);

                    blockIscAux(blocksInThisGroup(iBlock)) = thisBlockIsc;
                    if bpdMode%%Now add the I-V of the bpd
                        revFlag = false(size(blockVolt(:,iBlock)));
                        revFlag(iMinVolt+1:end) = true;%flags for which the voltage is negative
                        iaux = i4s(revFlag)-thisBlockIsc;
                        if strcmp(bpdType,'ideal')
                            blockVolt(revFlag,iBlock) = 0;
                        elseif ismatrix(bpdType)%user-specified BPD I-V curve [Voltage(V),Current(A)]
                            if isRevUndefined %the reverse I-V curve is defined only by the BPD
                                blockVolt(revFlag,iBlock) = -interp1(bpdType(:,2),bpdType(:,1),iaux,'pchip');
                            else %connect in parallel BPD and reverse cell I-V curve
                                voltBPD = -interp1(bpdType(:,2),bpdType(:,1),iaux,'pchip')';
                                voltRevCell = blockVolt(revFlag,iBlock);
                                currRevAux = i4s(revFlag)';
                                if currRevAux(1)>thisBlockIsc %extend the curves to vertical axis
                                    voltBPD = [0;voltBPD];
                                    voltRevCell = [0;voltRevCell];
                                    currRevAux = [thisBlockIsc;currRevAux];
                                end
                                isValidBPD = ~isnan(voltBPD);
                                isValidRevCell = ~isnan(voltRevCell);
                                vr1 = min(voltBPD(isValidBPD));
                                vr2 = min(voltRevCell(isValidRevCell));
                                voltRev = linspace(0,max(vr1,vr2),200)';%200 points are enough for this interpolation
                                currBPD = interp1(voltBPD(isValidBPD),currRevAux(isValidBPD),voltRev,'pchip')-thisBlockIsc;
                                currRevCell = interp1(voltRevCell(isValidRevCell),currRevAux(isValidRevCell),voltRev,'pchip')-thisBlockIsc;
                                currRev = currBPD+currRevCell+thisBlockIsc;%add currents and then shift up to Isc
                                blockVolt(revFlag,iBlock) = interp1(currRev,voltRev,i4s(revFlag),'pchip');
                            end
                        else % no bypass diodes (intentionally)
                            continue %only for clarity
                        end
                    end
                    if debugMode, plot(ax1,blockVolt(:,iBlock),i4s,'p-'); end
                    groupVolt(:,iGroup) = groupVolt(:,iGroup)+blockVolt(:,iBlock);%add each block in series to the corresponding group
                end

                if ~bpdMode && isRevUndefined% if there are no bypass diodes and the reverse I-V is not defined then the group's I-V curve is extended to the vertical axis
                    [~,myidx] = min(groupVolt(:,iGroup));
                    groupVolt(myidx,iGroup) = 0;
                end
                if debugMode, plot(ax1,groupVolt(:,iGroup),i4s,'b--'); end
            end
            blockIsc(iTime,:) = blockIscAux;
            [i4p,moduleVolt] = iv4parallel(i4s,groupVolt,parPoints);
            if debugMode, for iGroup=1:nCellGroups, plot(ax1,moduleVolt,i4p(iGroup,:),'x','MarkerSize',5); end, end
            moduleCurr = nansum(i4p,1);%connect groups in parallel to obtain the module's current
            if debugMode, plot(ax1,moduleVolt,moduleCurr,'r.-','LineWidth',2); end
            viModule(:,:,iTime) = [moduleVolt',moduleCurr'];

        end


    case {'total-cross-tied','tct','t-c-t'} %MUST BE REWRITTEN, IT IS NOT WORKING YET
        error('This code must be rewritten!');
        % wb = waitbar(0,'Starting simulation, interpolating I-V curves...');
        [i4p,v4p] = iv4parallel(viDataset,5000);%i4p is a matrix, v4p is a 1-D vector
        viModule = single(zeros(parPoints,2,nInstants));

        parfor iTime = 1:nInstants
            iCurve = iCellDataset(:,:,iTime);%Cells's I-V curve indexes
            i_g = single(zeros(size(i4p,2),nCellGroups));
            if debugMode f1=figure(); ax1=axes(f1); ax1.NextPlot='add'; end
            for iGroup=1:nCellGroups
                iv_cig =  iCurve(pvModConnMat==iGroup);%I-V curve index for each cell in the group
                i_g(:,iGroup) = sum(i4p(iv_cig,:),1);%Group's series voltage
                if debugMode plot(ax1,v4p,i_g(:,iGroup),'p'); end %plot the I-V curve of each group
            end

            %v_m and i_m are the I-V curve of the module
            [moduleCurr,v4s] = iv4series('i',i_g,'v',v4p,'points',parPoints);
            [~,iMinVolt] = min(v4s,[],2);
            v4s((iMinVolt-1)*size(v4s,1)+(1:length(iMinVolt))') = 0;

            if debugMode for iGroup=1:nCellGroups plot(ax1,v4s(iGroup,:),moduleCurr,'-','LineWidth',1.5); end, end %plot the resampled I-V curve of each group with solid lines

            v_gp = sum(v4s,1);
            [~,iMinVolt] = min(v_gp);
            v_gp(iMinVolt)=0;
            if debugMode
                figure();
                plot(v_gp,moduleCurr,'x');
            end
            viModule(:,:,iTime) = [v_gp',moduleCurr'];
            %             waitbar(iTime/tsteps,wb,['Simulation ',num2str(iTime),' of ',num2str(tsteps),' completed']);
        end
    case {'half cell','hc','butterfly'} %MUST BE REWRITTEN, IT IS NOT WORKING YET
        error('This code must be rewritten!');
        if ~bpdMode
            bpdType = 'ideal';
        end
        %         wb = waitbar(0,'Starting simulation, interpolating I-V curves...');
        [i4s,v4s] = iv4series('vi',viDataset,'points',serPoints);%i4s is a 1-D vector, v4s is a matrix
        viModule = single(zeros(parPoints,2,nInstants));
        parfor iTime = 1:nInstants
            iCurve = iCellDataset(:,:,iTime);%Cells's I-V curve indexes
            blockVolt = single(zeros(size(v4s,2),nCellGroups));
            if debugMode, f1=figure(); ax1=axes(f1); ax1.NextPlot='add'; end
            for iGroup=1:nCellGroups
                iv_cig =  iCurve(pvModConnMat==iGroup);%I-V curve index for each cell in the group
                blockVolt(:,iGroup) = sum(v4s(iv_cig,:),1);%Group's series voltage
                if debugMode, plot(ax1,blockVolt(:,iGroup),i4s,'p'); end
            end
            %v_m and i_m are the I-V curve of the module
            [i4p,v_gp] = iv4parallel(i4s,blockVolt,serPoints);
            if debugMode, for iGroup=1:nCellGroups,plot(ax1,v_gp,i4p(iGroup,:),'k.','MarkerSize',5); end, end
            %groups in parallel: 1+2, 3+4, 5+6
            i_g1 = sum(i4p([1 2],:));
            i_g2 = sum(i4p([3 4],:));
            i_g3 = sum(i4p([5 6],:));

            if debugMode, f2=figure();ax2=axes(f2); ax2.NextPlot = 'add'; plot(ax2,v_gp,i_g1,'o');plot(ax2,v_gp,i_g2,'o');plot(ax2,v_gp,i_g3,'o'); end

            [moduleCurr,v_gs] = iv4series('i',[i_g1',i_g2',i_g3'],'v',v_gp,'points',parPoints);

            for iGroup=1:size(v_gs,1)
                if isRevUndefined
                    [minVval,iMinVolt] = min(v_gs(:,iGroup));
                    if minVval >=0
                        v_gs(iMinVolt,iGroup) = 0;
                        revFlag = isnan(v_gs(:,iGroup)); %flags for which the voltage is negative
                    else
                        error('Interpolation Problem #AF201B.');
                    end
                else
                    revFlag = v_gs(iGroup,:)<=0 | isnan(v_gs(iGroup,:)); %flags for which the voltage is negative
                end
                [~,isc_ix] = min(abs(v_gs(iGroup,:)));
                thisBlockIsc = moduleCurr(isc_ix);%this is a good approximation
                iaux = moduleCurr(revFlag)-thisBlockIsc;
                if strcmp(bpdType,'ideal')
                    v_gs(iGroup,revFlag) = 0;
                end
            end


            if debugMode, plot(ax2,v_gs(1,:),moduleCurr,'k');plot(ax2,v_gs(2,:),moduleCurr,'k');plot(ax2,v_gs(3,:),moduleCurr,'k'); end

            moduleVolt = sum(v_gs,1);
            viModule(:,:,iTime) = [moduleVolt',moduleCurr'];
        end
    otherwise
        error('Invalid module type.');

end



if verbose
    fprintf('\n');
    disp('************************************************************');
    disp('******** Finished PV module I-V curve interpolation ********');
    disp('************************************************************');
    fprintf('\n');
end

end

%% Parsing optional arguments
function [debugMode,parallelMode,bpdMode,bpdType,verbose,parPoints,serPoints] = parseOptionalArg(varargin_)
debugMode = false;
parallelMode = false;
bpdMode = false;
bpdType = NaN;
verbose = true;
parPoints = 2000;
serPoints = 2000;

%Parsing optional input arguments
numargin = length(varargin_);
if ~rem(numargin,2)%there must be a value for each key
    for k=1:2:numargin
        val = varargin_{k+1};
        switch lower(varargin_{k})
            case 'debugging'
                if islogical(val)
                    debugMode = val;
                else
                    error('Debuggging value must be logical.');
                end
            case 'parallel mode'
                if islogical(val)
                    parallelMode = val;
                else
                    error('Parallel mode value must be logical.');
                end
            case 'verbose'
                if islogical(val)
                    verbose = val;
                else
                    error('Verbose value must be logical.');
                end
            case 'parallel points'
                if isscalar(val)
                    parPoints = round(val);
                else
                    error('Parallel points value must be a round scalar.');
                end
            case 'series points'
                if isscalar(val)
                    serPoints = round(val);
                else
                    error('Series points value must be a round scalar.');
                end
            case 'bypass'
                if (ischar(val) && strcmpi(val,'ideal')) || ismatrix(val)
                    if ismatrix(val) && isnumeric(val)
                        di_dv = diff(val(:,2))./diff(val(:,1));
                        if val(1,1)~= 0
                            error('The bypass diode I-V curve must start at V=0.');
                        elseif sum(di_dv<=0) || sum(isinf(di_dv))
                            error('The bypass diode I-V curve must be and biyective and increasing function.');
                        end
                    end
                    bpdType = val;
                    bpdMode = true;
                else
                    error('Bypass value error. Check the function documentation for valid values.');
                end
            otherwise
                try
                    error([varargin_{k},' is an invalid input key.']);
                catch %if you can't concatenate then the key is not a char
                    error('One of the keys is not a char.');
                end
        end
    end
else
    error('Invalid input argument format.');
end
end
