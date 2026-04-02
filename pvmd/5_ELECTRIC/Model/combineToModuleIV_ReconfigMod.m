function [V_module,reconfig_setting] = combineToModuleIV_ReconfigMod(V_,numCells,I,day,night,R,Rsw,algo,nc_r)
%%Combines all cell VI curves into one module IV curve for reconfigurable
%%modules
%Input:
%V_ -Cell VI curves
%numCells - number of cells
%I - current values ocrresponding to V_
%day - time instances with non-zero VI-cures
%night - time instances with  VI-cures=0
%R - Series resistance due to metallization and shading (this R value
%should be the total resistance for 16 cells)
%Inputs from user:
%Rsw - Resistance of the switches used in the switching circuit
%algo - Reconfiguration algorithm of choice
%nc_r - number of cell rows
%Ouput: Module VI curve
%Author: Malte Vogt, Devyani Salokhe
%Edited for PVMD toolbox by Karthik Ganapathi Subramanian

%Remove contact resistance from IV curve
V_=V_ - I*R;

%create all masks for the 6-block,96-cell reconfigurable module
blocks = 6; %total number of blocks
x = ones(4,4); %block configuration
[MASK,cfg] = createAllMasks(x,nc_r);

%calculate the block voltages
if nc_r == 8
    BlockMask =  [1*x 2*x 3*x; 4*x 5*x 6*x]; %for landscape orientations
else
   BlockMask =  [1*x 2*x 3*x; 4*x 5*x 6*x]'; %for portrait orientations
end
BlockMask = reshape(BlockMask, numCells, 1);
V_block = zeros((size(V_,1)/numCells)*blocks , length(I));
for h=1:size(day)
    V_cells = V_((h-1)*numCells+1:h*numCells , :);
    for b=1:blocks
        bx = (h-1)*blocks+b;
        V_block(bx,:) = sum(V_cells(BlockMask==b,:)) - I*(2*Rsw); 
    end
end
V_block(V_block<0) = nan; %ADDED BY ME!!! (GSK)

%calculate module IV as per chosen reconfiguration algorithm
 
if algo == 1
    [V__,~,chosenMASK,chosen_config] = SCCSReconfigAlgo(V_block,I,MASK,cfg,day);
else
    baseline = 0;
    deltaI_4allP = 0;
    noP = 0;
    deltaImpp = 0;
    rec_algo = 0;
    if algo == 2
        baseline = 1;
    elseif algo == 3
        baseline = 1;
        rec_algo = 1;
    elseif algo == 4
        noP = 1;
    elseif algo == 5
        noP = 1;
        rec_algo = 1;
    elseif algo == 6
        deltaI_4allP = 1;
    elseif algo == 7
        deltaI_4allP = 1;
        rec_algo = 1;
    elseif algo == 8
        deltaImpp = 1;
    elseif algo == 9
        deltaImpp = 1;
        rec_algo = 1;
    end
    [V__,~,chosenMASK,chosen_config] = NewReconfigAlgo(V_block,I,MASK,cfg,day,...
                                            baseline,deltaI_4allP,noP,deltaImpp,rec_algo);
end

%Add 0 curves
V_module=-999*ones(length(day)+length(night),size(I,2));
V_module(night,:)=zeros(length(night),length(I));
V_module(day,:)=V__;
V_module=max(V_module,0);
reconfig_setting.ChosenMASK = chosenMASK;
reconfig_setting.ChosenCfg = chosen_config;
end