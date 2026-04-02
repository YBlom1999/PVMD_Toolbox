function [V_mod2, I_mod2,chosenMASK2,chosen_config2] = SCCSReconfigAlgo(V_block,I,MASK,cfg,day)
%This function uses the reconfiguration algorithm developed by Andres Calcabrini to
%determine the optimum configuration for a 96-cell PV module.
%
%Input
%V_block - Cell block IV curves
%I - current values ocrresponding to V_block
%MASK, cfg - masks for the 27 possible configurations of the 96-cell 6-block PV module
%day - time instances with non-zero VI-cures
%Output
%V_mod2 - Module voltage
%I_mod2 - Module current
%chosenMASK2,chosen_config2 - masks corresponding to the chosen
%configuration for each time instant
%Author - Devyani Salokhe
 
%     rows = size(MASK(:,:,1),1);
%     cols = size(MASK(:,:,1),2);
    blocks = 6;
    Isc = zeros(blocks*length(day),1);
    %determine the short circuit current of each cell block
    for i=1:blocks*length(day)
        if sum(V_block(i,~isnan(V_block(i,:))))>0 
            if isnan(V_block(i,end))
               Isc(i) = I(find(isnan(V_block(i,:)),1)-1);
            else %LINES ADDED BY GSK, REQUIRES CONSULTATION
                V_block(i,end) = nan;
                Isc(i) = I(find(isnan(V_block(i,:)),1)-1);
            end
        end
    end
    Isc_block = reshape(Isc,blocks,length(day));

    %determine the best configuration
    T1 = 0.03;
    T2 = 0.05;
    T3 = 0.08;
    best_idx = estimateBestIdx(Isc_block,cfg,T1,T2,T3);

    chosen_config2 = {};
    chosenMASK2 = zeros(size(MASK));
    I_mod2 = zeros(length(day),length(I)); 
    p = zeros(length(day),1);
    V_mod2 = zeros(length(day),length(I));

    for ii=1:length(day)
        if (all(Isc_block(:,ii)))
            pointer = (ii-1)*blocks;
            chosen_config2{ii} = cfg{best_idx(ii)};
            p(ii) = size(chosen_config2{ii},1);
            chosenMASK2(:,:,ii) = MASK(:,:,best_idx(ii));
            
            %Resample the cell-block IV curves to get the module IV curves 
            [V_mod2(ii,:),I_mod2(ii,:)] = calculate_IVmod(V_block(pointer+1:pointer+6,:),I,chosenMASK2(:,:,ii),chosen_config2{ii});
        end    
    end
end
