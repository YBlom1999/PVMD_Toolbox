function [V_mod3, I_mod3,chosenMASK3,chosen_config3] = NewReconfigAlgo...
    (V_block,I,MASK,cfg,day,baseline,deltaI_4allP,noP,deltaImpp,rec_algo)
%this fn defines a reconfiguration algorithm for a 96 cell PV module with 6 cell blocks (reconfigurable units)
%Input
%V_block - Cell block IV curves
%I - current values ocrresponding to V_block
%MASK, cfg - masks for the 27 possible configurations of the 96-cell 6-block PV module
%day - time instances with non-zero VI-cures
%baseline,baseline,deltaI_4allP,noP,deltaImpp,rec_algo - flags which need
%to be set to choose the reconfiguration algorithm of choice. The flags will be
%automatically set based on the user's choice of algorithm
%Output
%V_mod3 - Module voltage
%I_mod3 - Module current
%chosenMASK3,chosen_config3 - masks corresponding to the chosen
%configuration for each time instant
%Author: Devyani Salokhe

    blocks = max(cfg{1}); %no. of reconfigurable units in the PV module
    
    V_mod3 = zeros(length(day),length(I));
    I_mod3 = zeros(length(day),length(I));
    chosenMASK3 = zeros(size(MASK));
    chosen_config3 = {};
    
    %the reconfiguration algorithm starts here
    sensor_block = 1;
    
    %setting threshold values based on user's algorithm of choice. These 
    %are the optimum threshold values that were calculated for 4th May
    if baseline == 1 && rec_algo == 0
        deltaV = 0.0044;
        deltaI = 0.0009;
    elseif baseline == 1 && rec_algo == 1
        deltaV = 0.0705;
        deltaI = 16.379;
    elseif noP == 1 && rec_algo == 0
        deltaV = 0.0026;
        deltaI = 0.0009;
        MASK = MASK(:,:,2:size(cfg,2));
        cfg = cfg(2:size(cfg,2));
    elseif noP == 1 && rec_algo == 1
        deltaV = 0.0705;
        deltaI = 16.379;
    elseif deltaI_4allP == 1 && rec_algo == 0
        deltaV = 0.0089;
        deltaI = 0.0063;
    elseif deltaI_4allP == 1 && rec_algo == 1
        deltaV = 0.0019;
        deltaI = 0.0063;
    elseif deltaImpp == 1 && rec_algo == 0
        deltaV = 0.0089;
        deltaI = 0.01;
    elseif deltaImpp == 1 && rec_algo == 1
        deltaV = 0.0019;
        deltaI = 0.08;
    end
    
    chosen_config3{1} = {[1,2,3,4,5,6]};
    chosenMASK3(:,:,1) = MASK(:,:,size(cfg,2));
    V_check = zeros(length(day),1);
    I_check = zeros(length(day),1);
    recfg_flag = zeros(length(day),1);
    for h= 1:length(day)
        pointer = (h-1)*blocks;
        V_temp = V_block(pointer+1:pointer+blocks,:);
        if(sum(V_temp(~isnan(V_temp))) > 0)
            strings = max(chosenMASK3(:,:,h),[],'all');
            blocks_per_string = blocks/strings;
            V_block1 = V_block(pointer+1:pointer+blocks,:);
            %measure module V,I and determine Vmpp and Impp
            [V_mod3(h,:),I_mod3(h,:)] = calculate_IVmod(V_block1,I,chosenMASK3(:,:,h),chosen_config3{h});
            P = V_mod3(h,:) .* I_mod3(h,:);
            [~,mpp] = max(P);
            Vmpp = V_mod3(h,mpp);
            Impp = I_mod3(h,mpp);

            %determine block-1 V and I at MPP
            V_block1 = V_block(pointer+sensor_block,:);
            not_nan = ~isnan(V_block1);
            if iscell(chosen_config3{h})
                chosen_config_ = cell2mat(chosen_config3{h});
            else
                chosen_config_ = chosen_config3{h};
            end
            I_string = Block1IV(V_block(pointer+1:pointer+6,:),I,Vmpp,chosen_config_);
            [r,~] = find(chosen_config_ == sensor_block);
            Iblock1_mpp = I_string(r);
            if length(Iblock1_mpp) == 1
                Iblock1_mpp = repelem(Iblock1_mpp,length(find(not_nan)),1);...
                    %ADDED BY GSK!!!! NEEDS CONSULTATION
            end
            if Iblock1_mpp == 0
                Vblock1_mpp = 0; %ADDED BY GSK, NEEDS CONSULTATION!!!!!!
            else
                Vblock1_mpp = interp1(I(not_nan),V_block1(not_nan),Iblock1_mpp,'linear','extrap');
            end
            if length(Vblock1_mpp) > 1
                Vblock1_mpp = Vblock1_mpp(1); %ADDED BY GSK!!!! Needs consultation
            end
            Iblock1_mpp = Iblock1_mpp(1);
            %perform voltage and current check 
            if Vmpp
                V_check(h) = abs((Vmpp - blocks_per_string*Vblock1_mpp)/Vmpp) > deltaV;
            end
            
            if Impp
                if baseline
                    I_check(h) = abs((Impp - strings*Iblock1_mpp)/Impp) < deltaI;
                elseif deltaI_4allP && strings == blocks
                    I_check(h) = abs((Impp - strings*Iblock1_mpp)/Impp) < deltaI;
                elseif noP
                    I_check(h) = abs((Impp - strings*Iblock1_mpp)/Impp) < deltaI;
                elseif deltaImpp && strings == blocks
                    I_string_sorted = sort(I_string);
                    x_1_2 = (max(I_string_sorted(1:2))-min(I_string_sorted(1:2)))/max(I_string_sorted(1:2));
                    x_3_4 = (max(I_string_sorted(3:4))-min(I_string_sorted(3:4)))/max(I_string_sorted(3:4));
                    x_5_6 = (max(I_string_sorted(5:6))-min(I_string_sorted(5:6)))/max(I_string_sorted(5:6));
                    I_check(h) = max([x_1_2,x_3_4,x_5_6]) < deltaI;
                end
            end
            
            %decide if reconfiguration is necessary
            if V_check(h) || I_check(h)
                recfg_flag(h) = 1;
            end
            
            %measure module V and I post choosing optimum configuration
            if recfg_flag(h)
                if rec_algo
                    Isc_block = zeros(blocks,1);
                    %determine the short circuit current of each cell block
                    for i=1:blocks
                        Isc_block(i) = I(find(isnan(V_temp(i,:)),1)-1);
                    end
                    if ~noP
                        T1 = 0.03;
                        T2 = 0.05;
                        T3 = 0.08;
                        best_idx = estimateBestIdx(Isc_block,cfg,T1,T2,T3);
                    else
                        T1 = 0.03;
                        T2 = 0.05;
                        T3 = 10;
                        best_idx = estimateBestIdx(Isc_block,cfg,T1,T2,T3);
                    end
                    chosen_config3{h} = cfg{best_idx};
                    chosenMASK3(:,:,h) = MASK(:,:,best_idx);
                else
                    for m=1:size(MASK,3)
                        chosenMASK = MASK(:,:,m);
                        chosen_cfg = cfg{m};
                        [Vm(m,:),Im(m,:)] = calculate_IVmod(V_block(pointer+1:pointer+6,:),I,chosenMASK,chosen_cfg);
                    end
                    P = Vm .* Im;
                    Pmpp = max(P,[],2);
                    [~,mpp] = max(Pmpp);
                    chosenMASK3(:,:,h) = MASK(:,:,mpp);
                    chosen_config3{h} = cfg{mpp};
                end
            end
            %resample cell block IV curves to get module IV curve
            [V_mod3(h,:),I_mod3(h,:)] = calculate_IVmod(V_block(pointer+1:pointer+6,:),I,chosenMASK3(:,:,h),chosen_config3{h});
        end
        chosenMASK3(:,:,h+1) = chosenMASK3(:,:,h);
        chosen_config3{h+1} = chosen_config3{h};
    end
    
end

