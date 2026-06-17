function [I_string2] = Block1IV(V_block2,I2,Vmpp2,chosen_config2)
%If the module is not in an all-series configuration, this function can be
%used to determine the current flowing through each parallel string.
%
%Input
%V_block2 - Cell block IV curves
%I2 - current values ocrresponding to V_block
%chosen_config2 - the configuration chosen by algorithm at each time step
%Vmpp2 - MPP voltage corresponding to the chosen configuration at each time
%step
%Output
%I_string2 - The current flowing through each parallel string
%Author - Devyani Salokhe
    
    strings = size(chosen_config2,1);
    V_string = zeros(strings,length(I2));
    I_string2 = zeros(strings,1);
    
    for s=1:strings
        blocks2 = chosen_config2(s,:);
        V_string(s,:) = sum(V_block2(blocks2,:),1);
        not_nan = ~isnan(V_string(s,:));
        if length(find(not_nan)) < 2
           I_string2(s) = 0; %LINES ADDED BY GSK (needs consultation)... 
            continue; 
        end
        if Vmpp2 >= min(V_string(s,not_nan))
            I_string2(s) = interp1(V_string(s,not_nan), I2(not_nan), Vmpp2,'linear','extrap');
        else
            I_string2(s) = interp1(V_string(s,not_nan), I2(not_nan), min(V_string(s,not_nan)),'linear','extrap');
        end
    end  
    
end