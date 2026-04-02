function [Vmod,Imod] = calculate_IVmod(V_block1,I1,chosenMASK1,chosen_config1)
%This function resamples the cell-block IV curves to calculate module IV
%curves.
%Input
%V_block,I - cell-block IV
%chosenMASK,chosen_config - masks corresponding to the chosen configuration
%R - Resistance
%Output
%Vmod, Imod - module IV

    nCells = numel(chosenMASK1);
    chosenMASK1 = reshape(chosenMASK1, nCells, 1);
    p = max(chosenMASK1); %no. of parallel strings in chosen configuration
    
    %calculate the voltage of each parallel strings
    if p == 1
        Vmod = sum(V_block1);
        Imod = I1;
    else 
        V_string = zeros(p,length(I1));
        min_V_string = zeros(p,1);
        for px=1:p
            V_string(px,:) = sum(V_block1(chosen_config1(px,:),:),1);
            min_V_string(px) = min(V_string(px,V_string(px,:)>=0));
        end

        Vmod = max(V_string);
        for v=1:length(I1)
            if any(Vmod(v) < min_V_string)
                Vmod(v) = nan;
            end
        end
        
        I_string = nan(p,length(I1));
        for kk=1:p
            %calculate each string's I wrt V_mod
            not_nan_Vstring = ~isnan(V_string(kk,:));
            pts = Vmod >= min(V_string(kk,not_nan_Vstring));
            if sum(not_nan_Vstring) < sum(pts) ||...
                    length(find(pts)) == 1
                I_string(kk,pts) = 0;
                continue; %ADDED BY GSK, Needs consultation.... 
            end
            I_string(kk,pts) = interp1(V_string(kk,not_nan_Vstring),...
                    I1(not_nan_Vstring),Vmod(pts),'linear','extrap');

            %extend the IV curve to y-axis
            nan_I = isnan(I_string(kk,:));
            not_nan_Vmod = ~isnan(Vmod);
            I_string(kk,nan_I&not_nan_Vmod) = max(I_string(kk,:)); 
        end
        Imod = sum(I_string,1);
        nanx = isnan(Imod);
        stepp = 5e-3*p;
        %this step is done to equate I-mod's length with I's length
        Imod(nanx) = max(Imod)+stepp:stepp:max(Imod)+stepp*sum(nanx);
    end
    Vmod = max(Vmod,0);
end
