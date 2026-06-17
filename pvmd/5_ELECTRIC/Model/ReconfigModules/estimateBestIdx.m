function best_idx = estimateBestIdx(isc_,cfg,T1,T2,T3)
% isc_ must have 6 rows and as many columns as time instants
%Author - Andres Calcabrini (added to PVMD toolbox by Devyani Salokhe)

for tix=1:size(isc_,2)
    [isc,idx] = sort(isc_(:,tix));
        aux = [idx(1) idx(2); idx(3) idx(4); idx(5) idx(6)];
    for k=2:16 %find the bst s2p3 configuration
        if isequal(sort(aux(1,:)),sort(cfg{k}(1,:)))
            if isequal(sort(aux(2,:)),sort(cfg{k}(2,:)))
                if isequal(sort(aux(3,:)),sort(cfg{k}(3,:)))
                    best_s2p3 = k;
                end
            elseif isequal(sort(aux(2,:)),sort(cfg{k}(3,:)))
                if isequal(sort(aux(3,:)),sort(cfg{k}(2,:)))
                    best_s2p3 = k;
                end
            end
        elseif isequal(sort(aux(1,:)),sort(cfg{k}(2,:)))
            if isequal(sort(aux(2,:)),sort(cfg{k}(1,:)))
                if isequal(sort(aux(3,:)),sort(cfg{k}(3,:)))
                    best_s2p3 = k;
                end
            elseif isequal(sort(aux(2,:)),sort(cfg{k}(3,:)))
                if isequal(sort(aux(3,:)),sort(cfg{k}(1,:)))
                    best_s2p3 = k;
                end
            end
        elseif isequal(sort(aux(1,:)),sort(cfg{k}(3,:)))
            if isequal(sort(aux(2,:)),sort(cfg{k}(2,:)))
                if isequal(sort(aux(3,:)),sort(cfg{k}(1,:)))
                    best_s2p3 = k;
                end
            elseif isequal(sort(aux(2,:)),sort(cfg{k}(1,:)))
                if isequal(sort(aux(3,:)),sort(cfg{k}(2,:)))
                    best_s2p3 = k;
                end
            end
            
        end
    end
      
    aux = [idx(1) idx(2) idx(3); idx(4) idx(5) idx(6)];
    for k=17:26 %find the bst s3p2 configuration
        if isequal(sort(aux(1,:)),sort(cfg{k}(1,:))) || isequal(sort(aux(2,:)),sort(cfg{k}(1,:)))
            best_s3p2 = k;
            break;
        end
    end
    
    x_1_6 = (max(isc(1:6))-min(isc(1:6)))/max(isc(1:6));
    x_1_3 = (max(isc(1:3))-min(isc(1:3)))/max(isc(1:3));
    x_4_6 = (max(isc(4:6))-min(isc(4:6)))/max(isc(4:6));
    x_1_2 = (max(isc(1:2))-min(isc(1:2)))/max(isc(1:2));
    x_3_4 = (max(isc(3:4))-min(isc(3:4)))/max(isc(3:4));
    x_5_6 = (max(isc(5:6))-min(isc(5:6)))/max(isc(5:6));
    if x_1_6<T1
        best_idx(tix) = 27;%s6p1 (all series)
    else
        if max([x_1_3,x_4_6])<T2
            best_idx(tix) = best_s3p2;
        else
            if max([x_1_2,x_3_4,x_5_6])<T3
                best_idx(tix) = best_s2p3;
            else
                best_idx(tix) = 1;%s1p6 (all parallel)
            end
        end
    end
    
end
end