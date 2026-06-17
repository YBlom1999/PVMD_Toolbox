function V_rec = includeLumCouplingTandem(V_emit,V_emitAdd,V_rec,I,effLC)
% includeLumCouplingTandem simulates the effect of luminensence coupling
% and updates the voltage of the receiving cell based on the light of the
% top cell
%
% The method is based on the method of K. Jager, et al., Solar RRL 5 (3), https://doi.org/10.1002/solr.202000628
%
%
% Parameters
% ----------
% V_emit : double
%   The voltage of the emitting cell 
% V_emitAdd : double
%   The voltage of the additional emitting cell (nan if not applicable)
% V_rec : double
%   The voltage of the receiving cell
% I : double
%   The current of both cells
% effLC : double
%   The efficiency of the luminenscene coupling
% 
% Returns
% -------
% V_rec : double
%   The voltage of the receiving cell
%
% Developed by Y. Blom

if isnan(V_emitAdd)
    [~,Isc_ind] = min(abs(V_emit),[],2);
    Isc_emit = I(Isc_ind);
    I_rec_emit = max(Isc_emit'-I,0);
    for i = 1:size(V_rec,1)
        V_rec(i,:) = interp1(I,V_rec(i,:),(I-effLC*I_rec_emit(i,:))',"linear","extrap");
    end
else
    [~,Isc_ind] = min(abs(V_emit),[],2);
    Isc_emit = I(Isc_ind);
    I_rec_emit = max(Isc_emit'-I,0);

    [~,Isc_ind] = min(abs(V_emitAdd),[],2);
    Isc_emitAdd = I(Isc_ind);
    I_rec_emitAdd = max(Isc_emitAdd'-I,0);

    for i = 1:size(V_rec,1)
        V_rec(i,:) = interp1(I,V_rec(i,:),(I-effLC(1)*I_rec_emit(i,:)-effLC(2)*I_rec_emitAdd(i,:))',"linear","extrap");
    end

end
end