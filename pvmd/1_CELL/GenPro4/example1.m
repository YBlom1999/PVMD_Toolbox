%EXAMPLE 1: simulation of flat c-Si wafer (press F5 to run)

clear Lay               %clear variable 'Lay'

%===LAYERS (material & thickness [um])===
Lay(1).med = 'air';     Lay(1).thi = inf;      %inc. medium
Lay(2).med = 'c-Si';    Lay(2).thi = 200;      %c-Si wafer     
Lay(3).med = 'air';     Lay(3).thi = inf;      %outgoing medium
%===
%Input data taken from 'air.nk' and 'c-Si.nk' in the nk-folder

GENPRO4(Lay);           %run simulation