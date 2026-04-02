%EXAMPLE 2: simulation of flat c-Si heterojunction solar cell

clear Lay Int            %clear variables 'Lay' and 'Int'

%===LAYERS===
Lay(1).med = 'air';              Lay(1).thi = inf;
Lay(2).med = 'c-Si';             Lay(2).thi = 200;
Lay(3).med = 'air';              Lay(3).thi = inf;

%===INTERFACES===        coatings are part of the interface
%interface 1: between layer 1 and 2 (air/c-Si)
Int(1).coat(1).med = 'ITO';      Int(1).coat(1).thi = 0.080;    %front
Int(1).coat(2).med = 'a-Si(p)';  Int(1).coat(2).thi = 0.005;
Int(1).coat(3).med = 'a-Si(i)';  Int(1).coat(3).thi = 0.005;
%interface 2: between layer 2 and 3 (c-Si/air)
Int(2).coat(1).med = 'a-Si(i)';  Int(2).coat(1).thi = 0.005;
Int(2).coat(2).med = 'a-Si(n)';  Int(2).coat(2).thi = 0.005;
Int(2).coat(3).med = 'Ag';       Int(2).coat(3).thi = 0.300;    %rear
%Note: interference occurs in 'coatings', not in 'layers'.
%===

[Lay,Int,out] = GENPRO4(Lay,Int);        %both 'Lay' and 'Int' are input