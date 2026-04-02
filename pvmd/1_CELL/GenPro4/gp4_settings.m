function S = gp4_settings

S.wav = 0.300:0.020:1.200;  %default wavelength range (used when wav-vector is not given as input)

S.nai = 18;         %number of angular intervals (5 = fast, 30 = accurate)
S.iai = 1;          %incident light angular interval (1 = 0 deg, nai = 90 deg)
S.parll = 1;        %parallel computation for speed (1 = yes, 0 = no)
S.plot = 2;         %plot results (0 = none, 1 = TAR, 2 = ray-texture, 3 = SM, 4 = rays) 3&4 only show when parll = 0
S.spec = 'AM1.5g';  %name of spectrum in spectrum.mat file
S.ELx = 1e-12;      %warning threshold energy non-conservation
S.sd = 1;           %skip dark interface (1 = yes, 0 = no)

%===ray-tracing settings===
S.n_rays = 500;    %nr of rays (10 = fast, 1000 = accurate)
S.mxx = 1000;       %maximum nr. of intersection for 1 ray
S.imc = 0.20/S.n_rays; %monte-carlo ray intensity threshold;
S.imn = 1e-6;       %minimum ray intensity for termination
S.ETx = 5e-2;       %warning threshold total terminated ray int.
S.SM = 1;           %give 'ScatterMatrix' data as output to reuse for speed (1 = yes, 0 = no)
S.AT = 0;           %give 'AngleTree' ray-data as output to reuse for speed (1 = yes, 0 = no) works with 'ray2' only

end