function S = lux_settings2

S.parll = 1;        %parallel computation (0 = no, 1 = yes)

%---settings for pruning the ray tree (speed vs. accuracy)---
S.Imin = 1e-4;      %minimum intensity of ray before termination (WEAK RAYS)
S.Imc = 0.5;        %minimum intensity of ray before switching to Monte-Carlo mode
S.mx_g = 100;       %maximum generation of ray before termination (OLD RAYS)

%---plot settings---
S.plt = 1;          %0 = plot b&w geometry,
                    %1 = plot output (calculated sensitivity) as color in geometry
                    %2 = plot b&w geometry with rays*
%* rays will not show in parallel computation
 
end