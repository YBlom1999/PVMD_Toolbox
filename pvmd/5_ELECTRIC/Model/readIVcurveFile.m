function [paramForwardBias,paramReverseBias,paramLightSoaking] = readIVcurveFile(IVcurvefile)
% readIVcurveFile reads the IV curve file and extracts the parameters for
% forward bias, reverse bias, and light soaking
%
% Parameters
% ----------
% IVcurvefile : string
%   The name of the IV curve file
%
% Returns
% -------
% paramForwardBias: struct
%   The parameters for forward bias
% paramReverseBias: struct
%   The parameters for reverse bias
% paramLightSoaking: struct
%   The parameters for light soaking
%
% Developed by Youri Blom

IVcurveFilePath = [IVcurvefile,'.xlsx'];

%read settings for forward bias
settingsForward = readcell(IVcurveFilePath,"Sheet","ForwardBias","Range","B2:F11");
paramForwardBias.Jph_J = flip(cell2mat(settingsForward(1,:)));
paramForwardBias.Jph_T = flip(cell2mat(settingsForward(2,:)));
paramForwardBias.J0_J = flip(cell2mat(settingsForward(3,:)));
paramForwardBias.J0_T = flip(cell2mat(settingsForward(4,:)));
paramForwardBias.n_J = flip(cell2mat(settingsForward(5,:)));
paramForwardBias.n_T = flip(cell2mat(settingsForward(6,:)));
paramForwardBias.Rsh_J = flip(cell2mat(settingsForward(7,:)));
paramForwardBias.Rsh_T = flip(cell2mat(settingsForward(8,:)));
paramForwardBias.Rs_J = flip(cell2mat(settingsForward(9,:)));
paramForwardBias.Rs_T = flip(cell2mat(settingsForward(10,:)));

%read settings for reverse bias
settingsReverse = readcell(IVcurveFilePath,"Sheet","ReverseBias","Range","B1:F16");
paramReverseBias.modelType = cell2mat(settingsReverse(1,1));
if paramReverseBias.modelType == 1
    paramReverseBias.Be = cell2mat(settingsReverse(4,1));
    paramReverseBias.phi_t = cell2mat(settingsReverse(5,1));
    paramReverseBias.V_b = cell2mat(settingsReverse(6,1));
    paramReverseBias.c = cell2mat(settingsReverse(7,1));
elseif paramReverseBias.modelType == 2
    paramReverseBias.J0_rev_J = flip(cell2mat(settingsForward(11,:)));
    paramReverseBias.J0_rev_T = flip(cell2mat(settingsForward(12,:)));
    paramReverseBias.Vs_rev_J = flip(cell2mat(settingsForward(13,:)));
    paramReverseBias.Vs_rev_T = flip(cell2mat(settingsForward(14,:)));
    paramReverseBias.Rsh_rev_J = flip(cell2mat(settingsForward(15,:)));
    paramReverseBias.Rsh_rev_T = flip(cell2mat(settingsForward(16,:)));
end

%read settings for light soaking
settingsLightSoaking = readcell(IVcurveFilePath,"Sheet","LightSoaking","Range","B2:B4");
paramLightSoaking.A = cell2mat(settingsLightSoaking(1,1));
paramLightSoaking.B = cell2mat(settingsLightSoaking(2,1));
paramLightSoaking.Ea = cell2mat(settingsLightSoaking(3,1));
end