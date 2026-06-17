function S = readDegradationSettings(S,SettingsFilePath)
% DegradationSettings generates a structure for the parameters of the
% analytical model.
%
% This function can be used to specify the parameters used for the lifetime
% calculation
% 
% Parameters
% ----------
%
% Returns
% -------
% S : structure
%   The settings for the analytical model
%
% Developed by by Youri Blom

set = readcell(SettingsFilePath,"Sheet","Degradation","Range","A2:B12");    %load settings file

for line = 1:size(set,1)
    S.(set{line,1}) = set{line,2};
end

% S.A_mois = 1e8;         % The pre-exponential constant of moisture induced degradation
% S.C_mois = 2;           % The exponent of moisture induced degradation
% S.Ea_mois = 0.809;      % The activation energy of moisture induced degradation
% S.A_dis = 0.12;         % The pre-exponential constant of discoloration
% S.Ea_dis = 0.39;        % The activation energy of discoloration
% S.factor_max = 0.5;     % The maximum change in minority carrier lifetimee due to LID
% S.C_LID = 1e7;          % The exponential decay constant for LID
% S.A_TC = 4e-5;          % The pre-exponential constant for thermal cyclying
% S.Ea_TC = 0.12;         % The activation energy for thermal cycling
% S.n_TC = 1.9;           % The exponent n of thermal cycling
% S.b_TC = 0.33;          % The exponent b of thermal cycling

end