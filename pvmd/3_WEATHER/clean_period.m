function [init_day,init_month,end_day,end_month] = ...
    clean_period(init_day,init_month,end_day,end_month,year)
%CLEAN_PERIOD Modify input period dates to make sure they are within the
%calendar limit
%
% Parameters
% ----------
% init_day : double
%   Initial day of the initial month for the simulations
% init_month : double
%   Initial month for the simulations
% end_day : double
%   Final day of the final month for the simulations
% end_month : double
%   Final month for the simulations
%
% Returns
% -------
% init_day : double
%   Initial day of the initial month for the simulations
% init_month : double
%   Initial month for the simulations
% end_day : double
%   Final day of the final month for the simulations
% end_month : double
%   Final month for the simulations

init_month = min(max(1,init_month),12);
init_day = min(max(1,init_day),eomday(year, init_month));
end_month = min(max(1,end_month),12);
end_day = min(max(1,end_day),eomday(year, end_month));
end

