function [outputCurrent,outputVoltages] = iv4series(varargin)
% Sintax:
% [outputCurrent,outputVoltages = iv4series(key1,value1,key2,value2,...)
% Description: This function is used to resample I-V curves to solve 
% series connections.
%
% Valid keys and values:
% key: 'points' >> value: integer (number of points used for interpolating
% the curves). Default value = 1000.
% key: 'vi' >> value: 3D double array. column 1: voltage, column 2:
% current. Each I-V curve goes in a different page.
% key: 'current limit' >> value: double. Default value: 20. If lower than
% the highest current in 'vi' it determines the interpolation range.
% key: 'v' >> value:
% key: 'i' >> value:
% Note: This function should receive with 'v' and 'i', or alternatively
% 'vi'. Mixing this optional arguments will throw an error.

nIVPairs = 1000;%default value
currentLim = 20;%default value
vis = [];
i_ = [];
v_ = [];
numargin = length(varargin);
if ~rem(numargin,2)%there must be a value for each key
    for iArg=1:2:numargin
        val = varargin{iArg+1};
        switch lower(varargin{iArg})
            case {'points'}
                if isnumeric(val) && val>1 && val==round(val)
                    nIVPairs = val;
                else
                    error('Points value must be a positive integer.');
                end
            case {'vi'}
                if isnumeric(val) && ndims(val) == 3
                    vis = val;
                else
                    error('vi must be a 3D array.');
                end
            case {'v'}
                if ismatrix(val)
                    v_ = val;
                else
                    error('v must be a 3D array.');
                end
            case {'i'}
                if ismatrix(val)
                    i_ = val;
                else
                    error('i must be a 3D array.');
                end
            case {'current limit'}
                if isnumeric(val) && val>0
                    currentLim = val;
                else
                    error('Current limit value must be a positive integer.');
                end
            otherwise
                try
                    error([varargin{iArg},' is an invalid input key.']);
                catch %if you can't concatenate then the key is not a char
                    error('One of the keys is not a char.');
                end
        end
    end
else
    error('Invalid input argument format.');
end

if ~isempty(vis) && isempty(i_) && isempty(v_)%This is the case when simulating SP modules
    maxI = min(max(max(vis(:,2,:))),currentLim);
    nCurves = size(vis,3);%number of urves to resample
    %95% of the points are used to sample the relevant part of the I-V
    %curve  in 1st and 2nd quadrants. The remaining 5% is used to
    %represent the part of the I-V curve in the 4th quadrant, which is
    %important to avoid problems when I-V curves are added together.
    outputCurrent = single([linspace(-maxI,-0.001,nIVPairs/20),linspace(0,maxI,nIVPairs/20*19)]);
    outputVoltages = single(zeros(nCurves,nIVPairs));
    for iCurve=1:nCurves
        curr = vis(:,2,iCurve);
        volt = vis(:,1,iCurve);
        [curr, index] = unique(curr); %to eliminate duplicated values
        outputVoltages(iCurve,:) = single(interp1(curr,volt(index),outputCurrent,'linear','extrap'));
    end
elseif ~isempty(i_) && ~isempty(v_) && isempty(vis)
    maxI = max(max(i_));
    nCurves = size(i_,2);%number of curves to resample
    outputCurrent = single(linspace(0,maxI,nIVPairs));
    outputVoltages = single(zeros(nCurves,nIVPairs));
    for iCurve=1:nCurves
        curr = i_(:,iCurve);
        if sum(curr.*v_')~=0 % don't interpolate if the P-V curve is all zeros
            [curr, index] = unique(curr); %to eliminate duplicated values
            volt = v_(index);
            outputVoltages(iCurve,:) = single(interp1(curr(~isnan(curr)),volt(~isnan(curr)),outputCurrent,'pchip'));
        else
            outputVoltages(iCurve,:) = 0;
        end
    end
else
    error('Invalid combination of inputs.');
end
end
