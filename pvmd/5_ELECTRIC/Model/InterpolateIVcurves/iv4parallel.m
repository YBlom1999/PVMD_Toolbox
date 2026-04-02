function [is,V] = iv4parallel(varargin)
%resample I-V curves to connect them in parallel.
%[is,V] = iv4parallel(ivs,POINTS)
%[is,V] = iv4parallel(I,vg,POINTS)

if length(varargin)==2 %one single matrix with all I-V curves
    iv = varargin{1};
    maxV = max(max(iv(:,1,:)));
    POINTS=varargin{2};%points in resampled curves
    NCURVES = size(iv,3);%number of urves to resample
    V = single(linspace(0,maxV,POINTS));
    is = single(zeros(NCURVES,POINTS));
    for k=1:NCURVES
        curr = iv(:,2,k);
        volt = iv(:,1,k);
        [volt, index] = unique(volt); %to eliminate duplicated values
        is(k,:) = interp1(volt,curr(index),V,'pchip');
    end
elseif length(varargin)==3 %input parameters (I,vg)
    I = varargin{1};
    vg = varargin{2};
    maxV = max(max(vg));
    minV = min(min(vg));
    POINTS=varargin{3};%points in resampled curves
    NCURVES = size(vg,2);%number of urves to resample
    if minV<0
        auxv = linspace(minV,0,POINTS/10+1);
        V = single([auxv(1:end-1),linspace(0,maxV,9/10*POINTS)]);%90% of the points for the forward characteristics
    else
        V = single(linspace(0,maxV,POINTS));
    end
    is = single(zeros(NCURVES,POINTS));
    for k=1:NCURVES
        volt = vg(:,k);
        if sum(volt.*I')~=0%don't interpolate if the P-V curve is all zeros
            [volt, index] = unique(volt); %to eliminate duplicated values
%             volt(1) = 0;% unique() eliminates points in the I-V curve in order of decreasing voltage, this line extends the range of the I-V curve back to 0 V 
            curr = I(index);
            is(k,:) = interp1(volt(~isnan(volt)),curr(~isnan(volt)),V,'pchip');
        else
            is(k,:) = 0;
        end
        
    end
else
    error('Invalid input')
end
end