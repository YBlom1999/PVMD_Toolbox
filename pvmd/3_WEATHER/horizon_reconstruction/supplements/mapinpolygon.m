% The gridded structure of the height map allows for the use of the 
% poly2mask() function from the Image Processing Toolbox which is much 
% faster than the more commonly used inpolygon() function, which is more 
% useful when the query points are randomly spaced. Special attention is 
% needed to understand the input requirements of poly2mak(), which was 
% created for identifying pixels inside a shape. The left most value has
% pixel index 1 and right most pixel index is equal to horizontal length.
% This means that the original x data needs to be scaled and shifted. The
% same holds for y data.

function [H,bbox,X,Y] = mapinpolygon(Hgrid,BBOXgrid,Xpoly,Ypoly,varargin)
% With the MxN gridded height map matrix (Hgrid) and its x and y limits 
% (BBOXgrid) as input, this function gives all elements outside the polygon 
% (Xpoly,Ypoly) a NaN value. In this version only a grid spacing of 0.25 or 
% 0.5 is allowed. If the polygon is a box, a fifth input (varargin) can be
% given a true boolean input for faster computation. The output H and bbox
% are cropped to the maximum and minimum values of x and y that are on the
% gridlines carried over from Hgrid. The outputs X and Y are the 
% coordinates in the same size as H, which are neglected when left empty.
% Author: Maarten Verkou
% Update: 16/06/2020


% If last element is true than the function is used for a box
if isempty(varargin)
    isBox = false;
elseif varargin{1}==true
    isBox = varargin{1};
else
    error('Undefined parameter');
end

% Make sure poly is of same size and remove if last value is NaN
if size(Xpoly)==size(Ypoly)
    if ~isnan(Xpoly(end)) && ~isnan(Ypoly(end))
        Xpoly(end+1) = NaN; Ypoly(end+1) = NaN;
    end
else
    error('X and Y points must have same dimensions!')
end

% Get X vector: east->west (Low->High)
xlin = linspace(BBOXgrid(1,1),BBOXgrid(2,1),size(Hgrid,2));
xstp = interp1([0.25,0.50],[0.25,0.50],abs(xlin(2)-xlin(1)),...
    'nearest','extrap');
if abs(xlin(1)-xlin(2)) ~= xstp 
%     warning('possible error in Height Map step size')
    xlin = BBOXgrid(1,1):xstp:BBOXgrid(2,1);
    if size(xlin,2)~=size(Hgrid,2)
%         Hgrid = Hgrid(:,1:end-1);
        xlin = xlin(1:end-1);
    end
end
    
% Get Y vector north->south (High->Low)
ylin = linspace(BBOXgrid(2,2),BBOXgrid(1,2),size(Hgrid,1));
ystp = interp1([0.25,0.50],[0.25,0.50],abs(xlin(2)-xlin(1)),...
    'nearest','extrap');
if abs(ylin(1)-ylin(2)) ~= ystp 
%     warning('possible error in Height Map step size')
    ylin = BBOXgrid(2,2):-ystp:BBOXgrid(1,2);
    if size(ylin,2)~=size(Hgrid,1)
%         Hgrid = Hgrid(1:end-1,:);
        ylin = ylin(1:end-1);
    end   
end
% Find border of polygon
x1 = find(xlin>min(Xpoly),1,'first'); x2 = find(xlin<max(Xpoly),1,'last');
y1 = find(ylin<max(Ypoly),1,'first'); y2 = find(ylin>min(Ypoly),1,'last');
if isempty(x1) || isempty(x2) || isempty(y1) || isempty(y2)
    warning('Polygon points out of bounds?')
    if isempty(x1); x1 = 1; end
    if isempty(x2); x2 = length(xlin); end
    if isempty(y1); y1 = 1;   end
    if isempty(y2); y2 = length(ylin); end
end
bbox = [xlin(x1),ylin(y2);xlin(x2),ylin(y1)];
if isBox
    H = Hgrid(y1:y2,x1:x2);
else
    % Goal is to get logical matrix with size of H with ones if the XY loc
    % is inside the polygon and zeros everywhere else
    in_H = zeros(length(ylin),length(xlin)); % Create empty logical matrix
    Xpoly2 = (Xpoly-xlin(1))/xstp+1; % Shift x values to use for poly2mask
    Ypoly2 = (ylin(1)-Ypoly)/ystp+1; % Shift y values to use for poly2mask
    % Find NaN value indices and loop through sections seperated by NaN
    nanIndices = find(isnan(Xpoly2));
    for i = 1:length(nanIndices)
        % Get index location of the subpart
        if i == 1  
            loc_part = 1:nanIndices(i)-1;
        else
            loc_part = nanIndices(i-1)+1:nanIndices(i)-1;
        end
        % Get sub part of polygon
        Xsub = Xpoly2(loc_part); Ysub = Ypoly2(loc_part);
        if length(Xsub)<=2 % To prevent errors if polygon is line
            Xsub = [Xsub(1); Xsub(2); Xsub(2); Xsub(1)];
            Ysub = [Ysub(1); Ysub(1); Ysub(2); Ysub(2)];
        end
        % Get matrix to find which points of H are in polygon subpart
        in_sub = poly2mask(Xsub,Ysub,length(ylin),length(xlin));
        in_H = in_H | in_sub;
    end
    Hnan=Hgrid; Hnan(~in_H) = NaN; 
    H = Hnan(y1:y2,x1:x2);
end
% If user wants to get X and Y matrix of similar shape as H
if nargout > 2
    [X,Y] = meshgrid(xlin(x1:x2),ylin(y1:y2));
end
end
