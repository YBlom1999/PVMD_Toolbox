function [SM,V_sky,F_sky,azimuth,zenith,A_sky] = BackwardTracer(vertices,faces,nCells,cellCorners,normalSolarCell,reflectivity,Scattering,RT)
%BackwardTracer performers a backward tracing, based on the method from A.
%Calcabrini
%
% This function calculates the sensitivity map of the PV system
%
% Parameters
% ----------
% vertices : double
%   array of all the vertices in the environment
% faces : double
%   array of all faces in the environment (specifies which vertices form a
%   face)
% nCells : double
%   number of cells in the module
% cellCorners: double
%   The corners of all cells
% normalSolarCell: double
%   The normal vector of the solar cell
% reflectivity: double
%   The reflection (albedo) of every face
% Scattering: double
%   Indication whether the face has a specular reflection or a scattered
%   reflection
% RT: struct
%   The optical properties of the solar cell
%
% Returns
% -------
% SM : double
%   Sensitivity map for all sky elements
% V_sky : double
%   Vertices for all sky elements
% F_sky: double
%   Faces for all sky elements
% azimuth: double
%   The azimuth of every sky element
% zenith: double
%   The zenith angle of every sky element
% A_sky: double
%   The area of each sky element [Solid Angle]
%
% Developed by A. Calcabrini, adjusted by Y. Blom

EQE = RT.RAT;
EQE_wav = RT.wav;
EQE_angles = RT.aoi;
nLayers = size(RT.lay,1)-1;

N_faces = size(faces,1);%number of triangles in the scene
coor_faces = single(reshape(vertices(faces',:)',3,3,N_faces));
%rows correspond to the x,y,z coordinates of each vertex, columns
%correspond to the each of the 3 vertices of each triangle
norm_faces = my_cross(coor_faces(:,2,:)-coor_faces(:,1,:),...
    coor_faces(:,3,:)-coor_faces(:,1,:));
vert1 = permute(coor_faces(:,1,:),[1,3,2]); %Coordinates of the first vertex of each triangle
vert2 = permute(coor_faces(:,2,:),[1,3,2]); %Coordinates of the second vertex of each triangle
vert3 = permute(coor_faces(:,3,:),[1,3,2]); %Coordinates of the third vertex of each triangle

%% Define solar cell(s)

centroidSolarCell = squeeze(mean(cellCorners,2));%

%Use column vectors
if size(centroidSolarCell,2)==3
    centroidSolarCell=centroidSolarCell';
end
if size(normalSolarCell,2)==3
    normalSolarCell=normalSolarCell';
end


tiltsSolarCell = acosd(normalSolarCell(3,:));
azimsSolarCell = atan2d(normalSolarCell(1,:),normalSolarCell(2,:));azimsSolarCell = azimsSolarCell+(azimsSolarCell<0)*360; %degrees East from North

%% Sky discretization: Define the number of rays for the simulations (MODIFY AT YOUR CONVENIENCE)
N_refinement = 3;
% Make sky dome
[V_sky,F_sky,A_sky,~,zenith,azimuth] = icohemisphere(N_refinement,1);

% Generate primary rays
[~,~,~,r1,~,~] = icohemisphere(N_refinement,0);
if size(r1,2)==3; r1=r1'; end
r1_upper = r1(3,:)>0;

r1 = single(r1./sqrt(sum(r1.^2)));%normalize rays
% r1 are the primary rays referenced to the origin [0;0;0]
% omega_sp is the solid angle of each sky patch

% Generate secondary rays for second generation with possibly different angular resolution
[~,~,~,r2,~,~] = icohemisphere(N_refinement,0);
if size(r2,2)==3; r2=r2'; end
r2_upper = r2(3,:)>0;
r2 = single(r2./sqrt(sum(r2.*r2)));%normalize rays

%% Ray-trace
SM = zeros(size(F_sky,1),nCells,nLayers,length(EQE_wav));
for iCell=1:nCells
    cellTilt = tiltsSolarCell(iCell);%sensor tilt
    cellAzimuth = azimsSolarCell(iCell);%cell azimuth
    %Find that index (or indices) of the plane(s) that make the solar cell or irradiance sensor
    cellCenter = centroidSolarCell(:,iCell);%cell center coordinates
    
    %Check which faces are inside the cell
    iCellContainingPlane = nan(1,N_faces,'single');
    for iTriangle = 1:N_faces
        for iTriCorner = 1:3
            if sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,1,iCell)).^2)<1 || ...
                    sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,2,iCell)).^2)<1 || ...
                    sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,3,iCell)).^2)<1 || ...
                    sum((coor_faces(:,iTriCorner,iTriangle)-cellCorners(:,4,iCell)).^2)<1
                iCellContainingPlane(iTriangle) = iTriangle;
            end
        end
    end
    
    cellNormal = ...
        [sind(cellTilt)*sind(cellAzimuth);
        sind(cellTilt)*cosd(cellAzimuth);
        cosd(cellTilt)];%versor normal to the cell (y=N, x=E)
    
    r1CellCosAoi = sum(r1.*cellNormal);%cosine of the angles between the normal to the cell and the direction of primary rays
    r1CellFwdFlag = r1CellCosAoi>0; %Create a flag for all rays that are in front of the cell (i.e. cos(AOI) > 0)
    nr1Fwd = sum(r1CellFwdFlag);%Number of primary rays to be cast
    r1CellFwd = r1(:,r1CellFwdFlag);%Primary rays towards the front of the solar cell
    r1CellFwdCosAoi = r1CellCosAoi(r1CellFwdFlag);
    
    % select data type based on size Scene
    if N_faces<= 2^16-1
        planeIndexDataType = 'uint16';
    elseif N_faces<= 2^32-1
        planeIndexDataType = 'uint32';
    elseif N_faces<= 2^64-1
        planeIndexDataType = 'uint64';
    else
        error('Too many triangles in the scene. Simplify your 3D model.')
    end
    
    %Calculate which rays have an intersection
    [intxnPtR1,iIntxnPlaneR1] = solveTrianglesRaysIntersections(vert1,vert2,vert3,cellCenter,cellCenter+r1CellFwd,iCellContainingPlane);
    iIntxnPlaneR1 = cast(iIntxnPlaneR1,planeIndexDataType);%to save memory. Here NaN values are automatically casted to 0 which in MATLAB is an invalid index value
    r1EscFlag = isnan(sum(intxnPtR1,1));%Flag that indicates which rays did not intersect with the scene
    
    indexInt = 1:length(r1EscFlag);
    indexInt = indexInt(~r1EscFlag);
    nInt1 = nr1Fwd-sum(r1EscFlag);%Number of intersections between primary rays and the scene
    r1SkyFlag = r1EscFlag & r1CellFwd(3,:)>0;  %Indicates which rays hit the sky
    %It is possible that some rays that point down do not hit anything if the ground is not sufficiently large
    if sum(r1EscFlag)-sum(r1SkyFlag) > 0
        warning([num2str(sum(r1EscFlag)-sum(r1SkyFlag)),' out of ',num2str(length(r1EscFlag)),' rays are lost due to the finite ground extension.']);
    end
    r1EscFlag_full = zeros(size(r1CellFwdFlag));
    r1EscFlag_full(r1CellFwdFlag) = r1EscFlag;
    [~,ind_AOI] = min(abs(acosd(r1CellCosAoi)-EQE_angles));
    
    EQE_full = squeeze(EQE(:,ind_AOI,:));
    EQE_full(:,:,1) = (1-EQE(:,ind_AOI,1)-EQE(:,ind_AOI,end));
    sensitivity_fullSphere = r1CellCosAoi.*EQE_full(:,:,1:nLayers);
    sensitivity_fullSphere(:,~r1EscFlag_full,:) = 0;
    
    SM(:,iCell,:,:) = permute(sensitivity_fullSphere(:,r1_upper,:),[2,3,1]);
    
    
    intxnPtR1 = intxnPtR1(:,~r1EscFlag);%Keep only the intersection points and filter out NaN values
    iIntxnPlaneR1 = iIntxnPlaneR1(~r1EscFlag);%Filter out NaN values
    intxnPlaneNormal = squeeze(norm_faces(:,1,iIntxnPlaneR1));%Normal to the intersected triangles
    intxnPlaneNormal = intxnPlaneNormal./sqrt(sum(intxnPlaneNormal.*intxnPlaneNormal));%Normalize
    %Some of the normals may pointing away from the cell depending on the
    %order of the vertices of the triangle. Those normals must be flipped
    cosAoiR1IntxnPlaneSkew = sum(-r1CellFwd(:,~r1EscFlag).*intxnPlaneNormal);%the minus sign is needed because we want the vector that points TO the cell
    intxnPlaneNormal(:,cosAoiR1IntxnPlaneSkew<0) = -intxnPlaneNormal(:,cosAoiR1IntxnPlaneSkew<0);%Flip the normals that point away from the cell
    
    % 3) Secondary ray tracing from the intersection points
    % Differentiate between lambertian and specular reflector    
    for iInt1 = 1:nInt1%for every intersection between a primary ray and the scene
        iInt1ContainingPlane = iIntxnPlaneR1(iInt1);%Index of the plane contating the primary intersection
        albedo = reflectivity(iIntxnPlaneR1(iInt1),:);
        r1cosAOI = r1CellFwdCosAoi(indexInt(iInt1));
        [~,ind_AOI] = min(abs(acosd(r1cosAOI)-EQE_angles));
        if Scattering(iIntxnPlaneR1(iInt1))%Secondary rays are cast only for diffuse surfaces
            r2WallCosAoi = sum(r2.*intxnPlaneNormal(:,iInt1));%angle between the secondary rays and the normal to the wall
            r2WallFwdFlag = r2WallCosAoi>0;
            r2WallFwd = r2(:,r2WallCosAoi>0);%keep only the ray that go forward from the lambertial reflector
            r2Origin = intxnPtR1(:,iInt1);
            
            %Perform hemispherical sampling to calculate diffuse irradiance
            [intxnPt,~] = solveTrianglesRaysIntersections(vert1,vert2,vert3,r2Origin,r2Origin+r2WallFwd,iInt1ContainingPlane);
            r2EscFlag = isnan(sum(intxnPt,1));%Secondary rays that do not intersect (i.e., escape from) with the scene
            
            r2EscFlag_full = zeros(size(r2WallFwdFlag));
            r2EscFlag_full(r2WallFwdFlag) = r2EscFlag;
            
            EQE_full = squeeze(EQE(:,ind_AOI,:));
            EQE_full(:,1) = (1-EQE(:,ind_AOI,1)-EQE(:,ind_AOI,end));
            
            sensitivity_fullSphere = ones(length(albedo),nLayers,length(r2))*r1cosAOI.*EQE_full(:,1:nLayers);
            sensitivity_fullSphere(:,:,~r2EscFlag_full) = 0;
            
            SM(:,iCell,:,:) = squeeze(SM(:,iCell,:,:))+permute(albedo'.*sensitivity_fullSphere(:,:,r2_upper)/sum(r2WallFwdFlag),[3,2,1]);
            
        else %Ideally specular reflections
            %Find the intersecting point on the mirror and the normal of the mirror
            rsp1 = cellCenter - intxnPtR1(:,iInt1);% Ray from the cell to specular intersection point
            rsp1 = rsp1./sqrt(sum(rsp1.^2));%normalize
            specNormal = intxnPlaneNormal(:,iInt1);
            rsp2 = 2*sum(rsp1.*specNormal).*specNormal-rsp1; %rsp2 is rsp1 reflected on the mirror (referenced to [0 0 0])
            %Only if rsp2 points to the sky we check if it can reach the sky
            if rsp2(3)>0 %If the ray points down we already know it won't reach the sky so there's no need to trace it
                [intp,~] = solveTrianglesRaysIntersections(vert1,vert2,vert3,intxnPtR1(:,iInt1),intxnPtR1(:,iInt1)+rsp2,iInt1ContainingPlane);
                %if the ray reaches the sky then we look for the closest
                %sky patch index using the PRIMARY sky discretization
                if sum(isnan(intp))==length(intp)%if the ray didn't hit any plane
                    x = rsp2'*r2;%scalar product between the reflected ray and the sky rays
                    [~,iNearest] = max(x);
                    Nearest = zeros(length(r2),1);
                    Nearest(iNearest) = 1;
                    Nearest = Nearest(r2_upper);
                    [~,iNearest] = max(Nearest);
                    
                    EQE_full = squeeze(EQE(:,ind_AOI,:));
                    EQE_full(:,1) = (1-EQE(:,ind_AOI,1)-EQE(:,ind_AOI,end));
                    
                    SM(iNearest,iCell,:,:) = squeeze(SM(iNearest,iCell,:,:))+albedo.*r1cosAOI.*EQE_full(:,1:nLayers)';
                    
                end
            end
        end
    end
end


function [int_point,tri_ix,min_dist] = solveTrianglesRaysIntersections(vert1Triangle,vert2Triangle,vert3Triangle,rayOrg,rayEnd,iInPlane)
% solveTrianglesRaysIntersections calculates the closest intersecting points 
% beteween many unitary rays % and a set of triangles. The function returns 
% the intersecting % points and the correspondix triangle indeces.
% If there is no intersection the function returns NaN.
% the vertices of the triangles must be column vectors. tri_v1, tri_v2 and
% tri_v3 are the first, second and third vertices of the triangle. The
% dimensions of these matrices must be 3 by N, where N is the number of
% triangles in the set.
% iOnPlane can be used to indicate the index of triangle(s) that contain
% the ray origin point. These triangles will be ignored in the ray-plane
% colission detection. To leave the containing plane unspecified use NaN.
% https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
%
% This function calculates the sensitivity map of the PV system
%
% Parameters
% ----------
% vert1Triangle : double
%   The first vertex of every triangle
% vert2Triangle : double
%   The second vertex of every triangle
% vert3Triangle : double
%   The third vertex of every triangle
% rayOrg : double
%   The origin of every ray
% rayEnd : double
%   The end point of every ray
% iInPlane : double
%   The index of which triangles are inside the cell
%
% Returns
% -------
% int_point : double
%   The coordinates of intersection
% tri_ix : double
%   The triangle (face) of intersection
% min_dist : double
%   The distance to the intersection point
%   
%    
% Example of use:
% Calculate the intersection between a ray in the y direction and a
% triangle with vertices (1;4;-1), (0;4;1) and (-1;4;-1)
%
% l_org = [0;0;0];
% l_end = [0;1;0];
% p1 = [1;4;-1];
% p2 = [0;4;1];
% p3 = [-1;4;-1];
% [hit_point,tri_ix,d] = solveTrianglesRaysIntersections(p1,p2,p3,l_org,l_end);
% patch('Faces',[1 2 3],'Vertices',[P1';P2';P3'],'FaceColor','green');
% view(50,25)
% hold on
% plot3(hit_point(1),hit_point(2),hit_point(3),'kx','MarkerSize',10,'LineWidth',2)
% disp(['Distance to intersecting point: ',num2str(d)])
% disp(['Index of intersected triangle: ',num2str(tri_ix)])
%
% Developed by A. Calcabrini, adjusted by Y. Blom

nRays = size(rayEnd,2);
tri_ix = nan(1,nRays,'single');
min_dist = nan(1,nRays,'single');
int_point = nan(3,nRays,'single');
u = vert2Triangle-vert1Triangle;
v = vert3Triangle-vert1Triangle;
npvector = my_cross(u,v);%normals to the planes (not normalized)
np_squared = sum(npvector.^2);
npversor = npvector./sqrt(np_squared);


d = rayEnd-rayOrg; %Careful: d should be normalized! Check your input arguments.

for iRay = 1:nRays
    dAux = d(:,iRay);
    distToPlane = sum((vert1Triangle-rayOrg).*npversor,1)./(dAux'*npversor);
    distToPlane(distToPlane<0) = NaN;
    projectedPoint = distToPlane.*dAux+rayOrg;%calculate points on the planes
    
    w = projectedPoint-vert1Triangle;
    
    gamma = sum(my_cross(u,w).*npvector)./np_squared;
    beta = sum(my_cross(w,v).*npvector)./np_squared;
    alpha = 1-gamma-beta;
    
    %Find if the intersecting point between the plan and the ray is insdide the triangle
    isInsideTriangle = alpha>=0 & alpha<=1 & beta>=0 & beta<=1 & gamma>=0 & gamma<=1;
    distToPlane(~isInsideTriangle) = NaN; %Remove distances that are not inside the triangle
    distToPlane(~isnan(iInPlane)) = NaN; %Remove distances which are in the plane
    
    %Default return values
    int_pointAux = nan(3,1,'single'); %intersection
    min_distAux = NaN; %distance from ray origin to intersection
    tri_ixAux = NaN; %index (column) of the intersected plane
    
    if sum(~isnan(distToPlane)) %if at least one intersection is found
        [min_distAux,tri_ixAux] = min(distToPlane);
        if ~isnan(min_distAux)
            int_pointAux = projectedPoint(:,tri_ixAux);
        end
    end
    
    tri_ix(iRay) = tri_ixAux;
    min_dist(iRay) = min_distAux;
    int_point(:,iRay) = int_pointAux;
end
end



function c = my_cross(a,b)
% Cross product for column vectors. (much simpler and faster than cross() in MATLAB's library)
c = zeros(size(a),'single');
c(3,:) = a(1,:).*b(2,:)-a(2,:).*b(1,:);
c(1,:) = a(2,:).*b(3,:)-a(3,:).*b(2,:);
c(2,:) = a(3,:).*b(1,:)-a(1,:).*b(3,:);
end

function h = flatplot3(V,F,C,h)
%2D plot of icohemishpere. Facet color C.  3D cartesian coordinates are
%converted to 2D cylincer coordinates first.

if isnumeric(h)                     %if handle is just a number
    Vcyl = cart2cyl(V);             %convert 3D to 2D coordinates
    figure(h);
    clf
    h = patch('Vertices',Vcyl,'Faces',F,...
        'FaceVertexCData',C,'FaceColor','flat');    %plot in 2D coordinates
    %---make it pretty---
    axis equal off
    shading flat
    colormap(parula(512))
    caxis([-0.1,1])
    hc = colorbar;
    set(get(hc,'Title'),'string','Sensitivity [-]')
    text(  0, 95,'North','HorizontalAlignment','center')
    text( 95,  0,'East' ,'HorizontalAlignment','center','Rotation',-90)
    text(  0,-95,'South','HorizontalAlignment','center')
    text(-95,  0,'West' ,'HorizontalAlignment','center','Rotation',90)
else                                %if h is a real handle of existing figure
    set(h,'FaceVertexCData',C);     %just update the color
end
drawnow
%..........................................................................
    function Vcyl = cart2cyl(Vcart)
        %convert 3D cartesian to 2D cylinder coordinates
        %used to convert icohemisphere vertex coordinates
        %(not the light source, which is at the center of each vertex)
        zeni = atand(sqrt(Vcart(:,1).^2 + Vcart(:,2).^2)./Vcart(:,3));
        azi = atan2d(-Vcart(:,1),-Vcart(:,2));
        
        Vcyl_x = zeni .* -sind(azi);
        Vcyl_y = zeni .* -cosd(azi);
        Vcyl_z = zeros(size(Vcyl_x));               %z-coordinate is set to 0
        Vcyl = [Vcyl_x,Vcyl_y,Vcyl_z];              %combine
    end

end


end