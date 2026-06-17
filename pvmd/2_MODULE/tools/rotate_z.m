function V = rotate_z(V0,ang)
%rotate vertex points around z-axis by angle ang
%Here positive angle is defined as clockwise (0 = South)
V(:,1) =  cosd(ang) * V0(:,1) + sind(ang) * V0(:,2);
V(:,2) = -sind(ang) * V0(:,1) + cosd(ang) * V0(:,2);
V(:,3) = V0(:,3);
end