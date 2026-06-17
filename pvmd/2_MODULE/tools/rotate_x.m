function V = rotate_x(V0,ang)
%rotate vertex points around x-axis by angle ang
V(:,1) = V0(:,1);
V(:,2) = cosd(ang) * V0(:,2) - sind(ang) * V0(:,3);
V(:,3) = sind(ang) * V0(:,2) + cosd(ang) * V0(:,3);
end