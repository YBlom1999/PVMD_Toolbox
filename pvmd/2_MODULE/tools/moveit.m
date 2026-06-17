function V = moveit(V,xyz)
%translate (move) vertex point in xyz-direction
V(:,1) = V(:,1) + xyz(1);
V(:,2) = V(:,2) + xyz(2);
V(:,3) = V(:,3) + xyz(3);
end