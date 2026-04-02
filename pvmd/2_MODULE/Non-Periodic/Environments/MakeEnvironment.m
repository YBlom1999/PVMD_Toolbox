clear all;
close all;

load('Example_house.mat')
index = 1:16;
F_new = F(index,:);
V_new = V(1:max(max(F_new)),:);
rgb_new = rgb(index,:);
opa_new = opa(index);
Scattering_new = Scattering(index);
materials_new = materials(index);
Albedo_new = Albedo(index);

new_triangles = {[[-1, -3, 3];[-1, -3, 5];[-1, 0, 5]],  [[1, -3, 3];[1, -3, 5];[1, 0, 5]], [[-1, -3, 3];[1, -3, 3];[1, -3, 5]],  [[-1, -3, 3];[-1, -3, 5];[1, -3, 5]],  [[-1, 0, 5];[-1, -3, 5];[1, -3, 5]],  [[-1, 0, 5];[1, 0, 5];[1, -3, 5]]};

for i = 1:length(new_triangles)
F_new(size(F_new,1)+1,:) = size(V_new,1)+(1:3);
V_new(size(V_new,1)+(1:3),:) = new_triangles{i};
rgb_new(size(rgb_new,1)+1,:) = rgb_new(index(end),:);
opa_new(size(opa_new,1)+1) = opa_new(index(end));
materials_new{length(materials_new)+1} = materials_new{index(end)};
Albedo_new(size(Albedo_new,1)+1) = Albedo_new(index(end));
Scattering_new(size(Scattering_new,1)+1) = Scattering_new(index(end));

end

p = patch('Faces',F_new,'Vertices',V_new); 
set(p,'FaceVertexCData',rgb_new,'CDataMapping','scaled',...
    'FaceColor','flat','FaceVertexAlphaData',opa_new,...
    'AlphaDataMapping','none','FaceAlpha','flat')

axis equal off tight                    %xyz aspect ratio equal
view(30,30)                            %standard view from side
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-10,10]);
ylim([-10,10]);
zlim([0,10]);

%% Make new environment
V = V_new;
F = F_new;
Scattering = Scattering_new;
Albedo = Albedo_new;
opa = opa_new;
rgb = rgb_new;
materials = materials_new;
save('Example_house_new.mat', 'Albedo', 'F', 'Scattering', 'V', 'materials', 'opa', 'rgb')