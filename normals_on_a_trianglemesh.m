% Load the Mesh
% Load the mesh from the off file   
[V, F] = load_mesh('sphere.off');

% Compute Per-Vertex Normals
N_vertex = per_vertex_normals(V, F);


% Display the Mesh with Vertex Normals
figure;
trisurf(F, V(:,1), V(:,2), V(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none');
hold on;
quiver3(V(:,1), V(:,2), V(:,3), N_vertex(:,1), N_vertex(:,2), N_vertex(:,3), 0.5, 'r');
axis equal;
lighting gouraud;
camlight;
title('Mesh with Per-Vertex Normals');
xlabel('X');
ylabel('Y');
zlabel('Z');
