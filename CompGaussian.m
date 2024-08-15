% Load the mesh
[V, F] = load_mesh('your_mesh_file.off');  % Replace 'your_mesh_file.off' with your actual file name

% Step 1: Compute Vertex Areas
M = massmatrix(V, F, 'barycentric');  % Barycentric area mass matrix
vertex_areas = diag(M);  % Extract diagonal elements to get vertex areas

% Step 2: Compute Angles in Each Triangle
angles = internalangles(V, F);  % Compute angles at each vertex of each triangle

% Step 3: Sum the Angles at Each Vertex
vertex_angles = accumarray(F(:), angles(:), [size(V, 1), 1]);

% Step 4: Compute Gaussian Curvature
gaussian_curvature = (2 * pi - vertex_angles) ./ vertex_areas;

% Visualization of Gaussian Curvature
figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), gaussian_curvature, 'EdgeColor', 'none');
colorbar;
title('Gaussian Curvature');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
