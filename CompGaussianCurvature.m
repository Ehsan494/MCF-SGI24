%Angle Deficit Method: This algorithm calculates the Gaussian curvature 
% at a vertex by measuring how much the sum of the angles around that vertex deviates from 
% 2pi (the sum of angles in a flat plane). The Gaussian curvature is then obtained by normalizing 
% this angle deficit by the area associated with the vertex (usually the barycentric area).

% Output:
% 
% 1. Gaussian Curvature Vector:
%    - gaussian_curvature: A vector of size [N, 1], where N is the number of vertices in the mesh.
%      This vector contains the Gaussian curvature value at each vertex of the mesh.
% 
% 2. Visualization Plot:
%    - A 3D surface plot showing the Gaussian curvature of the mesh.
%      - The mesh surface is colored according to the Gaussian curvature values.
%      - A colorbar is included to indicate the range of curvature values.
%      - The plot includes axis labels for X, Y, and Z coordinates and is titled "Gaussian Curvature".

% Load the mesh
[V, F] = load_mesh('sphere.off');  % Replace 'sphere.off' with your actual file name

% Step 1: Compute Vertex Areas
% The mass matrix provides the area associated with each vertex using barycentric areas
M = massmatrix(V, F, 'barycentric');  % Compute the barycentric area mass matrix
vertex_areas = diag(M);  % Extract the diagonal elements to get vertex areas

% Step 2: Compute Angles in Each Triangle
% Internal angles at each vertex of each triangle in the mesh are computed
angles = internalangles(V, F);  % Compute the internal angles of each triangle

% Step 3: Sum the Angles at Each Vertex
% For each vertex, sum the angles of all the triangles that share this vertex
vertex_angles = accumarray(F(:), angles(:), [size(V, 1), 1]);  % Accumulate angles at each vertex

% Step 4: Compute Gaussian Curvature
% Gaussian curvature is computed as the angle deficit (2π minus the sum of angles) normalized by the vertex area
gaussian_curvature = (2 * pi - vertex_angles) ./ vertex_areas;  % Calculate Gaussian curvature at each vertex

% Step 5: Visualization of Gaussian Curvature
% Plot the Gaussian curvature on the mesh using color to indicate curvature values
figure;
trisurf(F, V(:, 1), V(:, 2), V(:, 3), gaussian_curvature, 'EdgeColor', 'none');  % Plot mesh with Gaussian curvature
colorbar;  % Display a colorbar to show the scale of Gaussian curvature
title('Gaussian Curvature');  % Title of the plot
axis equal;  % Ensure equal scaling along all axes
xlabel('X');  % Label for the X-axis
ylabel('Y');  % Label for the Y-axis
zlabel('Z');  % Label for the Z-axis
