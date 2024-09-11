% Load a mesh
[V, F] = load_mesh('bear.off');  

% 1. Compute Mean Curvature and Normal Vector

% 1.1 Compute Cotangent Laplacian (L) and Mass Matrix (M)
L = cotmatrix(V, F);
M = massmatrix(V, F, 'barycentric');

% 1.2 Compute Mean Curvature Vector (H)
H = -M\ (L * V);

% 1.3 Compute the Magnitude of Mean Curvature (mean curvature at each vertex)
mean_curvature = sqrt(sum(H.^2, 2));

% 1.4 Compute the Normal Vector
normals = per_vertex_normals(V, F); 

% Normalize the mean curvature and newthing terms to avoid large variations
mean_curvature = mean_curvature / max(mean_curvature);  % Normalize to max value

% 2. Compute Gaussian Curvature (using Angle Deficit Method)
vertex_areas = diag(M);  % Extract the diagonal elements to get vertex areas

% 2.1 Compute Angles in Each Triangle
% Internal angles at each vertex of each triangle in the mesh are computed
angles = internalangles(V, F);  % Compute the internal angles of each triangle

% 2.2 Sum the Angles at Each Vertex
% For each vertex, sum the angles of all the triangles that share this vertex
vertex_angles = accumarray(F(:), angles(:), [size(V, 1), 1]);  % Accumulate angles at each vertex

% 2.3 Compute Gaussian Curvature
% Gaussian curvature is computed as the angle deficit (2π minus the sum of angles) normalized by the vertex area
gaussian_curvature = (2 * pi - vertex_angles) ./ vertex_areas;  % Calculate Gaussian curvature at each vertex

% 3. Compute Laplacian of Mean Curvature
laplacian_H =  L * mean_curvature;  

% 4. Compute the Gradient of Mean Curvature
% Use gptoolbox's gradient operator for scalar fields
grad_H = ehsan_vertex_gradients(H,F);
mag_grad_H = sqrt(sum(grad_H.^2, 2));  % Magnitude of the gradient

% 5. Compute the new term: "newthing"
H_squared = mean_curvature .^ 2;
newthing = laplacian_H + ((H_squared - 2 * gaussian_curvature) .* mean_curvature) .* normals + mean_curvature .* grad_H;
newthing = newthing / max(max(abs(newthing)));  % Normalize to avoid large values

% 6. Define small arc equation
h_values = linspace(-0.5, 0.5, 100);  % Smaller step h
num_vertices = size(V, 1);

% 7. Visualization of the original mesh
figure;
% Make the figure full-screen
set(gcf, 'Position', get(0, 'ScreenSize'));
quiver3(V(:,1), V(:,2), V(:,3), normals(:,1), normals(:,2), normals(:,3), 2, 'g');
hold on;
axis equal;
view(3);
title('Mesh with small arcs at vertices');
xlabel('X'); ylabel('Y'); zlabel('Z');
colormap jet;
lighting phong;
camlight;


% 8. Plot small arcs at each vertex
for i = 1:num_vertices
    % Current vertex position
    vertex = V(i, :);
    % Normal vector at vertex
    normal = normals(i, :);
    % New term at vertex
    new_term = newthing(i, :);
    
    % Initialize arc points
    arc_points = zeros(length(h_values), 3);
    
    % Calculate arc for each h value
    for j = 1:length(h_values)
        h = h_values(j);
        % Small arc equation
        arc_points(j, :) = vertex + mean_curvature(i) * normal * h + 0.5 * new_term * h^2;
    end
    
    % Plot the arc as a curve
    plot3(arc_points(:,1), arc_points(:,2), arc_points(:,3), 'r-', 'LineWidth', 1.5);
end
trisurf(F, V(:,1), V(:,2), V(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% 9. Display the original mesh with small arcs
% Set zoom mode and allow user interaction
zoom on; % Enables interactive zoom
disp('Click on a region to zoom into it.');
hold off;
