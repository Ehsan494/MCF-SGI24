%The purpose of the code is to demonstrate two different approaches for performing Mean Curvature Flow (MCF) on a 3D mesh. 
% Mean Curvature Flow is a geometric process where the vertices of the mesh are iteratively adjusted to smooth the surface. 
% This is done by moving each vertex in the direction of the mean curvature normal vector, effectively reducing surface irregularities over time.

% Input:
% - 'bear.off': A 3D mesh file representing the bear model.
% - noise_level: A parameter controlling the amount of noise added to the original mesh.
% - num_iterations: Number of iterations for the MCF process.
% - time_step: A small timestep used in the explicit and semi-implicit methods.

% Output:
% - A figure displaying the original noisy mesh, the result of explicit MCF, and the result of semi-implicit MCF.

[V, F] = load_mesh('bear.off');

% Set the noise level (e.g., 0.01)
noise_level = 0.01;

% Add noise to the mesh
V_noisy = add_noise_to_mesh(V, noise_level);

% Display the noisy mesh
figure;
trisurf(F, V_noisy(:,1), V_noisy(:,2), V_noisy(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
title('Noisy Mesh');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Define the number of iterations and timestep for smoothing
num_iterations = 400;
time_step = eps;

% 3. Mean Curvature Flow (MCF) - Explicit Method
V_explicit = V; % Copy the vertices for the explicit method

for iter = 1:num_iterations
    % Compute the Laplace-Beltrami operator and mean curvature normal
    L = cotmatrix(V_explicit, F);
    M = massmatrix(V, F, 'voronoi'); 
    %M = massmatrix(V, F, 'barycentric');
    HN = inv(M)* (L * V_explicit);
    
    % Update vertex positions
    V_explicit = V_explicit + time_step * HN;
end

% Display the smoothed mesh - Explicit Method
figure;
trisurf(F, V_explicit(:,1), V_explicit(:,2), V_explicit(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
title('MCF - Explicit Method');
xlabel('X');
ylabel('Y');
zlabel('Z');

% 3. Mean Curvature Flow (MCF) - Semi-Implicit Method
V_implicit = V; % Copy the vertices for the implicit method

for iter = 1:num_iterations
    % Compute the Laplace-Beltrami operator
     L = cotmatrix(V_explicit, F);
     M = massmatrix(V, F, 'voronoi');
     HN = inv(M)* (L * V_explicit);
    
    % Solve (I - time_step * L) * V_new = V_old
    V_implicit = (speye(size(V_implicit, 1)) - time_step * L) \ V_implicit;
end

% Display the smoothed mesh - Semi-Implicit Method
figure;
trisurf(F, V_implicit(:,1), V_implicit(:,2), V_implicit(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
title('MCF - Semi-Implicit Method');
xlabel('X');
ylabel('Y');
zlabel('Z');
