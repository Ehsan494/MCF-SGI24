% 1. Load a Mesh

% Load the mesh from the off file   
[V, F] = load_mesh('sphere.off');

% 2. Display the Mesh
figure;
trisurf(F, V(:,1), V(:,2), V(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
title('Original Mesh');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Define the number of iterations and timestep for smoothing
num_iterations = 200;
time_step = 0.01;

% 3. Mean Curvature Flow (MCF) - Explicit Method
V_explicit = V; % Copy the vertices for the explicit method

for iter = 1:num_iterations
    % Compute the Laplace-Beltrami operator and mean curvature normal
    L = cotmatrix(V_explicit, F);
    HN = -L * V_explicit;
    
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
    L = cotmatrix(V_implicit, F);
    
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
