% The purpose of the code is to demonstrate two different
% approaches (Euler Explicit + Desbrun's et al Semi-Implicit)
% for performing Mean Curvature Flow (MCF) on a 3D mesh. 
% Mean Curvature Flow is a geometric process where the vertices of the 
% mesh are iteratively adjusted to smooth the surface. 
% This is done by moving each vertex in the direction of the mean curvature normal vector,
% effectively reducing surface irregularities over time.

% Input:
% - 'bear.off': A 3D mesh file representing the bear model.
% - noise_level: A parameter controlling the amount of noise added to the original mesh.
% - num_iterations: Number of iterations for the MCF process.
% - time_step: A small timestep used in the explicit and semi-implicit methods.

% Output:
% - A figure displaying the original noisy mesh, the result of explicit MCF, and the result of semi-implicit MCF.

% Load your mesh
[V, F] = load_mesh('bear.off');

% Set the noise level (e.g., 0.01)
noise_level = 0.01;

% Regularization parameter
epsilon = eps;

% Add noise to the mesh
V_noisy = add_noise_to_mesh(V, noise_level);

% Display the noisy mesh
figure;
trisurf(F, V_noisy(:,1), V_noisy(:,2), V_noisy(:,3), 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
title('Noisy Mesh');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Define the number of iterations and timestep for smoothing
num_iterations = 100;
time_step = 0.0001;


% 3. Mean Curvature Flow (MCF) - Explicit Method
V_explicit = V_noisy; % Copy the vertices for the explicit method

for iter = 1:num_iterations
    % Compute the Laplace-Beltrami operator and mean curvature normal
    L = cotmatrix(V_explicit, F);
    M = massmatrix(V_explicit, F, 'barycentric'); 
    HN = inv(M)*(L * V_explicit);
    
    V_explicit = V_explicit + time_step * HN;
    
end

% Display the smoothed mesh - Explicit Method
figure;
trisurf(F, V_explicit(:,1), V_explicit(:,2), V_explicit(:,3), 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
title('MCF - Explicit Method');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Semi-Implicit Method by Desbrun et al.
V_semi = V_noisy; % Copy the vertices for the semi-implicit method

for iter = 1:num_iterations
    % Compute the Laplace-Beltrami operator
    L = cotmatrix(V_semi, F);
    M = massmatrix(V_semi, F, 'barycentric');
    c=time_step * ( inv(M)*L); d=speye(size(V_semi, 1));
    A = speye(size(c)) - c;
    % Regularize the matrix to ensure positive definiteness
    epsilon = eps; % Small positive value
    A = A + epsilon * speye(size(A)); % Regularize the matrix
        
    % Solve the system A * V_new = V_old
    
    %V_semi=A\ V_semi;
    U=zeros(size(V_semi));
    for i=1:3
        B=V_semi(:,i);
        % Solve the system A * X = B using Biconjugate Gradient (BiCG)
        [X, ~] = bicg(A, B, 1e-6, 1000);
        U(:,i)=X;
    end
    V_semi=U; 
      
        
end

% Display the smoothed mesh Semi-Implicit Method
figure;
trisurf(F, V_semi(:,1), V_semi(:,2), V_semi(:,3), 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
title('MCF - Semi-Implicit Method');
xlabel('X');
ylabel('Y');
zlabel('Z');
