% This MATLAB script demonstrates the application of the Semi-Implicit
% Mean Curvature Flow (MCF) smoothing technique on a 3D mesh. The goal is to iteratively
% adjust the vertices of the mesh to smooth the surface by reducing noise
% and surface irregularities. The script also visualizes the process
% by creating a GIF of the smoothing iterations and plots the change in
% volume of the mesh over time.

% ---------- Load and Prepare the Mesh ----------

% Load your mesh from an OFF file (bear.off), which contains vertex and
% face data of a 3D bear model.
[V, F] = load_mesh('bear.off');

% Set the noise level to be added to the mesh vertices. This simulates
% imperfections or irregularities in the mesh. Example: noise_level = 0.01
noise_level = 0.01;

% Add random noise to the mesh vertices. This creates a "noisy" version
% of the mesh, simulating surface irregularities that need to be smoothed.
V_noisy = add_noise_to_mesh(V, noise_level);

% Display the noisy mesh to visualize the initial state before smoothing.
figure;
trisurf(F, V_noisy(:,1), V_noisy(:,2), V_noisy(:,3), 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
axis equal; % Ensure equal scaling on all axes
lighting gouraud; % Apply smooth lighting to the mesh
camlight; % Add a light source to improve visual appearance
title('Noisy Mesh'); % Title of the figure
xlabel('X'); ylabel('Y'); zlabel('Z'); % Label the axes

% ---------- Initialize Smoothing Parameters ----------

% Define the number of iterations for the smoothing process. This controls
% how many times the smoothing operation will be applied to the mesh.
num_iterations = 10;

% Define the timestep for each iteration. This controls the "amount" of
% smoothing applied in each iteration. Larger timesteps lead to more aggressive smoothing.
time_step = 0.1;

% Initialize an array to store the volume of the mesh at each iteration.
% This helps in tracking how the volume changes during the smoothing process.
volumes_semi = zeros(num_iterations+1, 1);

% Store the initial volume of the noisy mesh as the first entry in the array.
volumes_semi(1, 1) = mesh_volume(V_noisy, F);

% ---------- Semi-Implicit Method by Desbrun et al. ----------

% Initialize the mesh vertices for the semi-implicit smoothing method.
V_semi = V_noisy; % Start with the noisy mesh

% Set up a filename for the GIF that will store the animation of the
% smoothing process.
gif_filename_semi_implicit = 'mcf_semi_implicit_method.gif';

% Iterate through the smoothing process for the specified number of iterations.
for iter = 1:num_iterations
    % Compute the Laplace-Beltrami operator, which captures the curvature
    % information of the mesh and is used to determine how the vertices
    % should be moved to smooth the mesh.
    L = cotmatrix(V_semi, F);
    
    % Compute the mass matrix, which is used to normalize the curvature
    % information. This ensures that the smoothing is done proportionally
    % to the local area around each vertex.
    M = massmatrix(V_semi, F, 'barycentric');
    
    % Compute the coefficient matrix for the semi-implicit method. This matrix
    % is used to solve for the new vertex positions in each iteration.
    c = time_step * (inv(M) * L);
    
    % Initialize the identity matrix. This matrix will be used in constructing
    % the system of equations that needs to be solved for each vertex.
    d = speye(size(V_semi, 1));
    A = speye(size(c)) - c; % Construct the system matrix
    
    % Regularize the matrix A by adding a small positive value to its diagonal.
    % This ensures that the matrix is positive definite and avoids numerical
    % instability during the solution process.
    epsilon = 0.001; % Small positive value for regularization
    A = A + epsilon * speye(size(A)); % Apply regularization
    
    % Solve the linear system A * V_new = V_old for each coordinate (X, Y, Z).
    U = zeros(size(V_semi));
    for i = 1:3
        B = V_semi(:, i); % Extract the current coordinate vector (X, Y, or Z)
        
        % Solve the system using the Biconjugate Gradient (BiCG) method.
        % BiCG is an iterative solver suitable for large, sparse matrices.
        [X, ~] = bicg(A, B, 1e-6, 1000);
        U(:, i) = X; % Store the solution (updated coordinates)
    end
    V_semi = U; % Update the vertex positions for the next iteration

    % Compute the volume of the updated mesh and store it in the array.
    volumes_semi(iter + 1) = mesh_volume(V_semi, F);

    % Display the smoothed mesh at the current iteration.
    trisurf(F, V_semi(:,1), V_semi(:,2), V_semi(:,3), 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
    axis equal;
    lighting gouraud;
    camlight;
    title(sprintf('MCF - Semi-Implicit Method, Iteration %d', iter));
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % Save the current frame as an image file for creating the GIF.
    frame_filename = sprintf('semi_implicit_frame_%02d.png', iter);
    saveas(gcf, frame_filename);
    
    % Read the saved image and convert it to indexed color format.
    img = imread(frame_filename);
    [img_ind, cmap] = rgb2ind(img, 256);
    
    % Write the indexed image to the GIF file. Append each frame to create an animation.
    if iter == 1
        imwrite(img_ind, cmap, gif_filename_semi_implicit, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(img_ind, cmap, gif_filename_semi_implicit, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

% ---------- Plot the Change
