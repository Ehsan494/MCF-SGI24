% This script demonstrates the application of the forward euler
% explicit 
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
num_iterations = 400;

% Define the timestep for each iteration. This controls the "amount" of
% smoothing applied in each iteration. Larger timesteps lead to more aggressive smoothing.
time_step = 0.0001;

% Initialize an array to store the volume of the mesh at each iteration.
% This helps in tracking how the volume changes during the smoothing process.
volumes_explicit = zeros(num_iterations+1, 1);

% Store the initial volume of the noisy mesh as the first entry in the array.
volumes_explicit(1, 1) = mesh_volume(V_noisy, F);

% ---------- Euler-Explicit Method ----------

% Initialize the mesh vertices for the semi-implicit smoothing method.
V_explicit = V_noisy; % Start with the noisy mesh

% Set up a filename for the GIF that will store the animation of the
% smoothing process.
gif_filename_euler_explicit = 'mcf_euler_explicit_method.gif';

% Iterate through the smoothing process for the specified number of iterations.
for iter = 1:num_iterations

    L = cotmatrix(V_explicit, F);
    M = massmatrix(V_explicit, F, 'barycentric');
    epsilon = 0.001; % Small positive value for regularization
    M = M + epsilon * speye(size(M)); % Apply regularization
    HN = M \ (L * V_explicit);  % More stable than inv(M)*

    % Update vertex positions
    V_explicit = V_explicit + time_step * HN;

    % Compute the volume of the updated mesh and store it in the array.
    volumes_explicit(num_iterations+ 1) = mesh_volume(V_explicit, F);

    % Display the smoothed mesh at the current iteration.
    trisurf(F, V_explicit(:,1), V_explicit(:,2), V_explicit(:,3), 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
    axis equal;
    lighting gouraud;
    camlight;
    title(sprintf('MCF - Euler_explicit Method, Iteration %d', iter));
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % Save the current frame as an image file for creating the GIF.
    frame_filename = sprintf('euler_explicit_frame_%02d.png', iter);
    saveas(gcf, frame_filename);
    
    % Read the saved image and convert it to indexed color format.
    img = imread(frame_filename);
    [img_ind, cmap] = rgb2ind(img, 256);
    
    % Write the indexed image to the GIF file. Append each frame to create an animation.
    if iter == 1
        imwrite(img_ind, cmap, gif_filename_euler_explicit, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(img_ind, cmap, gif_filename_euler_explicit, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

% ---------- Plot the Change in Volume During Smoothing ----------

% Create a new figure to plot the change in mesh volume over the iterations.
figure;
hold on; % Hold on to allow multiple plots on the same figure

% Plot the volume of the mesh at each iteration for the euler explicit method.
plot(0:num_iterations, volumes_explicit, '-x', 'DisplayName', 'euler_forward_explicit Method');
xlabel('Iteration'); 
ylabel('Volume'); 
title('Change in Volume During Mesh Smoothing'); % Title of the plot
legend; 
grid on; 
hold off; 
