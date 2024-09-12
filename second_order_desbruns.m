% The purpose of the code is to demonstrate 3 different approaches for performing Mean Curvature Flow (MCF) on a 3D mesh. 
% Mean Curvature Flow is a geometric process where the vertices of the mesh are iteratively adjusted to smooth the surface. 
% This is done by moving each vertex in the direction of the mean curvature normal vector, effectively reducing surface irregularities over time.

% Input:
% - 'bear.off': A 3D mesh file representing the bear model.
% - noise_level: A parameter controlling the amount of noise added to the original mesh.
% - num_iterations: Number of iterations for the MCF process.
% - time_step: A small timestep used in the explicit and semi-implicit methods.

% Output:
% - GIF files displaying the smoothing process for the explicit and semi-implicit MCF methods.

% Load your mesh
[V, F] = load_mesh('bear.off');

% Set the noise level (e.g., 0.01)
noise_level = 0.001;

% Add noise to the mesh
V_noisy = add_noise_to_mesh(V, noise_level);

% Expand limits slightly if needed

xLimits = xlim;
yLimits = ylim;
zLimits = zlim;

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
num_iterations = 50;
time_step = 0.00001;
% Initialize arrays to store volumes
volumes_semi = zeros(num_iterations+1, 1);
volumes_semi(1, 1) = mesh_volume(V_noisy,F);

% Set axis limits based on the noisy mesh
axis equal;


% Set camera view
view_angle = view;

% ---------- Semi-Implicit Method by Desbrun et al. ----------
V_semi = V_noisy; % Copy the vertices for the semi-implicit method

% Create GIF file for the semi-implicit method
gif_filename_semi_implicit = 'mcf_secondorder_semi_implicit_method.gif';

for iter = 1:num_iterations
    iter
    % Compute the Laplace-Beltrami operator
     L = cotmatrix(V_semi, F);
     M = massmatrix(V_semi, F, 'voronoi');
     LM=M\L;
     % First-order term (mean curvature flow step)
     first_order_term= time_step*LM;
     eigenvalues = eigs(first_order_term);

% Compute the spectral radius (largest absolute value of the eigenvalues)
spectral_radius = max(abs(eigenvalues));

% Display the spectral radius
disp(['Spectral Radius: ', num2str(spectral_radius)]);

     % Second-order term (refinement step based on second-order derivatives)
     second_order_term = time_step^2 * (LM) * (LM);
     A=speye(size(first_order_term));
    % Update vertices
    V_semi = (A + first_order_term + second_order_term)*V_semi;

        % Check for numerical stability
        if any(isnan(V_semi(:))) || any(isinf(V_semi(:)))
            stability_semi(idx) = NaN;
            break;
        end

    % Compute the volume of the updated mesh
    volumes_semi(num_iterations+1) = mesh_volume(V_semi, F);

    % Expand limits slightly if needed


    % Display the smoothed mesh - Semi-Implicit Method
    trisurf(F, V_semi(:,1), V_semi(:,2), V_semi(:,3), 'FaceColor', [1, 0, 1], 'EdgeColor', 'none');
    axis equal;
    lighting gouraud;
    camlight;
    title(sprintf('MCF - Second-order-Semi-Implicit Method, Iteration %d', iter));
    xlabel('X');
    ylabel('Y');
    zlabel('Z');

    
    % Save the current frame as an image
    frame_filename = sprintf('secondorder_semi_implicit_frame_%02d.png', iter);
    saveas(gcf, frame_filename);
    
    % Read the image and convert to indexed color
    img = imread(frame_filename);
    [img_ind, cmap] = rgb2ind(img, 256);
    
    % Write to GIF file
    if iter == 1
        imwrite(img_ind, cmap, gif_filename_semi_implicit, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(img_ind, cmap, gif_filename_semi_implicit, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
% Plot the volumes against the iterations
figure;
hold on;
plot(0:num_iterations, volumes_semi, '-x', 'DisplayName', 'Semi-Implicit Method');
xlabel('Iteration');
ylabel('Volume');
title('Change in Volume During Mesh Smoothing');
legend;
grid on;
hold off;


