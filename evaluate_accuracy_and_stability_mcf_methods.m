% Script to evaluate and compare Mean Curvature Flow (MCF) methods:
% Explicit and Desbrun et al.'s Semi-Implicit on a 3D mesh.
% This script calculates the smoothed mesh for different time steps and measures accuracy and computational efficiency.
% revised: Sept, 5 2024

% Load your target mesh
[V, F] = load_mesh('bear.off');

% Set the noise level (e.g., 0.01)
noise_level = 0.01;

% Add noise to the mesh
V_noisy = add_noise_to_mesh(V, noise_level);

epsilon = 0.001; % Small perturbation for numerical stability
error_threshold = 1;     % Manually set threshold for stability
% Define parameters
num_iterations = 50; % Number of iterations for smoothing
%time_steps = logspace(log10(eps), log10(2), 20); % Time steps from small to larger values
time_steps = eps:0.2:10; % Time steps from small to larger values


% Initialize storage for error metrics
num_time_steps = length(time_steps); % Number of time steps

% Initialize arrays to store errors, runtimes, and stability metrics
errors_explicit = zeros(num_time_steps, 1);
errors_semi_implicit = zeros(num_time_steps, 1);
runtime_explicit = zeros(num_time_steps, 1);
runtime_semi_implicit = zeros(num_time_steps, 1);
stability_explicit = zeros(num_time_steps, 1);
stability_semi_implicit = zeros(num_time_steps, 1);
% Initialize convergence tracking
convergence_explicit = zeros(num_iterations, length(time_steps));
convergence_semi_implicit = zeros(num_iterations, length(time_steps));

% Loop over each time step to evaluate accuracy and efficiency
for idx = 1:num_time_steps
    time_step = time_steps(idx);
    
    % Explicit Method
    V_explicit = V_noisy; % Start with noisy mesh
    tic; % Start timing
    for iter = 1:num_iterations
        L = cotmatrix(V_explicit, F);
        M = massmatrix(V_explicit, F, 'voronoi');
        M = M + epsilon * speye(size(M)); % Apply regularization
        HN = M \ (L * V_explicit);  % More stable than inv(M)*
        % Update vertex positions
        V_explicit = V_explicit + time_step * HN;
        % Track convergence rate
        convergence_explicit(iter, idx) = norm(V_explicit - V, 'fro');
        
        % Check for numerical stability
        if any(isnan(V_explicit(:))) || any(isinf(V_explicit(:)))
            stability_explicit(idx) = NaN;
            break;
        end
    end
    runtime_explicit(idx) = toc; % Stop timing
    
    % Compute error for Explicit Method
    errors_explicit(idx) = norm(V_explicit - V, 'fro');
    
    % Semi-Implicit Method
    V_semi_implicit = V_noisy; % Start with noisy mesh
    tic; % Start timing
    for iter = 1:num_iterations
        L = cotmatrix(V_semi_implicit, F);
        M = massmatrix(V_semi_implicit, F, 'voronoi');
        M = M + epsilon * speye(size(M)); % Apply regularization
        % Set up the linear system for implicit time integration
        A = M - time_step * L;
        b = M * V_semi_implicit;
    
    % Solve the linear system to find the new vertices
    new_vertices = A \ b;
     
        % Track convergence rate
        convergence_semi_implicit(iter, idx) = norm(V_semi_implicit - V, 'fro');
        
        % Check for numerical stability
        if any(isnan(V_semi_implicit(:))) || any(isinf(V_semi_implicit(:)))
            stability_semi_implicit(idx) = NaN;
            break;
        end
    end
    runtime_semi_implicit(idx) = toc; % Stop timing

    
    % Compute error for Semi-Implicit Method
    errors_semi_implicit(idx) = norm(V_semi_implicit - V, 'fro');
    
end

% Determine the maximum allowable time step for stability (Δt_max)
% By checking when the error exceeds the threshold
delta_t_max_explicit = time_steps(find(errors_explicit(:, end) <= error_threshold, 1, 'last'));
delta_t_max_semi_implicit = time_steps(find(errors_semi_implicit(:, end) <= error_threshold, 1, 'last'));

% Display Results
fprintf('Max Stable Time Step (Explicit Euler): %f\n', delta_t_max_explicit);
fprintf('Max Stable Time Step (Desbrun et al Semi-Implicit): %f\n', delta_t_max_semi_implicit);

% Plot stability vs. time step for each method
figure;
loglog(time_steps, errors_explicit(:, end), '-o', 'DisplayName', 'Explicit Euler');
hold on;
loglog(time_steps, errors_semi_implicit(:, end), '-s', 'DisplayName', 'Semi-Implicit');
xlabel('Time Step');
ylabel('Frobenius Norm of Error');
title('Stability of MCF Methods: Explicit Euler vs Desbrun et al Semi-Implicit');
legend('show');
grid on;

% Plot the relationship between time step and accuracy for each method
figure;
loglog(time_steps, errors_explicit, '-o', 'DisplayName', 'Explicit Method');
hold on;
loglog(time_steps, errors_semi_implicit, '-s', 'DisplayName', 'Desbrun et al Semi-Implicit Method');
xlabel('Time Step');
ylabel('Frobenius Norm of Error');
title('Frobenius Norm of Error vs. Time Step for Different MCF Methods');
legend('show');
grid on;

% Plot the runtime vs. time step for each method
figure;
loglog(time_steps, runtime_explicit, '-o', 'DisplayName', 'Explicit Method');
hold on;
loglog(time_steps, runtime_semi_implicit, '-s', 'DisplayName', 'Desbrun et al Semi-Implicit Method');
xlabel('Time Step');
ylabel('Runtime (seconds)');
title('Runtime vs. Time Step for Different MCF Methods');
legend('show');
grid on;


% Log-log plot to compare convergence rates
figure;
loglog(time_steps, errors_explicit, '-o', 'DisplayName', 'Explicit Method');
hold on;
loglog(time_steps, errors_semi_implicit, '-s', 'DisplayName', 'Desbrun et al Semi-Implicit Method');
xlabel('Time Step');
ylabel('Frobenius Norm of Error');
title('Convergence Rate Comparison for MCF Methods');
legend('show');
grid on;

% Fit lines to determine slopes and display the convergence rates
% Avoid using NaNs in the fitting process by removing any invalid data points
valid_idx_explicit = ~isnan(errors_explicit);
valid_idx_semi_implicit = ~isnan(errors_semi_implicit);

% Fit line for Euler Explicit Method
fit_explicit = polyfit(log(time_steps(valid_idx_explicit)), log(errors_explicit(valid_idx_explicit, end)), 1);
slope_explicit = exp(fit_explicit(1));

% Fit line for Desbrun et al. Semi-Implicit Method
fit_semi_implicit = polyfit(log(time_steps(valid_idx_semi_implicit)), log(errors_semi_implicit(valid_idx_semi_implicit, end)), 1);
slope_semi_implicit = exp(fit_semi_implicit(1));

% Display the slopes, which represent the convergence rates
fprintf('Convergence Rate for Explicit Method: %.2f\n', slope_explicit);
fprintf('Convergence Rate for Semi-Implicit Method: %.2f\n', slope_semi_implicit);

% (Optional) Display the fitted lines
plot(fit_explicit, log(1:num_iterations)', log(errors_explicit), 'k--');
plot(fit_semi_implicit, log(1:num_iterations)', log(errors_semi_implicit), 'r--');


