% Script to evaluate and compare Mean Curvature Flow (MCF) methods:
% Explicit, and Desbrun et al's Semi-Implicit on a 3D mesh.
% This script calculates the smoothed mesh for different time steps and measures accuracy and computational efficiency.

% Load your target mesh
[V, F] = load_mesh('bear.off');

% Set the noise level (e.g., 0.01)
noise_level = 0.01;

% Add noise to the mesh
V_noisy = add_noise_to_mesh(V, noise_level);

% Define the number of iterations for smoothing
num_iterations = 100;
time_steps = [1e-3, 1e-4, 1e-5, 1e-6]; % Example time steps

% Initialize arrays to store errors, runtimes, and convergence rates
errors_explicit = zeros(length(time_steps), 1);
errors_semi_implicit = zeros(length(time_steps), 1);

runtime_explicit = zeros(length(time_steps), 1);
runtime_semi_implicit = zeros(length(time_steps), 1);

% Initialize convergence tracking
convergence_explicit = zeros(num_iterations, length(time_steps));
convergence_semi_implicit = zeros(num_iterations, length(time_steps));

% Loop over each time step to evaluate accuracy and efficiency
for idx = 1:length(time_steps)
    time_step = time_steps(idx);
    
    % Explicit Method
    V_explicit = V_noisy; % Start with noisy mesh
    tic; % Start timing
    for iter = 1:num_iterations
        L = cotmatrix(V_explicit, F);
        M = massmatrix(V_explicit, F, 'voronoi');
        epsilon = 0.001; % Small positive value for regularization
        M = M + epsilon * speye(size(M)); % Apply regularization
        HN = M \ (L * V_explicit);  % More stable than inv(M)*

        % Update vertex positions
        V_explicit = V_explicit + time_step * HN;
        
        % Track convergence rate
        convergence_explicit(iter, idx) = norm(V_explicit - V, 'fro');
    end
    runtime_explicit(idx) = toc; % Stop timing
    
    % Compute error for Explicit Method
    errors_explicit(idx) = convergence_explicit(end, idx);
    
    % Desbrun et al Semi-Implicit Method
    V_semi_implicit = V_noisy; % Start with noisy mesh
    tic; % Start timing
    for iter = 1:num_iterations
        L = cotmatrix(V_semi_implicit, F);
        M = massmatrix(V_semi_implicit, F, 'voronoi');
        c = time_step * (inv(M) * L);
        d = speye(size(V_semi_implicit, 1));
        A = speye(size(c)) - c; % Construct the system matrix
        epsilon = 0.001; % Small positive value for regularization
        A = A + epsilon * speye(size(A)); % Apply regularization

        U=zeros(size(V_semi_implicit));
        for i=1:3
            B=V_semi_implicit(:,i);
            % Solve the system A * X = B using Biconjugate Gradient (BiCG)
            [X, ~] = bicg(A, B, 1e-6, 1000);
            U(:,i)=X;
        end
        V_semi_implicit=U; 
        % Track convergence rate
        convergence_semi_implicit(iter, idx) = norm(V_semi_implicit - V, 'fro');
    end
    runtime_semi_implicit(idx) = toc; % Stop timing
    
    % Compute error for Semi-Implicit Method
    errors_semi_implicit(idx) = convergence_semi_implicit(end, idx);
    
   
end

% Plot the relationship between time step and accuracy for each method
figure;
loglog(time_steps, errors_explicit, '-o', 'DisplayName', 'Explicit Method');
hold on;
loglog(time_steps, errors_semi_implicit, '-s', 'DisplayName', 'Semi-Implicit Method');
xlabel('Time Step');
ylabel('Accuracy (Frobenius Norm of Error)');
title('Accuracy vs. Time Step for Different MCF Methods');
legend('show');
grid on;

% Plot the runtime vs. time step for each method
figure;
loglog(time_steps, runtime_explicit, '-o', 'DisplayName', 'Explicit Method');
hold on;
loglog(time_steps, runtime_semi_implicit, '-s', 'DisplayName', 'Semi-Implicit Method');
xlabel('Time Step');
ylabel('Runtime (seconds)');
title('Runtime vs. Time Step for Different MCF Methods');
legend('show');
grid on;

% Plot convergence rate for each method
figure;
plot(1:num_iterations, convergence_explicit(:, 1), '-o', 'DisplayName', 'Explicit Method');
hold on;
plot(1:num_iterations, convergence_semi_implicit(:, 1), '-s', 'DisplayName', 'Semi-Implicit Method');
xlabel('Iterations');
ylabel('Frobenius Norm of Error');
title('Convergence Rate for Different MCF Methods (Time Step 1e-3)');
legend('show');
grid on;
