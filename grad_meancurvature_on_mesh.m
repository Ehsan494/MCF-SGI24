% Step 1: Load a Mesh
[V, F] = load_mesh('sphere.off');

% Step 2: Compute Laplace-Beltrami Operator (L) and Mass Matrix (M)
L = cotmatrix(V, F);  % cotangent Laplace-Beltrami operator
M = massmatrix(V, F, 'barycentric');

% Step 3: Compute the Magnitude of Mean Curvature (mean curvature at each vertex)
H =M\ (L * V);
mean_curvature = sqrt(sum(H.^2, 2));

% Step 5: Compute the Gradient of Mean Curvature
% Use gptoolbox's gradient operator for scalar fields
grad_H = ehsan_vertex_gradients(H, F);

Mag_grad_H= sqrt(sum(grad_H.^2, 2));

% Step 6: Visualization or Further Processing
figure;
trisurf(F, V(:,1), V(:,2), V(:,3), Mag_grad_H, 'EdgeColor', 'none');
axis equal;
lighting gouraud;
camlight;
colorbar;
title('Gradient of Mean Curvature');
xlabel('X');
ylabel('Y');
zlabel('Z');
