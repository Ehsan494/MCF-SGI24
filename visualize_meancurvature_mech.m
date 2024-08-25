% 1. Load a Mesh
[V, F] = load_mesh('sphere.off');

% 2. Compute Cotangent Laplacian (L) and Mass Matrix (M)
L = cotmatrix(V, F);
M = massmatrix(V, F, 'barycentric');

% 3. Compute Mean Curvature Vector (H)
H = -inv(M) * (L * V);

% 4. Compute the Magnitude of Mean Curvature (mean curvature at each vertex)
mean_curvature = sqrt(sum(H.^2, 2));

% Optional: Visualize the Mean Curvature on the Mesh
trisurf(F, V(:,1), V(:,2), V(:,3), mean_curvature, 'EdgeColor', 'none');
axis equal;
colorbar;
title('Mean Curvature on Mesh');
xlabel('X');
ylabel('Y');
zlabel('Z');
