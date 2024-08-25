% Load a mesh
[V, F] = load_mesh('sphere.off');  

% Compute Mean Curvature and Normal Vector
L = cotmatrix(V, F);
M = massmatrix(V, F, 'barycentric');
H = inv(M) * (L * V);
mean_curvature = sqrt(sum(H.^2, 2));
normals = per_vertex_normals(V, F); 

%gaussian_curvature = discrete_gaussian_curvature(V,F);
laplacian_H =  L * mean_curvature;  
grad_H = grad(V, F) * mean_curvature;

H_squared = mean_curvature .^ 2;
%grad_H=grad_H(1:299,1);

vertex_areas = diag(M);  % Extract the diagonal elements to get vertex areas

% Step 2: Compute Angles in Each Triangle
% Internal angles at each vertex of each triangle in the mesh are computed
angles = internalangles(V, F);  % Compute the internal angles of each triangle

% Step 3: Sum the Angles at Each Vertex
% For each vertex, sum the angles of all the triangles that share this vertex
vertex_angles = accumarray(F(:), angles(:), [size(V, 1), 1]);  % Accumulate angles at each vertex

% Step 4: Compute Gaussian Curvature
% Gaussian curvature is computed as the angle deficit (2π minus the sum of angles) normalized by the vertex area
gaussian_curvature = (2 * pi - vertex_angles) ./ vertex_areas;  % Calculate Gaussian curvature at each vertex


% Visualization
figure;

% Visualize the surface with color mapping based on mean curvature
trisurf(F, V(:,1), V(:,2), V(:,3), mean_curvature, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
hold on;

% Plot the normals as quiver plot
quiver3(V(:,1), V(:,2), V(:,3), normals(:,1), normals(:,2), normals(:,3), 0.5, 'g');

% Improve lighting and shading
light('Position',[-1 -1 0.5],'Style','infinite');
lighting gouraud;
material shiny;
colorbar; % Show the color scale

% Set axis properties
axis equal;
axis off;
title('Visualization of Normals and Guassian Mean Curvature Flow');
xlabel('X');
ylabel('Y');
zlabel('Z');
rotate3d on; % Enable rotation

% Optional: Add rotation animation
for angle = 0:1:360
    view(angle, 30);
    drawnow;
end

% Initialize the GIF
filename = 'mean_curvature_rotation.gif';

% Optional: Add rotation animation and save as GIF
for angle = 0:1:360
    view(angle, 30);
    drawnow;

    % Capture the plot as an image
    frame = getframe(gcf);
    img = frame2im(frame);
    [img_ind, colormap] = rgb2ind(img, 256);
    
    % Write to the GIF File
    if angle == 0
        imwrite(img_ind, colormap, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
    else
        imwrite(img_ind, colormap, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
end
