% Visualizing mc, per_vertex_normals, with color mapping based on mc


% Load a mesh
[V, F] = load_mesh('sphere.off');  

% Compute Mean Curvature and Normal Vector
L = cotmatrix(V, F);
M = massmatrix(V, F, 'barycentric');
H = -M \ (L * V);
mean_curvature = sqrt(sum(H.^2, 2));
normals = per_vertex_normals(V, F); 
laplacian_H =  L * mean_curvature;  
grad_H = grad(V, F) * mean_curvature;
H_squared = mean_curvature .^ 2;


% Initialize the figure
figure;
filename = 'Visualization_Curvature_Quantities.gif'; % Output GIF file name

% Visualize the surface with color mapping based on mean curvature
trisurf(F, V(:,1), V(:,2), V(:,3), mean_curvature, 'FaceAlpha', 0.9, 'EdgeColor', 'none');
colormap(jet);
colorbar;
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
title('Visualization of MC and Per Vertex Normals');
xlabel('X');
ylabel('Y');
zlabel('Z');
rotate3d on; % Enable rotation

% Capture frames for the GIF
for angle = 0:1:360
    view(angle, 30);
    drawnow;
    
    % Capture the current frame as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF file
    if angle == 0
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.03);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
    end
end
