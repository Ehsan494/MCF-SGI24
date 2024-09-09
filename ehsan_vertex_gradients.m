function grad_at_vertices = ehsan_vertex_gradients(V, F)
    % COMPUTE_VERTEX_GRADIENTS Compute the gradient at every vertex in a triangular mesh
    %
    % grad_at_vertices = compute_vertex_gradients(V, F)
    %
    % Inputs:
    %   V  #vertices by 3 list of mesh vertex positions (x, y, z)
    %   F  #faces by 3 list of mesh face indices
    % Outputs:
    %   grad_at_vertices  #vertices by 3 list of gradient values (x, y, z)
    
    num_vertices = size(V, 1);
    num_faces = size(F, 1);
    
    % Initialize gradient matrix
    grad_at_vertices = zeros(num_vertices, 3);
    area_sum = zeros(num_vertices, 1);
    
    % Compute the gradient for each face
    for i = 1:num_faces
        % Indices of the vertices of the face
        v1 = F(i, 1);
        v2 = F(i, 2);
        v3 = F(i, 3);
        
        % Coordinates of the vertices
        p1 = V(v1, :);
        p2 = V(v2, :);
        p3 = V(v3, :);
        
        % Compute the normal vector of the face
        edge1 = p2 - p1;
        edge2 = p3 - p1;
        normal = cross(edge1, edge2);
        area = norm(normal) / 2;
        normal = normal / norm(normal); % Unit normal
        
        % Compute gradients for the face
        grad_v1 = (p2 - p3) / (2 * area);
        grad_v2 = (p3 - p1) / (2 * area);
        grad_v3 = (p1 - p2) / (2 * area);
        
        % Update gradient at vertices
        grad_at_vertices(v1, :) = grad_at_vertices(v1, :) + grad_v1;
        grad_at_vertices(v2, :) = grad_at_vertices(v2, :) + grad_v2;
        grad_at_vertices(v3, :) = grad_at_vertices(v3, :) + grad_v3;
        
        % Accumulate the area for normalization
        area_sum(v1) = area_sum(v1) + area;
        area_sum(v2) = area_sum(v2) + area;
        area_sum(v3) = area_sum(v3) + area;
    end
    
    % Normalize gradients by the area sum
    for i = 1:num_vertices
        if area_sum(i) > 0
            grad_at_vertices(i, :) = grad_at_vertices(i, :) / area_sum(i);
        end
    end
end
