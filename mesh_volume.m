function volume = mesh_volume(V, F)
    % MESH_VOLUME computes the volume of a 3D mesh
    %
    % Inputs:
    %   V - #V by 3 matrix of vertex coordinates
    %   F - #F by 3 matrix of face indices into V
    %
    % Output:
    %   volume - scalar volume of the mesh

    % Initialize volume
    volume = 0;

    % Loop over all faces
    for i = 1:size(F, 1)
        % Get vertices of the current face
        v1 = V(F(i, 1), :);
        v2 = V(F(i, 2), :);
        v3 = V(F(i, 3), :);
        
        % Compute the volume of the tetrahedron formed by v1, v2, v3, and the origin
        % This is done using the scalar triple product
        tet_volume = dot(v1, cross(v2, v3)) / 6.0;
        
        % Accumulate the volume
        volume = volume + tet_volume;
    end

    % Take absolute value of the volume (since it can be negative)
    volume = abs(volume);
end
