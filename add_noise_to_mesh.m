function [V_noisy] = add_noise_to_mesh(V, noise_level)
    % Input:
    % V - Original vertices of the mesh (Nx3 matrix)
    % noise_level - The magnitude of noise to be added to each vertex
    
    % Generate random noise
    noise = noise_level * randn(size(V));
    
    % Add noise to the vertices
    V_noisy = V + noise;
end

