% Shrinking Sphere Simulation and Visualization
% This script simulates the shrinking of a sphere over time under mean
% curvature flow and visualizes the process. The radius of the sphere 
% decreases as time progresses, and the corresponding volume is computed. 
% The shrinking sphere is saved as an animated GIF.

% INPUT:
% R0 - Initial radius of the sphere
% num_steps - Number of time steps for the simulation
% dt - Time step size
% pause_time - Pause time between frames in the visualization (to control speed)
% gif_filename - Name of the output GIF file

% OUTPUT:
% A GIF showing the shrinking of the sphere over time, and two plots 
% displaying the evolution of the sphere's radius and volume over time.


% Parameters
R0 = 1;             % Initial radius of the sphere
num_steps = 500;    % Number of time steps
dt = 0.001;   % Time step size
pause_time = 0.05;
gif_filename = 'shrinking_sphere.gif'; % Name of the output GIF file

% Time array
t = linspace(0, (num_steps-1)*dt, num_steps);

% Initialize radius array
R = zeros(1, num_steps);
R(1) = R0;

% Create initial sphere for visualization
[XS, YS, ZS] = sphere(50); % 50 is the resolution of the sphere

% Open a figure for visualization
figure;
subplot(2, 2, [1 2]); % Use top half of the figure for 3D visualization
h = surf(R0*XS, R0*YS, R0*ZS); % Initial sphere surface
axis equal;
title('Shrinking Sphere');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
shading interp;
colormap jet;
lighting gouraud;
camlight headlight;

% Fix the axis limits to keep the frame constant
axis([-R0 R0 -R0 R0 -R0 R0]);

% Loop to calculate the radius over time and update the visualization
for step = 2:num_steps
    R(step) = sqrt(R0^2 - 4 * t(step)); % R(t) = sqrt(R0^2 - 4t)
    if R(step) <= 0
        R(step) = 0;
        break;
    end
   
    % Update the sphere visualization
    set(h, 'XData', R(step)*XS, 'YData', R(step)*YS, 'ZData', R(step)*ZS);
    title(sprintf('Shrinking Sphere - Time: %.4f', t(step)));
    drawnow;

    % Capture the frame for GIF
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    
    % Write to the GIF File
    if step == 2
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', pause_time);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', pause_time);
    end

    % Pause for a short time to slow down the visualization
    pause(pause_time);

    
end
set(h, 'XData', R(step)*XS, 'YData', R(step)*YS, 'ZData', R(step)*ZS);
    title(sprintf('Shrinking Sphere - Time: %.4f', t(step)));
    drawnow;
% Calculate the volume at each time step
V = (4/3) * pi * R.^3; % V(t) = (4/3) * pi * R^3(t)

% Plotting the radius over time
subplot(2, 2, 3);
plot(t, R, 'b-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Radius R(t)');
title('Radius of the Sphere Over Time');
grid on;

% Plotting the volume over time
subplot(2, 2, 4);
plot(t, V, 'r-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Volume V(t)');
title('Volume of the Sphere Over Time');
grid on;

% Capture the frame for GIF
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    
    % Write to the GIF File
    if step == 2
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 2);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 2);
    end
% Display final results
fprintf('Final radius: %.4f\n', R(end));
fprintf('Final volume: %.4f\n', V(end));

% If the sphere has completely shrunk, notify
if R(end) == 0
    disp('The sphere has shrunk to a point.');
end
