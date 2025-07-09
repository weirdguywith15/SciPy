% Clear previous data
clear;
clc;

% Parameters for the neutron star
M = 1e30; % Mass of the neutron star (arbitrary units)
r0 = [0, 0]; % Position of the neutron star (arbitrary units)
sigma = 0.2; % Standard deviation for Gaussian (arbitrary units)

% Create a grid for the heat map
x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);

% Calculate the Gaussian distribution for the neutron star
Z = M * exp(-((X - r0(1)).^2 + (Y - r0(2)).^2) / (2 * sigma^2));

% Plotting the 3D heat map
figure;
surf(X, Y, Z);
shading interp; % Smooth the surface
colorbar; % Add color bar to indicate density values
xlabel('X Position (arbitrary units)');
ylabel('Y Position (arbitrary units)');
zlabel('Density (arbitrary units)');
title('3D Heat Map of a Single Neutron Star');
view(3); % Set view to 3D
% Clear previous data
clear;
clc;

% Parameters for the neutron star
M = 1e30; % Mass of the neutron star (arbitrary units)
r0 = [0, 0]; % Position of the neutron star (arbitrary units)
sigma = 0.2; % Standard deviation for Gaussian (arbitrary units)

% Create a grid for the heat map
x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);

% Initialize the figure
figure;
h = surf(X, Y, zeros(size(X))); % Initialize the surface plot
shading interp; % Smooth the surface
colorbar; % Add color bar to indicate density values
xlabel('X Position (arbitrary units)');
ylabel('Y Position (arbitrary units)');
zlabel('Density (arbitrary units)');
title('3D Heat Map of a Single Neutron Star');
view(3); % Set view to 3D

% Animation loop
for t = 1:100
    % Update the Gaussian distribution for the neutron star
    Z = M * exp(-((X - r0(1)).^2 + (Y - r0(2)).^2) / (2 * sigma^2));
    
    % Add a time-dependent perturbation (optional)
    Z = Z .* (1 + 0.1 * sin(2 * pi * t / 50)); % Example: oscillating density
    
    % Update the surface plot
    set(h, 'ZData', Z, 'CData', Z); % Update Z and color data
    
    % Pause to control the animation speed
    pause(0.05);
end
