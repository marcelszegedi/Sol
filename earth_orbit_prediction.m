%% Earth-Sun model
clear; clc; close all;

%% Parameters
n = 200; % number of time steps
t = linspace(0, 2*pi, n); % time vector
a = 1; % semi-major axis of Earth's orbit (AU)
b = 0.9833; % semi-minor axis of Earth's orbit (AU)
e = sqrt(1 - (b/a)^2); % eccentricity of Earth's orbit
T = 1; % period of Earth's orbit (yr)
w = 2*pi/T; % angular velocity of Earth's orbit (rad/yr)

%% Earth's orbit
% Define a time vector for Earth's orbit
t_orbit = linspace(0, 2*pi, n);

% Calculate the x, y, and z coordinates of Earth's orbit
x_orbit = a*(cos(t_orbit) + e);
y_orbit = b*sin(t_orbit);
z_orbit = zeros(size(x_orbit));

% Find the indices in the time vector that correspond to April 1st 2023 and April 1st 2024
t_apr_2023_idx = round(interp1(linspace(0,2*pi,n),1:n,(2/365.25)*2*pi));
t_apr_2024_idx = round(interp1(linspace(0,2*pi,n),1:n,(3/365.25)*2*pi));

%% Visualization
figure;
set(gcf, 'color', 'w');
hold on;
grid on;
axis equal;
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
zlim([-1.5, 1.5]);
xlabel('X (AU)');
ylabel('Y (AU)');
zlabel('Z (AU)');
view(30, 30);
light('Position', [1 1 0]);

% Sun
r_sun = 0.1; % radius of Sun (AU)
sun_color = [1 1 0]; % yellow
[xs, ys, zs] = sphere(100);
xs = xs * r_sun;
ys = ys * r_sun;
zs = zs * r_sun;
surf(xs, ys, zs, 'FaceColor', sun_color, 'EdgeColor', 'none');

% Earth
r_earth = 0.02; % radius of Earth (AU)
earth_color = [0 0 1]; % blue
for i = 1:length(x_orbit)
    earth_pos = [x_orbit(i); y_orbit(i); z_orbit(i)]; % position of Earth
    xe = earth_pos(1) + r_earth*cos(linspace(0, 2*pi, 100)); % x-coordinates of Earth points
    ye = earth_pos(2) + r_earth*sin(linspace(0, 2*pi, 100)); % y-coordinates of Earth points
    ze = earth_pos(3)*ones(size(xe)); % z-coordinates of Earth points
    plot3(xe, ye, ze, '.', 'Color', earth_color);
    pause(0.05);
    drawnow;
end

% Position of Earth on April 1st, 2023
n_apr_2023 = round(n * (datetime('01-Apr-2023') - datetime('01-Jan-2023')) / T);
x_apr_2023 = a*(cos(t(n_apr_2023)) + e);
y_apr_2023 = b*sin(t(n_apr_2023));
z_apr_2023 = 0;
plot3(x_apr_2023, y_apr_2023, z_apr_2023, 'ro', 'MarkerSize', 50, 'LineWidth', 2);

% Position of Earth on April 1st, 2024
n_apr_2024 = round(n * (datetime('01-Apr-2024') - datetime('01-Jan-2023')) / T);
x_apr_2024 = a*(cos(t(n_apr_2024)) + e);
y_apr_2024 = b*sin(t(n_apr_2024));
z_apr_2024 = 0;
plot3(x_apr_2024, y_apr_2024, z_apr_2024, 'ro', 'MarkerSize', 50, 'LineWidth', 2);


% Red markers for April 1st 2023 and April 1st 2024
marker_size = 50;
plot3(x_orbit(t_apr_2023_idx), y_orbit(t_apr_2023_idx), z_orbit(t_apr_2023_idx), '.', 'Color', 'r', 'MarkerSize', marker_size);
earth_pos = [x(i); y(i); z(i)];

