
%% Define constants
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
M = 1.989e30; % Mass of the Sun (kg)
c = 299792458; % Speed of light (m/s)
A = 147098290000; % Semi-major axis of Earth's orbit (m)
AU = 149.6e9; % Astronomical unit [m]

a_earth = 1.00000011*AU; % Semi-major axis [m]
e_earth = 0.0167; % Eccentricity
b_earth = a_earth*sqrt(1-e_earth^2);
c_earth = b_earth; % Assume circular orbit for simplicity

% Earth's orbit parameters
a = 1.496e11; % Semi-major axis [m]
e = 0.0167; % Eccentricity  

% Define the angles for April 2023 and May 2023
theta_apr = 193.46712;
theta_may = 221.47994;
theta_jun = 252.41449;

%%

% Calculate the distance between the Earth and Sun at each angle
r_apr = 1/((G*M)/c^2)-A*cosd(theta_apr);
r_may = 1/((G*M)/c^2)-A*cosd(theta_may);
r_jun = 1 /((G*M)/c^2)-A*cosd(theta_jun);

% x, y, and z coordinates of the Earth at each angle
x_apr = r_apr*cosd(theta_apr);
y_apr = r_apr*sind(theta_apr);
z_apr = 0;

x_may = r_may*cosd(theta_may);
y_may = r_may*cosd(theta_may);
z_may = 0;

x_jun = r_jun*cosd(theta_jun);
y_jun = r_jun*sind(theta_jun);
z_jun = 0;

% Position of each body
figure;
hold on;
plot3(x_apr, y_apr, z_apr, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot3(x_may, y_may, z_may, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(x_jun, y_jun, z_jun, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot3(0, 0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k') % plot Sun at origin

% angle for earth orbit
theta_ellipse = linspace(0, 2*pi, 500);

% earth orbit xyz
x_ellipse = a_earth*cos(theta_ellipse);
y_ellipse = b_earth*sin(theta_ellipse);
z_ellipse = zeros(size(theta_ellipse));

% earth orbit
plot3(x_ellipse, y_ellipse, z_ellipse, 'g', 'LineWidth', 1.5)

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Earths Position in April and May 2023');
grid on;
axis equal;
view(3);

legend('4/5/2023','5/5/2023', 'Sun (Sol)')

set(gca,'color',[0 0 0])

