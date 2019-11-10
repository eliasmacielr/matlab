clear all %#ok<CLALL>
close all
clc

T = 5; % Time of the simulation (s)
dt = 0.01; % Time step (s)
N = T/dt;

x = zeros(N, 1);
y = zeros(N, 1);
z = zeros(N, 1);
phi = zeros(N, 1);
theta = zeros(N, 1);

% Disk
R = 0.5; % Radius (m)
sides = 48;
angle = linspace(0, 2*pi, sides)';
circ1 = R*cos(angle);
circ2 = R*sin(angle) + R;
xDisk = zeros(sides,N);
yDisk = zeros(sides,N);
zDisk = zeros(sides,N);

% Initial conditions
x0 = 0;
y0 = 0;
phi0 = pi/2;
theta0 = 0;
omega = pi/4;
Omega = 4*sqrt(2);
phi(1) = phi0;
for t = 1:N
    x(t) = (Omega/omega)*R*sin(omega*t*dt + phi0) + x0;
    y(t) = -(Omega/omega)*R*cos(omega*t*dt + phi0) + y0;
    phi(t) = omega*t*dt + phi0;
    theta(t) = Omega*t*dt + theta0;
    % TODO: apply rotation matrices
    xDisk(:,t) = x(t) + zeros(sides, 1);
    yDisk(:,t) = y(t) + circ1;
    zDisk(:,t) = z(t) + circ2;
end

% Setup video
figure
set(gcf, 'color', 'w')
plot3(x(1), y(1), 0);
xlabel('\itx (m)')
ylabel('\ity (m)')
zlabel('\itz (m)            ', 'rotation', 0)
axis equal
xlim([min(x)-R, max(x)+R])
ylim([min(y)-R, max(y)+R])
zlim([0, 1.2*(2*R)])
grid on

path = line('xdata', x(1:1), 'ydata', y(1:1), 'zdata', z(1:1), 'color', 'b', 'linewidth', 2);
disk = patch('xdata', xDisk(:,1), 'ydata', yDisk(:,1), 'zdata', zDisk(:,1), 'facecolor', [1, 0.8, 0.6], 'linewidth', 2);

animation = VideoWriter('vertical-rolling-disk.avi');
animation.FrameRate = 1/dt;
open(animation);

for t = 1:N
    set(path, 'xdata', x(1:t), 'ydata', y(1:t), 'zdata', z(1:t));
    set(disk, 'xdata', xDisk(:,t), 'ydata', yDisk(:,t), 'zdata', zDisk(:,t));
    drawnow
    writeVideo(animation, getframe(gcf));
end 

close(animation);
