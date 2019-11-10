clear all %#ok<CLALL>
close all
clc

T = 10; % Time of the simulation (s)
dt = 0.01; % Time step (s)
N = T/dt;

% TODO: comment all these
R = 0.5; % Radius of the disk (m)
x = zeros(N,1);
y = zeros(N,1);
z = zeros(N,1);
phi = zeros(N,1);
theta = zeros(N,1);

% Initial conditions
x0 = 0;
y0 = 0;
phi0 = pi/2;
theta0 = 0;
omega = pi/9;
Omega = 4*sqrt(2);
phi(1) = phi0;
for t = 1:N
    x(t) = (Omega/omega)*R*sin(omega*t*dt + phi0) + x0;
    y(t) = -(Omega/omega)*R*cos(omega*t*dt + phi0) + y0;
    phi(t) = omega*t*dt + phi0;
    theta(t) = Omega*t*dt + theta0;
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
% disk = patch('xdata', xcirc(:,1), 'ydata', ycirc(:,1), 'zdata', zcirc(:,1), 'facecolor', [1, 0.8, 0.6], 'linewidth', 2);
% pointA = line('xdata', xA(1), 'ydata', yA(1), 'zdata', zA(1), 'marker', 'o', 'color', 'r', 'markerfacecolor', 'r', 'linewidth', 3);

animation = VideoWriter('vertical-rolling-disk.avi');
animation.FrameRate = 1/dt;
open(animation);

for t = 1:N
    set(path, 'xdata', x(1:t), 'ydata', y(1:t), 'zdata', z(1:t));
%     set(disk, 'xdata', xcirc(:,t), 'ydata', ycirc(:,t), 'zdata', zcirc(:,t));
%     set(pointA, 'xdata', xA(t), 'ydata', yA(t), 'zdata', zA(t));
    drawnow
    writeVideo(animation, getframe(gcf));
end 

close(animation);
