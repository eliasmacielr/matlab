% Planar Pendulum Symplectic Integrator
clear
clc
clf

g = 9.81;   % Gravity (ms-2)
l = 4;      % Pendulum length (m)
T = 20;     % Total time (s)
h = 0.1;    % Timestep (s)
N = T/h;

q = zeros(N,1);
v = zeros(N,1);
% Initial conditions
q(1) = pi/4;    % pendulum's angle with the vertical
for k = 1 : N-1
    v(k+1) = v(k) - h * g/l * sin(q(k));
    q(k+1) = q(k) + h * v(k+1);
end

% Angle (q) to Cartesian coordinates
x = l*sin(q);
y = -l*cos(q);

for i = 1 : N
    subplot(2,1,1)
    plot(x(i),y(i),'ko',[0 x(i)],[0 y(i)],'r-')
    title(['Planar pendulum simulation (q = ' num2str(q(i)) ')'],'fontsize',12)
    xlabel('x[m]','fontsize',12)
    ylabel('y[m]','fontsize',12)
    axis equal
    axis([-l l -l 0])

    subplot(2,1,2)
    plot(q(i),v(i),'bo')
    hold on
    title('Planar pendulum phase portrait','fontsize',12)
    xlabel('q','fontsize',12)
    ylabel('q.','fontsize',12)
    axis equal
    axis([-q(1) q(1) min(v) max(v)])

    pause(h)  % Shows results at each time interval
end
