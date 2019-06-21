clear
clc
clf

g = 9.81;   % Gravity (ms-2)
l = 2;      % Pendulum length (m)
T = 5;      % Total time (s)
h = 0.01;   % Timestep (s)
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
    subplot(1,1,1)
    plotarrayx = [0 x(i)];
    plotarrayy = [0 y(i)];
    plot(x(i),y(i),'ko',plotarrayx,plotarrayy,'r-')
    title(['Planar pendulum simulation (q = ' num2str(q(i)) ')'],'fontsize',12)
    xlabel('x[m]','fontsize',12)
    ylabel('y[m]','fontsize',12)
    axis([-l l -l 0])
    
    pause(h)  % Shows results at each time interval
end
