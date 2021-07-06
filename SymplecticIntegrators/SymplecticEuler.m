% Planar Pendulum Symplectic Integrator
clear
clc
clf

% Symbols, variables and values
syms theta theta_t omega_0

g = 9.81;   % Gravity (ms-2)
l = 4;      % Pendulum length (m)
T = 10;     % Total time (s)
h = 0.1;    % Timestep (s)
N = T/h;

gValue = 9.81;
rValue = l;
omega_0Value = sqrt(gValue/rValue);

q = zeros(N,1);
v = zeros(N,1);
% Initial conditions
q(1) = 1;    % pendulum's angle with the vertical
v(1) = 0;
for k = 1 : N-1
    v(k+1) = v(k) - h * g/l * sin(q(k));
    q(k+1) = q(k) + h * v(k+1);
end

E(theta, theta_t, omega_0) = (1/2)*(theta_t^2+(2*omega_0*sin(theta/2))^2);
Eplot(theta, theta_t) = subs(E,omega_0,omega_0Value);

axis equal
plot(q,v,'bo')
xlim([-2.5 2.5])
ylim([-3 3])
hold on
plot(q(1),v(1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b',...
    'MarkerSize',8)
hold on
fc = fcontour(Eplot(theta, theta_t), [-2.5 2.5 -3 3],'LineColor','black');
ax = gca;
ax.FontSize = 16;
xlabel('$q$','FontSize',24,'Interpreter','latex')
ylabel('$\dot{q}$','FontSize',24,'Interpreter','latex')
