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
q(1) = -2;    % pendulum's angle with the vertical
v(1) = 2;
for k = 1 : N-1
    v(k+1) = v(k) - h * g/l * sin(q(k));
    q(k+1) = q(k) + h * v(k+1);
end

E(theta, theta_t, omega_0) = (1/2)*(theta_t^2+(2*omega_0*sin(theta/2))^2);
Eplot(theta, theta_t) = subs(E,omega_0,omega_0Value);

axis equal
plot(q,v,'bo')
xlim([-4 10])
ylim([-5 5])
hold on
plot(q(1),v(1),'ko','MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
fc = fcontour(Eplot(theta, theta_t), [-4 10 -5 5], 'LineWidth', 1);
xlabel('$q$','fontsize',24,'Interpreter','latex')
ylabel('$\dot{q}$','fontsize',24,'Interpreter','latex')
