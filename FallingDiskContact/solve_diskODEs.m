R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;
alpha = 0; % dissipation parameter

t0 = 0;
tf = 10;
T = tf - t0;
h = 0.1;
% N = int32(T/h);
tol = 1e-6;

q0 = [0; 0; 0; 0; 0];
qdot0 = [0; 0; 0; 0; pi];
[y0,yp0] = decic(@diskODEs,0,[q0;qdot0],[0 0 1 0 0 0 0 0 0 1], ...
    [qdot0;zeros(5,1)],[0 0 0 0 0 0 0 0 0 0]);
[t,y] = ode15i(@diskODEs,t0:h:tf,y0,yp0,odeset('RelTol',tol));

X = y(:,1);
Y = y(:,2);
theta = y(:,3);
phi = y(:,4);
psi = y(:,5);

figure

subplot(2,1,1);
plot(t,X,t,Y,t,theta,t,phi,t,psi)
legend({'$X$','$Y$','$\theta$','$\phi$','$\psi$'},'Interpreter','latex')
xlabel('Tiempo (s)')
title(strcat('Configuraci{\''o}n del sistema, $q(0) = ',...
    latex(sym(y0(1:5)')), '$'), 'Interpreter', 'latex')

subplot(2,1,2);
plot(t,[y0(6);diff(X)/h],t,[y0(7);diff(Y)/h],t,[y0(8);diff(theta)/h], ...
    t,[y0(9);diff(phi)/h],t,[y0(10);diff(psi)/h]);
legend({'$\dot{X}$','$\dot{Y}$','$\dot{\theta}$','$\dot{\phi}$', ...
    '$\dot{\psi}$'}, 'Interpreter','latex')
xlabel('Tiempo (s)')
title(strcat('Velocidad del sistema, $\dot{q}(0) = ', ...
    latex(sym(yp0(1:5)')), '$'), ...
    'Interpreter', 'latex')
