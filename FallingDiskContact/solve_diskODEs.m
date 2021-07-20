R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;

t0 = 0;
tf = 25;
h = 0.01;
tol = 1e-6;

span = [.8 1.2];

theta0 = 20*(pi/180);
phidot0 = -0.15*(2*pi);
psidot0 = ((I_T-I_A-m*R^2)*sin(theta0)*phidot0^2-m*g*R)/((I_A+m*R^2)*tan(theta0)*phidot0);
q0 = [0; 0; theta0; 0; 0];
qdot0 = [pi/2; 0; 0; phidot0; psidot0];
[y0,yp0] = decic(@diskODEs,0,[q0;qdot0],[0 0 0 0 0 0 0 0 0 0], ...
    [qdot0;zeros(5,1)],[0 0 0 0 0 0 0 0 0 0]);
[t,y] = ode15i(@diskODEs,t0:h:tf,y0,yp0,odeset('RelTol',tol));

% Get coordinates
X = y(:,1);
Y = y(:,2);
theta = y(:,3);
phi = y(:,4);
psi = y(:,5);

% Compute velocities from coordinates
Xdot = [y0(6);diff(X)/h];
Ydot = [y0(7);diff(Y)/h];
thetadot = [y0(8);diff(theta)/h];
phidot = [y0(9);diff(phi)/h];
psidot = [y0(10);diff(psi)/h];

% figure
% subplot(2,1,1)
% plot(t,X,t,Y,t,theta,t,phi,t,psi)
% legend({'$X$','$Y$','$\theta$','$\phi$','$\psi$'},'Interpreter','latex')
% xlabel('Tiempo (s)')
% title(strcat('Configuraci{\''o}n del sistema, $q(0) = ',...
%     latex(sym(y0(1:5)')), '$'), 'Interpreter', 'latex')
% subplot(2,1,2)
% plot(t,Xdot,t,Ydot,t,thetadot,t,phidot,t,psidot)
% legend({'$\dot{X}$','$\dot{Y}$','$\dot{\theta}$','$\dot{\phi}$', ...
%     '$\dot{\psi}$'}, 'Interpreter','latex')
% xlabel('Tiempo (s)')
% title(strcat('Velocidad del sistema, $\dot{q}(0) = ', ...
%     latex(sym(yp0(1:5)')), '$'), ...
%     'Interpreter', 'latex')

E = 1/2*m*(Xdot.^2 + Ydot.^2 + R^2*sin(theta).*thetadot.^2) + ...
    1/2*(I_A*(psidot - phidot.*sin(theta)).^2 + ...
        I_T*(thetadot.^2 + phidot.^2.*(cos(theta).^2))) + ...
    m*g*R*cos(theta);

% figure
% set(gcf, 'color', 'w')
% plot(t, E, '-b', 'linewidth', 2)
% xlabel('Tiempo (s)')
% ylabel('Energía mecánica total (J)')
% ylim([min(min(E)*span), max(max(E)*span)])

q = [X, Y, theta, phi, psi];

save(strcat('res-ode15i-h',num2str(h),'-alpha',num2str(0.1),'.mat'),'t0','tf','h','q','y0');
