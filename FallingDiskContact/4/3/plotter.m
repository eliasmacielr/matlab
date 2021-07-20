clearvars

figure

y0 = zeros(10,1);
load res-contact2-h0.1-alpha0.1.mat

t = t0:h:tf;
X = q(1,:);
Y = q(2,:);
theta = q(3,:);
phi = q(4,:);
psi = q(5,:);
subplot(2,1,1)
plot(t,X,t,Y,t,theta,t,phi,t,psi)
xlim([t0 tf])
legend({'$X$','$Y$','$\theta$','$\phi$','$\psi$'},'Interpreter','latex')
xlabel('Tiempo (s)')

Xdot = [y0(6),diff(X)/h];
Ydot = [y0(7),diff(Y)/h];
thetadot = [y0(8),diff(theta)/h];
phidot = [y0(9),diff(phi)/h];
psidot = [y0(10),diff(psi)/h];
subplot(2,1,2)
plot(t,Xdot,t,Ydot,t,thetadot,t,phidot,t,psidot)
xlim([t0 tf])
legend({'$\dot{X}$','$\dot{Y}$','$\dot{\theta}$','$\dot{\phi}$', ...
    '$\dot{\psi}$'}, 'Interpreter','latex')
xlabel('Tiempo (s)')

R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;

Econtacto = 1/2*m*(Xdot.^2 + Ydot.^2 + R^2*sin(theta).*thetadot.^2) + ...
    1/2*(I_A*(psidot - phidot.*sin(theta)).^2 + ...
        I_T*(thetadot.^2 + phidot.^2.*(cos(theta).^2))) + ...
    m*g*R*cos(theta);

savefig(strcat('contacto.fig'));

% 

% clearvars -except Econtacto
% 
% figure
% 
% y0 = zeros(10,1);
% load res-ode15i-h0.1-alpha0.mat
% 
% t = t0:h:tf;
% X = transpose(q(:,1));
% Y = transpose(q(:,2));
% theta = transpose(q(:,3));
% phi = transpose(q(:,4));
% psi = transpose(q(:,5));
% t = (t0+h)*(0:max(size(X))-1);
% subplot(2,1,1)
% plot(t,X,t,Y,t,theta,t,phi,t,psi)
% xlim([t0 tf])
% legend({'$X$','$Y$','$\theta$','$\phi$','$\psi$'},'Interpreter','latex')
% xlabel('Tiempo (s)')
% 
% Xdot = [y0(6),diff(X)/h];
% Ydot = [y0(7),diff(Y)/h];
% thetadot = [y0(8),diff(theta)/h];
% phidot = [y0(9),diff(phi)/h];
% psidot = [y0(10),diff(psi)/h];
% subplot(2,1,2)
% plot(t,Xdot,t,Ydot,t,thetadot,t,phidot,t,psidot)
% xlim([t0 tf])
% legend({'$\dot{X}$','$\dot{Y}$','$\dot{\theta}$','$\dot{\phi}$', ...
%     '$\dot{\psi}$'}, 'Interpreter','latex')
% xlabel('Tiempo (s)')
% 
% R = 0.5;
% m = 5;
% I_A = 1/2*m*R^2;
% I_T = 1/4*m*R^2;
% g = 9.81;
% 
% Eode15i = 1/2*m*(Xdot.^2 + Ydot.^2 + R^2*sin(theta).*thetadot.^2) + ...
%     1/2*(I_A*(psidot - phidot.*sin(theta)).^2 + ...
%         I_T*(thetadot.^2 + phidot.^2.*(cos(theta).^2))) + ...
%     m*g*R*cos(theta);
% 
% savefig(strcat('ode15i.fig'));

% 

clearvars -except Econtacto Eode15i

figure

y0 = zeros(10,1);
load res-ode15i-h0.01-alpha0.1.mat

t = t0:h:tf;
X = transpose(q(:,1));
Y = transpose(q(:,2));
theta = transpose(q(:,3));
phi = transpose(q(:,4));
psi = transpose(q(:,5));
t = (t0+h)*(0:max(size(X))-1);
subplot(2,1,1)
plot(t,X,t,Y,t,theta,t,phi,t,psi)
xlim([t0 tf])
legend({'$X$','$Y$','$\theta$','$\phi$','$\psi$'},'Interpreter','latex')
xlabel('Tiempo (s)')

Xdot = [y0(6),diff(X)/h];
Ydot = [y0(7),diff(Y)/h];
thetadot = [y0(8),diff(theta)/h];
phidot = [y0(9),diff(phi)/h];
psidot = [y0(10),diff(psi)/h];
subplot(2,1,2)
plot(t,Xdot,t,Ydot,t,thetadot,t,phidot,t,psidot)
xlim([t0 tf])
legend({'$\dot{X}$','$\dot{Y}$','$\dot{\theta}$','$\dot{\phi}$', ...
    '$\dot{\psi}$'}, 'Interpreter','latex')
xlabel('Tiempo (s)')

R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;

Eref = 1/2*m*(Xdot.^2 + Ydot.^2 + R^2*sin(theta).*thetadot.^2) + ...
    1/2*(I_A*(psidot - phidot.*sin(theta)).^2 + ...
        I_T*(thetadot.^2 + phidot.^2.*(cos(theta).^2))) + ...
    m*g*R*cos(theta);

savefig(strcat('bench.fig'));

%

clearvars -except Econtacto Eref t0 tf

h = 0.1;
t = t0:h:tf;
% tode15i = (t0+h)*(0:max(size(Eode15i))-1);
th10 = (t0+(h/10))*(0:max(size(Eref))-1);

figure
set(gcf, 'color', 'w')
plot(t, Econtacto, ':bs')
% hold on
% plot(tode15i, Eode15i, '-r*')
hold on
plot(th10, Eref, '-k')
xlim([t0 tf])
legend({'Contacto (orden 2)','Referencia'}, 'Interpreter','latex')
xlabel('Tiempo (s)')
ylabel('Energía mecánica total (J)')

savefig(strcat('energia.fig'));
