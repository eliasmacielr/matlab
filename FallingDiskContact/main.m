equations;

load sym_equations
%% Simulation parameters and Integration

R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;

alpha = 0; % dissipation parameter

t0 = 0;
tf = 50;
T = tf - t0;
h = 0.1;
N = int32(T/h) + 1;
tol = 1e-6;

span = [.8 1.2];

F = vpa(subs([eq_theta; eq_phi; eq_psi; Omega_d_1; Omega_d_2]));

% Initial conditions
% theta0 = 20*(pi/180);
% phidot0 = -0.15*(2*pi);
% psidot0 = ((I_T-I_A-m*R^2)*sin(theta0)*phidot0^2-m*g*R)/((I_A+m*R^2)*tan(theta0)*phidot0);
q0 = [0; 0; 0; 0; 0];
qdot0 = [0; 0; 0; 0; 0];

q = zeros(5, N);
q(:,1) = vpa(q0);
% Get q(:,2) using the disk's differential equations
[y0,yp0] = decic(@diskODEs,0,[q0;qdot0],[0 0 0 0 0 0 0 0 0 0], ...
    [qdot0;zeros(5,1)],[0 0 0 0 0 0 0 0 0 0]);
[~,y] = ode15i(@diskODEs,[t0,t0+h/2,t0+h],y0,yp0,odeset('RelTol',tol));
q(:,2) = vpa(transpose(y(end,1:5)));

for j = 2:N-1
    tj = t0 + (j-1)*h;

    X_i     = q(1,j-1);
    Y_i     = q(2,j-1);
    theta_i = q(3,j-1);
    phi_i   = q(4,j-1);
    psi_i   = q(5,j-1);

    X_j     = q(1,j);
    Y_j     = q(2,j);
    theta_j = q(3,j);
    phi_j   = q(4,j);
    psi_j   = q(5,j);

    [q(:,j+1), i] = newton_n_dim(q(:,j), q_k, subs(F), tol, 10);
%     fprintf("%d\n", i);
end

%% Animation and state space portraits
% animate_rolling_disk(q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),R,h);

% Get coordinates
X = q(1,:);
Y = q(2,:);
theta = q(3,:);
phi = q(4,:);
psi = q(5,:);

% Compute velocities from coordinates
Xdot = [y0(6),diff(X)/h];
Ydot = [y0(7),diff(Y)/h];
thetadot = [y0(8),diff(theta)/h];
phidot = [y0(9),diff(phi)/h];
psidot = [y0(10),diff(psi)/h];

t = t0:h:tf;

figure
subplot(2,1,1)
plot(t,X,t,Y,t,theta,t,phi,t,psi)
legend({'$X$','$Y$','$\theta$','$\phi$','$\psi$'},'Interpreter','latex')
xlabel('Tiempo (s)')
title(strcat('Configuraci{\''o}n del sistema, $q(0) = ',...
    latex(sym(y0(1:5)')), '$'), 'Interpreter', 'latex')
subplot(2,1,2)
plot(t,Xdot,t,Ydot,t,thetadot,t,phidot,t,psidot)
legend({'$\dot{X}$','$\dot{Y}$','$\dot{\theta}$','$\dot{\phi}$', ...
    '$\dot{\psi}$'}, 'Interpreter','latex')
xlabel('Tiempo (s)')
title(strcat('Velocidad del sistema, $\dot{q}(0) = ', ...
    latex(sym(yp0(1:5)')), '$'), ...
    'Interpreter', 'latex')

E = 1/2*m*(Xdot.^2 + Ydot.^2 + R^2*sin(theta).*thetadot.^2) + ...
    1/2*(I_A*(psidot - phidot.*sin(theta)).^2 + ...
        I_T*(thetadot.^2 + phidot.^2.*(cos(theta).^2))) + ...
    m*g*R*cos(theta);

figure
set(gcf, 'color', 'w')
plot(t, E, '-b', 'linewidth', 2)
xlabel('Tiempo (s)')
ylabel('Energía mecánica total (J)')
ylim([min(min(E)*span), max(max(E)*span)])

% savefig(strcat('resultados-',num2str(tol),'-',num2str(T),'s.fig'));

%% Save last simulation results
save(strcat('res-contact2-h',num2str(h),'-alpha',num2str(alpha),'.mat'),'t0','tf','h','q','y0');
