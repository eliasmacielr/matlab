equations;

load sym_equations
%% Simulation parameters and Integration

R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;
alpha = 0;

t0 = 0;
tf = 5;
T = tf - t0;
h = 0.1;
N = int32(T/h);
tol = 1e-6;

F = vpa(subs([eq_theta; eq_phi; eq_psi; Omega_d_1; Omega_d_2]));

% Initial conditions
q0 = [0; 0; 0; 0; 0];
qdot0 = [0; 0; 0; 0; pi];

q = zeros(5, N);
q(:,1) = vpa(q0);
% Get q(:,2) using the disk's differential equations (it does not work)
[y0,yp0] = decic(@diskODEs,0,[q0;qdot0],[0 0 1 0 0 0 0 0 0 1],...
    [qdot0;zeros(5,1)],[0 0 0 0 0 0 0 0 0 0]);
[~,y] = ode15i(@diskODEs,[0,h/2,h],y0,yp0);
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
    fprintf("%d\n", i);
end

%% Animation and state space portraits
animate_rolling_disk(q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),R,h);

figure

subplot(2,1,1);
t = linspace(0,T,N);
plot(t,q(1,:),t,q(2,:),t,q(3,:),t,q(4,:),t,q(5,:))
legend({'$X$','$Y$','$\theta$','$\phi$','$\psi$'},'Interpreter','latex')
xlabel('Tiempo (s)')
title(strcat('Posici{\''o}n del sistema, $q(0) = ', latex(sym(q0')), '$'), ...
    'Interpreter', 'latex')

subplot(2,1,2);
% TODO: ver como sacar las velocidades
t = linspace(0,T,N-1);
plot(t,diff(q(1,:))/h,t,diff(q(2,:))/h,t,diff(q(3,:))/h, ...
    t,diff(q(4,:))/h,t,diff(q(5,:))/h);
legend({'$\dot{X}$','$\dot{Y}$','$\dot{\theta}$','$\dot{\phi}$', ...
    '$\dot{\psi}$'}, 'Interpreter','latex')
xlabel('Tiempo (s)')
title(strcat('Velocidad del sistema, $\dot{q}(0) = ', ...
    latex(sym(qdot0')), '$'), ...
    'Interpreter', 'latex')

% savefig(strcat('resultados-',num2str(tol),'-',num2str(T),'s.fig'));

%% Save last simulation results
% save(strcat('vars-',num2str(tol),'-',num2str(T),'s.mat'),'q','R','h');
