% TODO: obtain the equations and then set the values for the simulation
% TODO: write documentation and describe the problem and where it came from

clearvars
close all
clc

%% Symbols
% q_{j-1}
syms X_i Y_i theta_i phi_i psi_i z_i
% q_j
syms X_j Y_j theta_j phi_j psi_j z_j
% q_{j+1}
syms X_k Y_k theta_k phi_k psi_k z_k
% constants
syms R m I_A I_T g alpha
% "Control" forces
u_psi = 0;
% Lagrange multipliers
syms lambda_1 lambda_2
% time and step size
syms tj h

assume((R > 0) & (m > 0) & (I_A > 0) & (I_T > 0) & (g > 0) & ...
    (alpha > 0) & (h > 0))

%% Discrete Lagrangians and nonholonomic constraints
% Write L(q_j,q_{j+1},z_j,z_{j+1})
% q (X, Y and psi are not used here)
theta = theta_j;
phi = phi_j;
% q_dot
X_dot = (X_k-X_j)/h;
Y_dot = (Y_k-Y_j)/h;
theta_dot = (theta_k-theta_j)/h;
phi_dot = (phi_k-phi_j)/h;
psi_dot = (psi_k-psi_j)/h;

f_j = [0 0 0*tj 0*tj 0*tj];
f_k = [0 0 0*(tj+h) 0*(tj+h) 0*(tj+h)];
q_j = [X_j Y_j theta_j phi_j psi_j]';
q_k = [X_k Y_k theta_k phi_k psi_k]';

L_j = 1/2*m*(X_dot^2 + Y_dot^2 + R^2*(sin(theta))^2*theta_dot^2) + ...
    1/2*(I_A*(psi_dot - phi_dot*sin(theta))^2 + ...
         I_T*(theta_dot^2 + phi_dot^2*(cos(theta))^2)) - ...
    m*g*R*(cos(theta_j)+cos(theta_k))/2 - alpha*z_j + ...
    (f_j*q_j+f_k*q_k)/2;

% Write constraints (uses q_j and q_{j+1})
Omega_d_1 = X_dot + R*cos(theta)*sin(phi)*theta_dot + ...
    R*sin(theta)*cos(phi)*phi_dot - R*cos(phi)*psi_dot; % (1)
Omega_d_2 = Y_dot - R*cos(theta)*cos(phi)*theta_dot + ...
    R*sin(theta)*sin(phi)*phi_dot - R*sin(phi)*psi_dot; % (2)

% Write L(q_{j-1},q_j,z_{j-1},z_j)
% q
X = X_i;
Y = Y_i;
theta = theta_i;
phi = phi_i;
psi = psi_i;
% q_dot
X_dot = (X_j-X_i)/h;
Y_dot = (Y_j-Y_i)/h;
theta_dot = (theta_j-theta_i)/h;
phi_dot = (phi_j-phi_i)/h;
psi_dot = (psi_j-psi_i)/h;

f_i = [0 0 0*(tj-h) 0*(tj-h) 0*(tj-h)];
f_j = [0 0 0*tj 0*tj 0*tj];
q_i = [X_i Y_i theta_i phi_i psi_i]';
q_j = [X_j Y_j theta_j phi_j psi_j]';

L_i = 1/2*m*(X_dot^2 + Y_dot^2 + R^2*(sin(theta))^2*theta_dot^2) + ...
    1/2*(I_A*(psi_dot - phi_dot*sin(theta))^2 + ...
         I_T*(theta_dot^2 + phi_dot^2*(cos(theta))^2)) - ...
    m*g*R*(cos(theta_i)+cos(theta_j))/2 - alpha*z_i + ...
    (f_i*q_i+f_j*q_j)/2;

% Precompute the term that is multiplied to D_2 L(q_{j-1},q_j,...)
partial_z = (1 + h*diff(L_j,z_j))/(1 - h*diff(L_i,z_j));

% Write system of equations in the form F(q) = 0, which with the
% constraint equations (1) and (2) make up the whole system to be
% solved.

lambda_1 = diff(L_j,X_j) + diff(L_i,X_j)*partial_z;
lambda_2 = diff(L_j,Y_j) + diff(L_i,Y_j)*partial_z;

eq_theta = diff(L_j,theta_j) + diff(L_i,theta_j)*partial_z - ...
    lambda_1*R*cos(theta_j)*sin(phi_j) + ...
    lambda_2*R*cos(theta_j)*cos(phi_j);
eq_phi = diff(L_j,phi_j) + diff(L_i,phi_j)*partial_z - ...
    lambda_1*R*sin(theta_j)*cos(phi_j) - ...
    lambda_2*R*sin(theta_j)*sin(phi_j);
eq_psi = diff(L_j,psi_j) + diff(L_i,psi_j)*partial_z + ...
    lambda_1*R*cos(phi_j) + lambda_2*R*sin(phi_j);

%% Initial conditions and Integration
R = 1;
m = 1;
I_A = 1;
I_T = 1;
g = 9.8;
alpha = 0.1;

t0 = 0;
tf = 3;
T = tf - t0;
h = 0.2;
N = int32(T/h);
tol = 1e-4;

eq_theta = subs(eq_theta);
eq_phi = subs(eq_phi);
eq_psi = subs(eq_psi);
Omega_d_1 = subs(Omega_d_1);
Omega_d_2 = subs(Omega_d_2);

q = zeros(5, N);

% Get integrator equations with the values for q

% q(:,j) = [X(j);Y(j);theta(j);phi(j);psi(j)]
q(:,1) = [0; 0; 0; 0; 0];
% TODO: obtener q(:,2) a partir de una aproximacion
q(:,2) = [R*pi/12; 0; 0; 0; pi/12];

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

    [q(:,j+1), i] = newton_n_dim(q(:,j), [X_k; Y_k; theta_k; phi_k; psi_k], ...
        subs([eq_theta; eq_phi; eq_psi; Omega_d_1; Omega_d_2]), tol, 10);
    fprintf("%d\n", i);
end

%% Animation and state space portraits
animate_rolling_disk(q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),R,h);

figure

subplot(2,1,1);
t = linspace(0,T,N);
plot(t,q(1,:),t,q(2,:),t,q(3,:),t,q(4,:),t,q(5,:))
legend({'X','Y','\theta','\phi','\psi'},'Interpreter','tex')
xlabel('Tiempo (s)')
title('Posición del sistema, q(0) = (0,0,\pi/6,0,0)','Interpreter','tex')

subplot(2,1,2);
% Todo, ver como sacar las velocidades
t = linspace(0,T,N-1);
plot(t,diff(q(1,:))/h,t,diff(q(2,:))/h,t,diff(q(3,:))/h, ...
    t,diff(q(4,:))/h,t,diff(q(5,:))/h);
legend({'X^{''}','Y^{''}','\theta^{''}','\phi^{''}','\psi^{''}'}, ...
    'Interpreter','tex')
xlabel('Tiempo (s)')
title('Velocidad del sistema, q''(0) = (0,0,\pi/6,0,0)','Interpreter', ...
    'tex')

savefig(strcat('resultados-',num2str(tol),'-',num2str(T),'s'));

%% Save last simulation results
save(strcat('vars-',num2str(tol),'-',num2str(T),'s.mat'),'q','R','h');
