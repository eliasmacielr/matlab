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
q_j = [X_j; Y_j; theta_j; phi_j; psi_j];
q_k = [X_k; Y_k; theta_k; phi_k; psi_k];

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
q_i = [X_i; Y_i; theta_i; phi_i; psi_i];
q_j = [X_j; Y_j; theta_j; phi_j; psi_j];

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

save('equations','q_k',...
    'eq_theta','eq_phi','eq_psi','Omega_d_1','Omega_d_2');
