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
% Dissipation
syms alpha
% External forces
% F_{j-1}
syms F_X_i F_Y_i F_theta_i F_phi_i F_psi_i
% F_j
syms F_X_j F_Y_j F_theta_j F_phi_j F_psi_j
% F_{j+1}
syms F_X_k F_Y_k F_theta_k F_phi_k F_psi_k
% Constants
syms R m I_A I_T g
% Lagrange multipliers
syms lambda_1 lambda_2
% step size
syms h

assume((R > 0) & (m > 0) & (I_A > 0) & (I_T > 0) & (g > 0) & ...
    (alpha > 0) & (h > 0))

%% Discrete Lagrangians and nonholonomic constraints
% Write L(q_j,q_{j+1},z_j,z_{j+1})
% q (X, Y and psi are not used here)
theta = (theta_k+theta_j)/2;
phi = (phi_k+phi_j)/2;
% q_dot
X_dot = (X_k-X_j)/h;
Y_dot = (Y_k-Y_j)/h;
theta_dot = (theta_k-theta_j)/h;
phi_dot = (phi_k-phi_j)/h;
psi_dot = (psi_k-psi_j)/h;

F_j = [F_X_j F_Y_j F_theta_j F_phi_j F_psi_j];
F_k = [F_X_k F_Y_k F_theta_k F_phi_k F_psi_k];
q_j = [X_j; Y_j; theta_j; phi_j; psi_j];
q_k = [X_k; Y_k; theta_k; phi_k; psi_k];

L_j = 1/2*m*(X_dot^2 + Y_dot^2 + R^2*sin(theta)^2*theta_dot^2) + ...
    1/2*(I_A*(psi_dot - phi_dot*sin(theta))^2 + ...
         I_T*(theta_dot^2 + phi_dot^2*cos(theta)^2)) - ...
    m*g*R*(cos(theta_j)+cos(theta_k))/2 - alpha*(z_j+z_k)/2 + ...
    (F_j*q_j+F_k*q_k)/2;

% Write constraints (uses q_j and q_{j+1})
Omega_d_1 = X_dot + R*cos(theta)*sin(phi)*theta_dot + ...
    R*sin(theta)*cos(phi)*phi_dot - R*cos(phi)*psi_dot; % (1)
Omega_d_2 = Y_dot - R*cos(theta)*cos(phi)*theta_dot + ...
    R*sin(theta)*sin(phi)*phi_dot - R*sin(phi)*psi_dot; % (2)

% Write L(q_{j-1},q_j,z_{j-1},z_j)
% q
X = (X_j+X_i)/2;
Y = (Y_j+Y_i)/2;
theta = (theta_j+theta_i)/2;
phi = (phi_j+phi_i)/2;
psi = (psi_j+psi_i)/2;
% q_dot
X_dot = (X_j-X_i)/h;
Y_dot = (Y_j-Y_i)/h;
theta_dot = (theta_j-theta_i)/h;
phi_dot = (phi_j-phi_i)/h;
psi_dot = (psi_j-psi_i)/h;

F_i = [F_X_i F_Y_i F_theta_i F_phi_i F_psi_i];
q_i = [X_i; Y_i; theta_i; phi_i; psi_i];

L_i = 1/2*m*(X_dot^2 + Y_dot^2 + R^2*sin(theta)^2*theta_dot^2) + ...
    1/2*(I_A*(psi_dot - phi_dot*sin(theta))^2 + ...
         I_T*(theta_dot^2 + phi_dot^2*cos(theta)^2)) - ...
    m*g*R*(cos(theta_i)+cos(theta_j))/2 - alpha*(z_i+z_j)/2 + ...
    (F_i*q_i+F_j*q_j)/2;

% Precompute the term that is multiplied to D_2 L(q_{j-1},q_j,...)
partial_z = (1 + h*diff(L_j,z_j))/(1 - h*diff(L_i,z_j));

% Write system of equations in the form F(q) = 0, which with the
% constraint equations (1) and (2) corresponds to the whole system to be
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

save sym_equations.mat q_k eq_theta eq_phi eq_psi Omega_d_1 Omega_d_2;
