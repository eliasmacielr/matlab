clear all
close all
clc

%% Symbols
% q_{k-1}
syms X_j Y_j theta_j phi_j psi_j
% q_k
syms X_k Y_k theta_k phi_k psi_k
% q_{k+1}
syms X_l Y_l theta_l phi_l psi_l
% constants
syms m R I_A I_T g
% Lagrange multipliers
syms lambda_1 lambda_2
% step size
syms h

%% Expressions for the Discrete (unconstrained) Lagrangian
L_1 = h*( ...
    (1/2)*m*(((X_l-X_k)/h)^2 + ((Y_l-Y_k)/h)^2 + R*(sin((theta_k+theta_l)/2))^2*((theta_l-theta_k)/h)^2) + ...
    (1/2)*(I_A*(((psi_l-psi_k)/h) - ((phi_l-phi_k)/h)*sin((theta_k+theta_l)/2))^2 + I_T*(((theta_l-theta_k)/h)^2 + ((phi_l-phi_k)/h)^2*(cos((theta_k+theta_l)/2))^2)) - ...
    m*g*R*cos((theta_k+theta_l)/2) ...
    );

L_2 = h*( ...
    (1/2)*m*(((X_k-X_j)/h)^2 + ((Y_k-Y_j)/h)^2 + R*(sin((theta_j+theta_k)/2))^2*((theta_k-theta_j)/h)^2) + ...
    (1/2)*(I_A*(((psi_k-psi_j)/h) - ((phi_k-phi_j)/h)*sin((theta_j+theta_k)/2))^2 + I_T*(((theta_k-theta_j)/h)^2 + ((phi_k-phi_j)/h)^2*(cos((theta_j+theta_k)/2))^2)) - ...
    m*g*R*cos((theta_j+theta_k)/2) ...
    );

%% Disk equations
lambda_1 = diff(L_1,X_k) + diff(L_2,X_k);
lambda_2 = diff(L_1,Y_k) + diff(L_2,Y_k);
eq_theta = diff(L_1,theta_k) + diff(L_2,theta_k) - lambda_1*R*cos((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2) + lambda_2*R*cos((theta_k+theta_l)/2)*cos((phi_k+phi_l)/2) == 0;
eq_phi = diff(L_1,phi_k) + diff(L_2,phi_k) - lambda_1*R*sin((theta_k+theta_l)/2)*cos((phi_k+phi_k)/2) - lambda_2*R*sin((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2) == 0;
eq_psi = diff(L_1,psi_k) + diff(L_2,psi_k) + lambda_1*R*cos((phi_k+phi_l)/2) + lambda_2*R*sin((phi_k+phi_l)/2) == 0;

% Constraints
eq_X = ((X_l-X_k)/h) + R*cos((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2)*((theta_l-theta_k)/h) + R*sin((theta_k+theta_l)/2)*cos((phi_k+phi_l)/2)*((phi_l-phi_k)/h) - R*cos((phi_k+phi_l)/2)*((psi_l-psi_k)/h) == 0;
eq_Y = ((Y_l-X_k)/h) - R*cos((theta_k+theta_l)/2)*cos((phi_k+phi_l)/2)*((theta_l-theta_k)/h) + R*sin((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2)*((phi_l-phi_k)/h) - R*sin((phi_k-phi_l)/2)*((psi_l-psi_k)/h) == 0;

%% Initial conditions and Integration
T = 1;
h = 0.1;
N = T/h;

m = 0.5;
R = 0.15;
I_A = 2;
I_T = 1;
g = 9.8;

q = zeros(5, N);

X_j = 0;
Y_j = 0;
theta_j = 0;
phi_j = 0;
psi_j = 0;

X_k = 0;
Y_k = 0;
theta_k = 0;
phi_k = pi/6;
psi_k = pi/6;

sol = vpasolve( ...
    [subs(eq_theta) subs(eq_phi) subs(eq_psi) subs(eq_X) subs(eq_Y)], ...
    [theta_l phi_l psi_l X_l Y_l]);
q(:,1) = [sol.X_l; sol.Y_l; sol.theta_l; sol.phi_l; sol.psi_l];

X_j = X_k;
Y_j = Y_k;
theta_j = theta_k;
phi_j = phi_k;
psi_j = psi_k;

X_k = q(1,1);
Y_k = q(2,1);
theta_k = q(3,1);
phi_k = q(4,1);
psi_k = q(5,1);

sol = vpasolve( ...
    [subs(eq_theta) subs(eq_phi) subs(eq_psi) subs(eq_X) subs(eq_Y)], ...
    [theta_l phi_l psi_l X_l Y_l]);
q(:,2) = [sol.X_l; sol.Y_l; sol.theta_l; sol.phi_l; sol.psi_l];

for i = 2:N-1
    X_j = q(1,i-1);
    Y_j = q(2,i-1);
    theta_j = q(3,i-1);
    phi_j = q(4,i-1);
    psi_j = q(5,i-1);

    X_k = q(1,i);
    Y_k = q(2,i);
    theta_k = q(3,i);
    phi_k = q(4,i);
    psi_k = q(5,i);

    sol = vpasolve( ...
        [subs(eq_theta) subs(eq_phi) subs(eq_psi) subs(eq_X) subs(eq_Y)], ...
        [theta_l phi_l psi_l X_l Y_l]);
    q(:,i+1) = [sol.X_l; sol.Y_l; sol.theta_l; sol.phi_l; sol.psi_l];
end

%% Solve analytically
% syms X(t) Y(t) theta(t) phi(t) psi(t)
% 
% lambda_1 = diff(m*X, 2);
% lambda_2 = diff(m*Y, 2);
% ode1 = diff((m*R^2*(sin(theta))^2 + I_T)*diff(theta)) - m*R^2*diff(theta)^2*sin(theta)*cos(theta) + ...
%     I_A*(diff(psi) - diff(phi)*sin(theta))*diff(phi)*cos(theta) + I_T*diff(phi)^2*sin(theta)*cos(theta) - ...
%     m*g*R*sin(theta) == lambda_1*R*cos(theta)*sin(phi) - lambda_2*R*cos(theta)*cos(phi);
% ode2 = diff(-I_A*(diff(psi) - diff(phi)*sin(theta))*sin(theta) + I_T*diff(phi)*(cos(theta))^2) == ...
%     lambda_1*R*sin(theta)*cos(phi) + lambda_2*R*sin(theta)*sin(phi);
% ode3 = diff(I_A*(diff(psi) - diff(phi)*sin(theta))) == -lambda_1*R*cos(phi) - lambda_2*R*sin(phi);
% ode4 = diff(X) == R*cos(phi)*diff(psi) - R*sin(theta)*cos(phi)*diff(phi) - R*cos(theta)*sin(phi)*diff(theta);
% ode5 = diff(Y) == R*sin(phi)*diff(psi) - R*sin(theta)*sin(phi)*diff(phi) + R*cos(theta)*cos(phi)*diff(theta);
% odes = [ode1; ode2; ode3; ode4; ode5];
% 
% cond1 = X(0) == 0;
% cond2 = Y(0) == 0;
% cond3 = theta(0) == 0;
% cond4 = phi(0) == 0;
% cond5 = psi(0) == 0;
% conds = [cond1; cond2; cond3; cond4; cond5];
% 
% S = dsolve(odes,conds);
