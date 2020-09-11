% TODO: write documentation and describe the problem and where it came from

clear var
close all
clc

%% Symbols
% q_{k-1}
syms X_j Y_j theta_j phi_j psi_j z_j
% q_k
syms X_k Y_k theta_k phi_k psi_k z_k
% q_{k+1}
syms X_l Y_l theta_l phi_l psi_l z_l
% constants
m = 1;
R = 0.15;
I_A = 2;
I_T = 1;
g = 9.8;
alpha = 0.1;
% "Control" forces
u_psi = 0;
% Lagrange multipliers
syms lambda_1 lambda_2
% step size
h = 0.1;

%% Disk discrete equations
% TODO: fix the system, there's something wrong, or the integrator does not
% work
lambda_1 = m*(-X_j + 2*X_k - X_l)/h^2 - m*alpha*(-X_j + X_k)/h;
lambda_2 = m*(-Y_j + 2*Y_k - Y_l)/h^2 - m*alpha*(-Y_j + Y_k)/h;
eq_theta = m*(R^2*(sin(theta_k)*cos(theta_k)*((theta_l-theta_k)/h)^2 ...
    - (sin(theta_k^2))^2*(theta_l-theta_k)/h^2 + ...
    (sin(theta_j))^2*(theta_k-theta_j)/h^2)) + ...
    I_A*((psi_l-psi_k)/h - ...
        (phi_l-phi_k)/h*sin(theta_k))*(-(phi_l-phi_k)/h*cos(theta_k)) + ...
    I_T*((theta_k-theta_l)/h^2 - ...
        ((phi_l-phi_k)/h)^2*cos(theta_k)*sin(theta_k) + ...
        (theta_k-theta_j)/h^2) + m*g*R*sin(theta_k) - ...
    lambda_1*R*cos(theta_k)*sin(phi_k) + ...
    lambda_2*R*cos(theta_k)*cos(phi_k);
eq_phi = I_A*(((psi_l-psi_k)/h - ...
    (phi_l-phi_k)/h*sin(theta_k))*sin(theta_k)/h - ...
    ((psi_k-psi_j)/h - (phi_k-phi_j)/h*sin(theta_j))*sin(theta_j)/h) + ...
    I_T*((cos(theta_k))^2*(phi_l-phi_k)/h^2 + (cos(theta_j))^2/h) - ...
    lambda_1*R*sin(theta_k)*cos(phi_k) - ...
    lambda_2*R*sin(theta_k)*sin(phi_k);
eq_psi = I_A*((-psi_l+psi_k)/h^2 + (phi_l-phi_k)/h*sin(theta_k) + ...
    (psi_k-psi_j)/h^2 - (phi_k-phi_j)/h^2*sin(theta_j)) + u_psi + ...
    lambda_1*R*cos(phi_k) + lambda_2*R*sin(phi_k);

% Constraints
eq_X = ((X_l-X_k)/h) + R*cos(theta_k)*sin(phi_k)*(theta_l-theta_k)/h + ...
    R*sin(theta_k)*cos(phi_k)*(phi_l-phi_k)/h - ...
    R*cos(phi_k)*(psi_l-psi_k)/h;
eq_Y = ((Y_l-X_k)/h) - R*cos(theta_k)*cos(phi_k)*(theta_l-theta_k)/h + ...
    R*sin(theta_k)*sin(phi_k)*(phi_l-phi_k)/h - ...
    R*sin(phi_k)*(psi_l-psi_k)/h;

%% Initial conditions and Integration
T = 1;
h = 0.1;
N = int32(T/h);

q = zeros(5, N);

q(:,1) = [0; 0; 0; 0; 0];
q(:,2) = [R*pi/12; 0; 0; 0; pi/12];

for k = 2:N-1
    X_j     = q(1,k-1);
    Y_j     = q(2,k-1);
    theta_j = q(3,k-1);
    phi_j   = q(4,k-1);
    psi_j   = q(5,k-1);

    X_k     = q(1,k);
    Y_k     = q(2,k);
    theta_k = q(3,k);
    phi_k   = q(4,k);
    psi_k   = q(5,k);

    sol_k = newton_n_dim(1e-4, q(:,k)', ...
        [X_l Y_l theta_l phi_l psi_l], ...
        [subs(eq_theta) subs(eq_phi) subs(eq_psi) subs(eq_X) ...
        subs(eq_Y)]', 1000);
    q(:,k+1) = sol_k';
end

%% Animation
animate_rolling_disk(q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),R,h);

%% Save last simulation results
% TODO: save q
