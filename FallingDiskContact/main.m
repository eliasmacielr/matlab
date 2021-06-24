% Elias Maciel, National University of Asuncion

clearvars
close all
clc

% TODO: check if discrete_equations exists and if not, execute
% discrete_equations.m
load discrete_equations
%% Simulation parameters and Integration

R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;

% Create disk struct
disk.R = R;
disk.m = m;
disk.I_A = I_A;
disk.I_T = I_T;

g = 9.81;

t0 = 0;
tf = 5;
T = tf - t0;
h = 0.1;
N = int32(T/h) + 1;
tol = 1e-6;

S = vpa(subs([eq_theta; eq_phi; eq_psi; Omega_d_1; Omega_d_2]));

% Initial conditions
% theta0 = 20*(pi/180);
% phidot0 = -0.15*(2*pi);
% psidot0 = ((I_T-I_A-m*R^2)*sin(theta0)*phidot0^2-m*g*R)/((I_A+m*R^2)*tan(theta0)*phidot0);
q0 = [0; 0; 0; 0; 0];
qdot0 = [0; 0; 0; 0; 0];
% Dissipation parameter
alpha = 0;
% External forces
F = {@(t) 0, @(t) 0, @(t) 0, @(t) 0, @(t) 1/2};

dimq = numel(q0);

q = zeros(dimq, N);
q(:,1) = vpa(q0);
% Get q(:,2) using the disk's differential equations
fixed_q0 = [0 0 0 0 0];
fixed_qdot0 = [0 0 0 0 0];
[~,y] = solve_diskODEs(disk,q0,fixed_q0,qdot0,fixed_qdot0,alpha,F,t0,...
    t0+h,h/2,tol);
q(:,2) = vpa(transpose(y(end,1:dimq)));

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

    F_X_i     = F{1}(tj-h);
    F_Y_i     = F{2}(tj-h);
    F_theta_i = F{3}(tj-h);
    F_phi_i   = F{4}(tj-h);
    F_psi_i   = F{5}(tj-h);

    F_X_j     = F{1}(tj);
    F_Y_j     = F{2}(tj);
    F_theta_j = F{3}(tj);
    F_phi_j   = F{4}(tj);
    F_psi_j   = F{5}(tj);

    [q(:,j+1), i] = newton_n_dim(q(:,j), q_k, subs(S), tol, 10);
%     fprintf("%d\n", i);
end

t = t0:h:tf;

% Get coordinates
X = q(1,:);
Y = q(2,:);
theta = q(3,:);
phi = q(4,:);
psi = q(5,:);

% Compute velocities from coordinates
Xdot = [y(1,6),diff(X)/h];
Ydot = [y(1,7),diff(Y)/h];
thetadot = [y(1,8),diff(theta)/h];
phidot = [y(1,9),diff(phi)/h];
psidot = [y(1,10),diff(psi)/h];

E = 1/2*m*(Xdot.^2 + Ydot.^2 + R^2*sin(theta).*thetadot.^2) + ...
    1/2*(I_A*(psidot - phidot.*sin(theta)).^2 + ...
        I_T*(thetadot.^2 + phidot.^2.*(cos(theta).^2))) + ...
    m*g*R*cos(theta);

%% Get reference solution
h_ref = h / 10;

[t_ref,y] = solve_diskODEs(disk,q0,fixed_q0,qdot0,fixed_qdot0,alpha,...
    F,t0,tf,h_ref,tol);

% Get coordinates
X_ref = y(:,1);
Y_ref = y(:,2);
theta_ref = y(:,3);
phi_ref = y(:,4);
psi_ref = y(:,5);

% Compute velocities from coordinates
Xdot_ref = [y(1,6);diff(X_ref)/h_ref];
Ydot_ref = [y(1,7);diff(Y_ref)/h_ref];
thetadot_ref = [y(1,8);diff(theta_ref)/h_ref];
phidot_ref = [y(1,9);diff(phi_ref)/h_ref];
psidot_ref = [y(1,10);diff(psi_ref)/h_ref];

E_ref = 1/2*m*(Xdot_ref.^2 + Ydot_ref.^2 + R^2*sin(theta_ref).*thetadot_ref.^2) + ...
    1/2*(I_A*(psidot_ref - phidot_ref.*sin(theta_ref)).^2 + ...
        I_T*(thetadot_ref.^2 + phidot_ref.^2.*(cos(theta_ref).^2))) + ...
    m*g*R*cos(theta_ref);

%% Animation and state space portraits
% animate_rolling_disk(q(1,:),q(2,:),q(3,:),q(4,:),q(5,:),R,h);

% TODO: repeat this same plot instructions for all coordinates and
% velocities and the energy function
figure
plot(t,X,'--x','MarkerIndices',1:10:length(X))
hold on
plot(t_ref,X_ref)
xlim([t0 tf])
legend('Contacto (orden 2)','Referencia')
xlabel('Tiempo (s)')
ylabel('X (m)')

figure
plot(t,E,'--x','MarkerIndices',1:10:length(X))
hold on
plot(t_ref,E_ref)
xlim([t0 tf])
legend('Contacto (orden 2)','Referencia')
xlabel('Tiempo (s)')
ylabel('Energía (J)')

% savefig(strcat('resultados-',num2str(tol),'-',num2str(T),'s.fig'));

%% Save last simulation results
% save(strcat('res-contact2-h',num2str(h),'-alpha',num2str(alpha),'.mat'),'t0','tf','h','q','y0');
