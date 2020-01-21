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

%eq_X = simplify(diff(L_1,X_k) + diff(L_2,X_k)) == lambda_1;
eq_X = diff(L_1,X_k) + diff(L_2,X_k) == lambda_1;
eq_Y = diff(L_1,Y_k) + diff(L_2,Y_k) == lambda_2;
eq_theta = diff(L_1,theta_k) + diff(L_2,theta_k) == lambda_1*R*cos((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2) - lambda_2*R*cos((theta_k+theta_l)/2)*cos((phi_k+phi_l)/2);
eq_phi = diff(L_1,phi_k) + diff(L_2,phi_k) == lambda_1*R*sin((theta_k+theta_l)/2)*cos((phi_k+phi_k)/2) + lambda_2*R*sin((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2);
eq_psi = diff(L_1,psi_k) + diff(L_2,psi_k) == -lambda_1*R*cos((phi_k+phi_l)/2) - lambda_2*R*sin((phi_k+phi_l)/2);

% Constraints
omega_1 = ((X_l-X_k)/h) + R*cos((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2)*((theta_l-theta_k)/h) + R*sin((theta_k+theta_l)/2)*cos((phi_k+phi_l)/2)*((phi_l-phi_k)/h) - R*cos((phi_k+phi_l)/2)*((psi_l-psi_k)/h);
omega_2 = ((Y_l-X_k)/h) - R*cos((theta_k+theta_l)/2)*cos((phi_k+phi_l)/2)*((theta_l-theta_k)/h) + R*sin((theta_k+theta_l)/2)*sin((phi_k+phi_l)/2)*((phi_l-phi_k)/h) - R*sin((phi_k-phi_l)/2)*((psi_l-psi_k)/h);

%% Integrator


%% Lagrangian function (was going to use this instead of writing the complete equations)
function Lagrangian = L(X,Y,theta,phi,psi,dX,dY,dtheta,dphi,dpsi)
    Lagrangian = (1/2)*m*(dX^2 + dY^2 + R*(sin(theta))^2*dtheta^2) + ...
        (1/2)*(I_A*(dpsi - dphi*sin(theta))^2 + I_T*(dtheta^2 + dphi^2*(cos(theta))^2)) - ...
        m*g*R*cos(theta);
end
