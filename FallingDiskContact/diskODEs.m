function res = diskODEs(t,y,yp)
% y := (q,q')

% TODO: obtener estos valores desde el archivo de simulacion --------------
R = 0.5;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.81;

alpha = 0; % dissipation parameter

FX = 0*t;
FY = 0*t;
Ftheta = 0*t;
Fphi = 0*t;
Fpsi = 0*t;
% -------------------------------------------------------------------------

theta = y(3); phi = y(4);
DX = yp(1); DY = yp(2); Dtheta = yp(3); Dphi = yp(4); Dpsi = yp(5);
DDX = yp(6); DDY = yp(7); DDtheta = yp(8); DDphi = yp(9); DDpsi = yp(10);

lambda_1 = m*DDX - FX + alpha*m*DX;
lambda_2 = m*DDY - FY + alpha*m*DY;

res = [
    y(6) - DX;
    y(7) - DY;
    y(8) - Dtheta;
    y(9) - Dphi;
    y(10) - Dpsi;

    m*R^2*2*sin(theta)*cos(theta)*Dtheta^2 + ...
    (m*R^2*sin(theta)^2 + I_T)*DDtheta - ...
    m*R^2*Dtheta^2*sin(theta)*cos(theta) + ...
    I_A*(Dpsi - Dphi*sin(theta))*Dphi*cos(theta) + ...
    I_T*Dphi^2*cos(theta)*sin(theta) - m*g*R*sin(theta) - Ftheta + ...
    (m*R^2*sin(theta)^2*Dtheta + I_T*Dtheta)*alpha - ...
    lambda_1*R*cos(theta)*sin(phi) + lambda_2*R*cos(theta)*cos(phi);

    -I_A*((DDpsi - DDphi*sin(theta) - Dphi*cos(theta)*Dtheta)*sin(theta) + ...
    (Dpsi - Dphi*sin(theta))*cos(theta)*Dtheta) + ...
    I_T*(DDphi*cos(theta)^2 - 2*Dphi*cos(theta)*sin(theta)*Dtheta) - ...
    Fphi + alpha*(-I_A*(Dpsi - Dphi*sin(theta))*sin(theta) + ...
    I_T*Dphi*cos(theta)^2) - ...
    lambda_1*R*sin(theta)*cos(phi) - lambda_2*R*sin(theta)*sin(phi);

    I_A*(DDpsi - DDphi*sin(theta) - Dphi*cos(theta)*Dtheta) - Fpsi + ...
    alpha*(I_A*Dpsi - I_A*Dphi*sin(theta)) + ...
    lambda_1*R*cos(phi) + lambda_2*R*sin(phi);

    DX + R*cos(theta)*sin(phi)*Dtheta + R*sin(theta)*cos(phi)*Dphi - ...
    R*cos(phi)*Dpsi;

    DY - R*cos(theta)*cos(phi)*Dtheta + R*sin(theta)*sin(phi)*Dphi - ...
    R*sin(phi)*Dpsi;
    ];
