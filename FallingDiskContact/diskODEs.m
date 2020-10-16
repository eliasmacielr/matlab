function res = diskODEs(t,y,yp)
% y := (q,q')

% TODO: obtener estos valores desde el archivo de simulacion --------------
R = 1;
m = 5;
I_A = 1/2*m*R^2;
I_T = 1/4*m*R^2;
g = 9.8;
alpha = 0.1;

fx = 0*t;
fy = 0*t;
ftheta = 0*t;
fphi = 0*t;
fpsi = 0*t;
% -------------------------------------------------------------------------

theta = y(3); phi = y(4);
DX = yp(1); DY = yp(2); Dtheta = yp(3); Dphi = yp(4); Dpsi = yp(5);
DDX = yp(6); DDY = yp(7); DDtheta = yp(8); DDphi = yp(9); DDpsi = yp(10);

lambda_1 = m*DDX - fx + m*DX*alpha;
lambda_2 = m*DDY - fy + m*DY*alpha;

res = [
    y(6) - DX;
    y(7) - DY;
    y(8) - Dtheta;
    y(9) - Dphi;
    y(10) - Dpsi;

    m*R^2*2*sin(theta)*cos(theta)*Dtheta^2 + ...
    (m*R*sin(theta)^2 + I_T)*DDtheta - ...
    m*R^2*Dtheta^2*sin(theta)*cos(theta) + ...
    I_A*(Dpsi - Dphi*sin(theta))*Dphi*cos(theta) + ...
    I_T*Dphi^2*cos(theta)*sin(theta) - m*g*R*sin(theta) + ftheta + ...
    (m*R^2*sin(theta)^2*Dtheta + I_T*Dtheta)*alpha - ...
    lambda_1*R*cos(theta)*sin(theta) + lambda_2*R*cos(theta)*cos(phi);

    I_A*((DDpsi - DDphi*sin(theta) - Dphi*cos(theta))*(-sin(theta)) - ...
    (Dpsi - Dphi*sin(theta))*cos(theta)) + ...
    I_T*(DDpsi*cos(theta)^2 - Dphi*2*cos(theta)*sin(theta)*Dtheta) - ...
    fphi + (I_A*(Dpsi - Dphi*sin(theta))*(-sin(theta))*(-sin(theta)) + ...
    I_T*Dphi*cos(theta)^2)*alpha - ...
    lambda_1*R*sin(theta)*cos(phi) - lambda_2*R*sin(theta)*sin(phi);

    I_A*(DDpsi - DDphi*sin(theta) - Dphi*cos(theta)) - fpsi + ...
    I_A*(Dpsi - Dphi*sin(theta))*alpha + ...
    lambda_1*R*cos(theta) + lambda_2*R*sin(theta);

    DX - R*cos(phi)*Dpsi + R*sin(theta)*cos(phi)*Dphi + ...
    R*cos(theta)*sin(phi)*Dtheta;

    DY - R*sin(phi)*Dpsi + R*sin(theta)*sin(phi)*Dphi - ...
    R*cos(theta)*cos(phi)*Dtheta;
    ];
