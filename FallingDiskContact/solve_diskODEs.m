function [t,y] = solve_diskODEs(disk, q0, fixed_q0, qdot0, fixed_qdot0, ...
    alpha, F, t0, tf, h, tol)

g = 9.81;

R = disk.R;
m = disk.m;
I_A = disk.I_A;
I_T = disk.I_T;

dimq = numel(q0);

[y0,yp0] = decic(@disk_odes,0,[q0;qdot0],[fixed_q0,fixed_qdot0], ...
    [qdot0;zeros(dimq,1)],[fixed_qdot0,zeros(1,dimq)]);

[t,y] = ode15i(@disk_odes,t0:h:tf,y0,yp0,odeset('RelTol',tol));

    function res = disk_odes(t,y,yp)
    % y := (q,q')

    FX = F{1}(t);
    FY = F{2}(t);
    Ftheta = F{3}(t);
    Fphi = F{4}(t);
    Fpsi = F{5}(t);

    theta = y(3); phi = y(4);
    DX = yp(1); DY = yp(2); Dtheta = yp(3); Dphi = yp(4); Dpsi = yp(5);
    DDX = yp(6); DDY = yp(7); DDtheta = yp(8); DDphi = yp(9);
    DDpsi = yp(10);

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

        -I_A*((DDpsi - DDphi*sin(theta) - ...
        Dphi*cos(theta)*Dtheta)*sin(theta) + ...
        (Dpsi - Dphi*sin(theta))*cos(theta)*Dtheta) + ...
        I_T*(DDphi*cos(theta)^2 - ...
        2*Dphi*cos(theta)*sin(theta)*Dtheta) - Fphi + ...
        alpha*(-I_A*(Dpsi - Dphi*sin(theta))*sin(theta) + ...
        I_T*Dphi*cos(theta)^2) - lambda_1*R*sin(theta)*cos(phi) - ...
        lambda_2*R*sin(theta)*sin(phi);

        I_A*(DDpsi - DDphi*sin(theta) - Dphi*cos(theta)*Dtheta) - ...
        Fpsi + alpha*(I_A*Dpsi - I_A*Dphi*sin(theta)) + ...
        lambda_1*R*cos(phi) + lambda_2*R*sin(phi);

        DX + R*cos(theta)*sin(phi)*Dtheta + ...
        R*sin(theta)*cos(phi)*Dphi - R*cos(phi)*Dpsi;

        DY - R*cos(theta)*cos(phi)*Dtheta + ...
        R*sin(theta)*sin(phi)*Dphi - R*sin(phi)*Dpsi;
        ];
    end
end
