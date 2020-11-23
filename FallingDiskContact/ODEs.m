syms X(t) Y(t) theta(t) phi(t) psi(t) DX(t) DY(t) Dtheta(t) Dphi(t) Dpsi(t) lambda1 lambda2
syms m g R I_A I_T alpha
syms F_X(t) F_Y(t) F_theta(t) F_phi(t) F_psi(t)
assume((m > 0) & (g > 0) & (R > 0) & (I_A > 0) & (I_T > 0) & (alpha > 0))

fun1 = diff(m*DX) - F_X + alpha*m*DX - lambda1;
fun2 = diff(m*DY) - F_Y + alpha*m*DY - lambda2;
fun3 = diff((m*R^2*sin(theta)^2+I_T)*Dtheta) - m*R^2*sin(theta)*cos(theta)*Dtheta^2 + I_A*(Dpsi-Dphi*sin(theta))*Dphi*cos(theta) + I_T*Dphi^2*cos(theta)*sin(theta) - m*g*R*sin(theta) - F_theta + alpha*(m*R^2*sin(theta)^2+I_T)*Dtheta - lambda1*R*cos(theta)*sin(phi) + lambda2*R*cos(theta)*cos(phi);
partialLphi = -I_A*(Dpsi - Dphi*sin(theta))*sin(theta) + I_T*Dphi*cos(theta)^2;
fun4 = diff(partialLphi) - F_phi + alpha*partialLphi - lambda1*R*sin(theta)*cos(theta) + lambda2*R*sin(theta)*sin(phi);
fun5 = diff(I_A*Dpsi - I_A*Dphi*sin(theta)) - F_psi + alpha*(I_A*Dpsi - I_A*Dphi*sin(theta)) + lambda1*R*cos(phi) + lambda2*R*sin(phi);
fun6 = DX + R*cos(theta)*sin(theta) + R*sin(theta)*cos(theta)*Dphi - R*cos(phi)*Dpsi;
fun7 = DY + R*cos(theta)*cos(phi)*Dtheta + R*sin(theta)*sin(phi)*Dphi - R*sin(phi)*Dpsi;
