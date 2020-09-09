function [x,err,f] = nrndn(x0,tol,iprint)
%
% inputs:
% x0 = initial guess of x, column vector of length n
% tol = RMS tolerance on relative error on x, scalar
% iprint = value of 1 requests iteration information
% f(x) entered in the function "funkeval" at the bottom of this file
%
% outputs:
%
% x = converged solution, column vector of length n
% f = RMS value of f, scalar
% err = RMS value of relative error on x, scalar
%
maxit = 1000;
n = max(size(x0));
Residual = zeros(n,1);
Jacobian = zeros(n,n);
InvJ = zeros(n,n);
dx = zeros(n,1);
x = zeros(n,1);
xold = zeros(n,1);
xeval = zeros(n,1);
xp = zeros(2);
fp = zeros(n,n,2);
dxcon = zeros(n,1);
dxcon(1:n) = 0.01;
x = x0;
err = 100.0;
iter = 0;
while ( err > tol )
    for j = 1:1:n
        dx(j) = min(dxcon(j)*x(j),dxcon(j));
        i2dx(j) = 1.0/(2.0*dx(j));
    end
    Residual = funkeval(x);
    for j = 1:1:n
        for i = 1:1:n
            xeval(i) = x(i);
        end
        xp(1) = x(j) - dx(j);
        xp(2) = x(j) + dx(j);
        for k = 1:1:2
            xeval(j) = xp(k);
            fp(:,j,k) = funkeval(xeval);
        end
    end
    for i = 1:1:n
        for j = 1:1:n
            Jacobian(i,j) = i2dx(j)*( fp(i,j,2) - fp(i,j,1) );
        end
    end
    xold = x;
    invJ = inv(Jacobian);
    deltax = -invJ*Residual;
    for j = 1:1:n
        x(j) = xold(j) + deltax(j);
    end
    iter = iter +1;
    err = sqrt( sum(deltax.^2) /n ) ;
    f = sqrt(sum(Residual.*Residual)/n);
    if (iprint == 1)
        fprintf (1,'iter = %4i, err = %9.2e f = %9.2e \n ', iter, err, f);
    end
    if (iter > maxit)
        Residual
        error ('maximum number of iterations exceeded');
    end
end

function f = funkeval(x)
n = max(size(x));
f = zeros(n,1);
f(1) = x(1)^2 + x(2)^2 - 4;
f(2) = x(1)^2 - x(2) + 1;
