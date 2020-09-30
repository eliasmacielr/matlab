function [x,i] = newton_n_dim(x0,sym_vars,sym_eqns,tol,maxit)
%% newton method for solving a system of n nonlinear equations for n variables
% Given n equations, the function performs the newton method, converging
% to the exact solution.
%
%input:     tol:        maximum tolerable RSS of errors in output vector
%           x0:         row vector of initial estimate
%           sym_vars:   row vector of n symbolic variables
%           sym_eqns:   column vector of >=n symbolic equations
%
%output:    solution:               row vector of solution.
%
%assumptions:   1.  Input sym functions are differentiable
%               2.  Convergence is dependent on the functions.
%                  -check convergence constraints.
%                    -http://en.wikipedia.org/wiki/Newton's_method
%%   Example:
%        syms a b
%        F1 = a - 15;             %(15 = a)
%        F2 = b^2 - 10;           %(10 = b^2)
%        tol = .01;
%        x0 = [10;1];
%with n equations and n unknowns:
%        sol = newton_n_dim(tol,x0,[a;b],[F1;F2]);
%Kyle J. Drerup
%Ohio University EECS
%11-9-2010
%Elias Maciel (adaptation)
%National University of Asuncion
%2020-09-29
%%  the code...
J = jacobian(sym_eqns,sym_vars);
x = x0;

i = 1;
while 1
    F = subs(sym_eqns,sym_vars,x);
    F_prime = subs(J,sym_vars,x);
    if ~isnumeric(F_prime)
        F_prime = eval(F_prime);
    end
    d_x = F_prime \ F;
    x = x - d_x;
    if sqrt(sum(d_x.^2)) <= tol
        break
    end
    if i >= maxit
        error("Error: Maximum number of iterations exceeded.")
    end
    i = i + 1;
end
end

% You're free to use and modify this code.  Sell it if you can.
