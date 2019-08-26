function [sol, val]=gaDemo1Eval(sol,~)
x=sol(1);
val=x+10*sin(5*x)+7*cos(4*x);
