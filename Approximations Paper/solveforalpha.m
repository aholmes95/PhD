function [myalpha, fval] = solveforalpha(e,b,g,N)
N=100;
fun = @(a) solvefun(a,e,b,g,N);
% myalpha = fzero(fun,[0 1]);
options = optimset('MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-14,'TolX',1e-14);
[myalpha,fval,exitflag,output]  = fminbnd(fun,0,1,options);


