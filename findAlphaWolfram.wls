#!/usr/bin/env wolframscript
(* ::Package:: *)

eps = ToExpression[$ScriptCommandLine[[2]]];
b = ToExpression[$ScriptCommandLine[[3]]];
g = ToExpression[$ScriptCommandLine[[4]]];
NN = ToExpression[$ScriptCommandLine[[5]]];

Ax = a*(2 - 2 x - y) x;
Bx = g*y;
Cx = (a/2 (2 - 2 x - y) + b) y;
Dx = 2 g (1 - x - y);
f1 = -(Ax - Bx);
f2 = -(Bx - Ax + Cx - Dx);
D11  = Ax + Bx;
D12 = -Ax - Bx;
D21 = -Ax - Bx;
D22 = Ax + Bx + Cx + Dx;
dD11 = D[D11, x];
dD12 = D[D12, y];
dD21 = D[D21, x];
dD22 = D[D22, y];
x = g^2/((eps + g)^2 + b*eps);
y = 2 g*eps/((eps + g)^2 + b*eps);
alpha1 = f1 - (1/(2 NN)) (dD11 + dD12);
alpha2 = f2 - (1/(2 NN)) (dD21 + dD22);
a1=NumberForm[Solve[alpha1 == alpha2, {a}], 12]
mya = {a}/.a1[[1]][[1]];
Export["myFPalpha.csv",mya,"CSV"];


