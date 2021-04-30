#!/usr/bin/env wolframscript
(* ::Package:: *)

params = Import["myParams.csv"];
a=params[[1]][[1]];
b=params[[2]][[1]];
g=params[[3]][[1]];
NN=params[[4]][[1]];
Ax = a*(2-2x-y)x;
Bx = g*y;
Cx = (a/2(2-2x-y)+b)y;
Dx = 2g(1-x-y);
f1=-(Ax-Bx);
f2=-(Bx-Ax+Cx-Dx);
D11  = Ax+Bx;
D12 = -Ax-Bx;
D21 = -Ax-Bx;
D22 = Ax+Bx+Cx+Dx;
dD11 = D[D11, x];
dD12 = D[D12, y];
dD21 = D[D21, x];
dD22 = D[D22, y];
alpha1 = f1-(1/(2NN))(dD11+dD12);
alpha2 = f2-(1/(2NN))(dD21+dD22);
alpha1 = Collect[Expand[alpha1],{x[t],y[t],x[t]^2,y[t]^2,x[t]*y[t]}];
alpha2 = Collect[Expand[alpha2],{x[t],y[t],x[t]^2,y[t]^2,x[t]*y[t]}];
roots = Solve[{alpha1==0,alpha2==0},{x,y}];
myRoot = {x,y}/.roots[[2]];
a1=myRoot[[1]];
a2=myRoot[[2]];
dalpha11=D[alpha1,x];
dalpha12=D[alpha1,y];
dalpha21=D[alpha2,x];
dalpha22=D[alpha2,y];
f11[x_,y_] = Evaluate[dalpha11];
f12[x_,y_] = Evaluate[dalpha12];
f21[x_,y_] = Evaluate[dalpha21];
f22[x_,y_] = Evaluate[dalpha22];
myEvals=Eigenvalues[{{f11[a1,a2],f12[a1,a2]},{f21[a1,a2],f22[a1,a2]}}];
myEvecs=Eigenvectors[{{f11[a1,a2],f12[a1,a2]},{f21[a1,a2],f22[a1,a2]}}];
evec1 = myEvecs[[1]];
evec2=myEvecs[[2]];
Export["myEvec1.csv",N[evec1,16],"CSV"];
Export["myEvec2.csv",N[evec2,16],"CSV"];
Export["mySS.csv",{N[a1,16],N[a2,16]},"CSV"];



eqn1 = -a*(p1+2p2)*p0+g*p1;
eqn2 = a*(p1+2p2)*p0-((a/2)*(p1+2p2)+b+g)*p1+2g*p2;
eqn3 =((a/2)*(p1+2p2)+b)*p1-2g*p2;
sols=NSolve[{eqn1==0,eqn2==0,eqn3==0,p0+p1+p2==1,p0>=0,p1>=0,p2>=0},{p0,p1,p2},Reals];
p11=p1/.sols[[2]];
p22=p2/.sols[[2]];
Export["mySSSHH.csv",{N[p11,16],N[p22,16]},"CSV"]



mmomdot = -2*a*m+g*n+2*a/NN*m2+a/NN*mn;
nmomdot = (2*a-2*g)*m+(-a-b-3*g)*n-2*a/NN*m2+a/(2*NN)*n2 + 2*g*NN;
m2momdot = 2*a*m+g*n+(-4*a-2*a/NN)*m2+(2*g-a/NN)*mn+2*a/NN*m2n+4*a/NN*m3;
n2momdot = (2*a-2*g)*m+(4*a-4*g-2*a/NN)*mn+(-2*a-2*b-6*g-a/(2*NN))*n2-2*a/NN*m2-4*a/NN*m2n +a/NN*n3+2*g*NN+(a+b-g+4*g*NN)*n;
mnmomdot = -g*n+g*n2+(-3*a-b-3*g+a/NN)*mn+(2*a-2*g+2*a/NN)*m2+2*a/NN*m2n-2*a/NN*m3+3*a/(2*NN)*mn2+(2*g*NN-2*a)*m;
clos1 = m3 - 3*m2*m + 2*m^3;
clos2 = n3 - 3*n2*n + 2*n^3;
clos3 = m2n - m2*n - 2*mn*m + 2*m^2*n;
clos4 = mn2 - m*n2 - 2*mn*n + 2*m*n^2;
sol = Solve[{mmomdot == 0, nmomdot == 0, m2momdot == 0, n2momdot == 0, mnmomdot == 0, clos1 == 0, clos2 == 0, clos3 == 0, clos4 == 0, m > 0, n > 0,m2 >= 0, m3 >= 0, n2 >= 0, n3 >= 0, mn >= 0, m2n >= 0, mn2 >= 0, 
m + n <= NN}, {m, m2, m3, n, n2, n3, mn, m2n, mn2}, Reals] // RootReduce // N[#, 10] &;
msol = m/.sol[[1]];
nsol = n/.sol[[1]];
Export["mySSMC.csv",{N[msol,16],N[nsol,16]},"CSV"]

