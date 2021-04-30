function kldiv = solvefun(a,e,b,g,N)
epsSS = [g^2/((e+g)^2+b*e), 2*g*e/((e+g)^2+b*e), (e^2+b*e)/((e+g)^2+b*e)];
retDist = findSSDistME(a,b,g,N);
kldiv = abs(CalKLDiv(epsSS,retDist));
% kldiv = sum(abs(epsSS-retDist));
end

