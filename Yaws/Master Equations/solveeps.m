function myFunc = solveeps(eps1)

[a,A,c,P] = master(4, eps1, 160);
infStead = 0;
for i=1:length(P)
    infStead = infStead + P(i)*A(i,2);
end

myFunc = infStead/6000 - eps1;
end

