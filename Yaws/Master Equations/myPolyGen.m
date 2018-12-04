function y = myPolyGen(alphaa)
[Q, A, initCond, P] = master(4, alphaa, 160);
y = 0;
for i=1:15
    y = y + A(i,2)*P(i);
end

y = y*1500;
y = y-145.5729683888405;