function [A, states] = solveME(a,b,g,N)
% a=0.1;
% b=0.1;
% g=0.15;
% N=250;
numStates = 0.5*(N+1)*(N+2);
A = sparse(numStates,numStates);
states = zeros(numStates,3);
counter = 1;
for i=1:N+1
    for j=1:N+2-i
        states(numStates+1-counter,:) = [i-1, j-1, N-i-j+2];
        counter = counter+1;
    end
end
% states
x = zeros(10*numStates,1);
y = zeros(10*numStates,1);
val = zeros(10*numStates,1);
endpos = 0;

m1old = states(:,1);
m2old = states(:,2);
m3old = states(:,3);

v11 = a/N*(m2old+2*m3old).*(m1old);
v12 = g*(m2old);
v13 = (a/(2*N)*(m2old+2*m3old)+b).*m2old;
v14 = 2*g*m3old;
v15 = -1*(a/N*(m2old+2*m3old).*(m1old)+g*(m2old)+(a/(2*N)*(m2old+2*m3old)+b).*m2old+2*g*m3old);
    
for i=1:numStates
    
    m1new = states(i,1);
    m2new = states(i,2);
    m3new = states(i,3);
    
    v21 = ((m1new==m1old-1) & (m2new == m2old+1));
    v22 = ((m1new==m1old+1) & (m2new == m2old-1));
    v23 = ((m2new==m2old-1) & (m3new == m3old+1));
    v24 = ((m2new==m2old+1) & (m3new == m3old-1));
    v25 = ((m1new==m1old) & (m2new==m2old));
    
    v1=v11.*v21;
    v2=v12.*v22;
    v3=v13.*v23;
    v4=v14.*v24;
    v5=v15.*v25;
    v = v1+v2+v3+v4+v5;
    
    vv = (v ~= 0);
    numelem = sum(vv);
    startpos = endpos + 1;
    endpos = startpos + numelem - 1;
    tempx = i;
    tempy = find(vv);
    tempval = v(vv);
    x(startpos:endpos) = tempx;
    y(startpos:endpos) = tempy(:);
    val(startpos:endpos) = tempval(:);
end
x = x(1:endpos);
y = y(1:endpos);
val = val(1:endpos);
A = sparse(x,y,val);