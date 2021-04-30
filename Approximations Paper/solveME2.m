function [A, states] = solveME2(a,b,g,N)
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

parfor i=1:numStates
    for j=1:numStates
        m1old = states(j,1);
        m2old = states(j,2);
        m3old = states(j,3);
        m1new = states(i,1);
        m2new = states(i,2);
        m3new = states(i,3);
        old = [m1old m2old m3old];
        new = [m1new m2new m3new];
        if (m1new==m1old-1) && (m2new == m2old+1)
            A(i,j) = a/N*(m2old+2*m3old)*(m1old);
        elseif (m1new==m1old+1) && (m2new == m2old-1)
            A(i,j) = g*(m2old);
        elseif (m2new==m2old-1) && (m3new == m3old+1)
            A(i,j) = (a/(2*N)*(m2old+2*m3old)+b)*m2old;
        elseif (m2new==m2old+1) && (m3new == m3old-1)
            A(i,j) = 2*g*m3old;
        elseif (m1new==m1old) && (m2new==m2old)
            A(i,j) = -1*(a/N*(m2old+2*m3old)*(m1old)+g*(m2old)+(a/(2*N)*(m2old+2*m3old)+b)*m2old+2*g*m3old);
        end
    end
end
% A(1:numStates,:)=A(1:numStates,:)';
% A(end,:) = ones(numStates,1);
% solV = zeros(numStates,1);
% solV = [solV; 1];
% x = A\solV