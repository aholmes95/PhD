function [returnDist] = findSSDistME(a,b,g,N)
% a=0.25;
% b=0.1;
% g=0.15;
% N=250;
% syms a b g
[A, states] = solveMEMat(a,b,g,N);
% disp('HELLO')
numStates =  length(states);
% A = A';
% sum(A)
A = A(2:numStates,2:numStates);
% [U,S,V] = svds(A,1,'smallestnz');
[U,S] = eigs(A,1,'smallestabs');
% [U,S,V] = eigs(A);
v = abs(U)/sum(abs(U));
% infdist = states(:,2)+2*states(:,3);
% numInf = sum(v.*infdist)
state1 = sum(v.*states(2:end,1));
state2 = sum(v.*states(2:end,2));
state3 = sum(v.*states(2:end,3));
returnDist = [state1, state2, state3];
returnDist = returnDist/sum(returnDist);