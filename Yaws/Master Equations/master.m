function [Q, A, initCond] = master(n, eps1, initInfNum, beta, delta, rho)
A = zeros((n+1)*(n+2)/2,3);
a=1;
delta1 = delta;
delta2 = delta;
lambda = 0.185;
initCond = zeros(n,1);
for s=0:n
    for i=0:n-s
        A(a,1) = s;
        A(a,2) = i;
        A(a,3) = n-s-i;
        if s==n-1 && i==1
            initCond(a) = initInfNum/1500;
        end
        if s==n
            initCond(a) = 1-initInfNum/1500;
        end
        a=a+1;
    end
end


numStates = 0.5*(n+1)*(n+2);       
Q = zeros(numStates);
%tic;
for j=1:numStates
    for k=1:numStates
        if n > 1
            if A(j,1) == A(k,1) && A(j,2) == A(k,2)
                Q(j,k) = -1*((eps1+beta*A(j,2)/(n-1))*A(j,1) + delta1*A(j,2) + delta2*A(j,3) + rho*(n-A(j,1)-A(j,2)) + lambda*A(j,2));
            else
                if A(j,1) == A(k,1) + 1 && A(j,2) == A(k,2) - 1
                    Q(j,k) = (eps1+beta*A(j,2)/(n-1))*(A(j,1));
                elseif A(j,1) == A(k,1) - 1 && A(j,2) == A(k,2)
                    Q(j,k) = delta2*(n-A(k,1)-A(k,2)+1);
                elseif A(j,1) == A(k,1) - 1 && A(j,2) == A(k,2) +1
                    Q(j,k) = delta1*A(j,2);
                elseif A(j,1) == A(k,1) && A(j,2) == A(k,2) - 1
                    Q(j,k) = rho*(n-A(k,1)-A(k,2)+1);
                elseif A(j,1) == A(k,1) && A(j,2) == A(k,2) + 1
                    Q(j,k) = lambda*A(j,2);
                else
                    Q(j,k) = 0;
                end
            end
        else
            if A(j,1) == A(k,1) && A(j,2) == A(k,2)
                Q(j,k) = -1*((eps1)*A(j,1) + delta1*A(j,2) + delta2*A(j,3) + rho*(n-A(j,1)-A(j,2)) + lambda*A(j,2));
            else
                if A(j,1) == A(k,1) + 1 && A(j,2) == A(k,2) - 1
                    Q(j,k) = (eps1)*(A(j,1));
                elseif A(j,1) == A(k,1) - 1 && A(j,2) == A(k,2)
                    Q(j,k) = delta2*(n-A(k,1)-A(k,2)+1);
                elseif A(j,1) == A(k,1) - 1 && A(j,2) == A(k,2) +1
                    Q(j,k) = delta1*A(j,2);
                elseif A(j,1) == A(k,1) && A(j,2) == A(k,2) - 1
                    Q(j,k) = rho*(n-A(k,1)-A(k,2)+1);
                elseif A(j,1) == A(k,1) && A(j,2) == A(k,2) + 1
                    Q(j,k) = lambda*A(j,2);
                else
                    Q(j,k) = 0;
                end
            end
        end
    end
end
%toc;
% [V, D, W] = eig(Q);
% pos = find(abs(diag(D))<0.00001);
% P = W(:,pos);
% P = P / sum(P); %eigenvector corresponding to zero eigenvalue
% zeroVec = zeros(numStates,1);
% ansVec = [zeroVec; 1];
Q=Q';
% Q(numStates+1,:) = ones(1, numStates);
% P=linsolve(Q, ansVec); %Solution to master equation at steady state
            
                
        