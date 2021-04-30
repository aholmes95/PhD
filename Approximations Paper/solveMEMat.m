function [A, states] = solveMEMat(a,b,g,N)
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
% disp('Yay')
m1old = sparse(states(:,1));
m2old = sparse(states(:,2));
m3old = sparse(states(:,3));
% m1old = states(:,1);
% m2old = states(:,2);
% m3old = states(:,3);
% disp('Yayy')
m1new = sparse(states(:,1));
m2new = sparse(states(:,2));
m3new = sparse(states(:,3));
% m1new = states(:,1);
% m2new = states(:,2);
% m3new = states(:,3);
% disp('Yayyy')
v11 = a/N*(m2old+2*m3old).*(m1old);
v12 = g*(m2old);
v13 = (a/(2*N)*(m2old+2*m3old)+b).*m2old;
v14 = 2*g*m3old;
v15 = -1*(a/N*(m2old+2*m3old).*(m1old)+g*(m2old)+(a/(2*N)*(m2old+2*m3old)+b).*m2old+2*g*m3old);
% disp('Size of v11');
% disp(size(v11));
% disp('Yayyyy')
% minew = {m1new,m2new,m3new};
% miold = {m1old,m2old,m3old};
% inds = [1,1,2,2,3,3];
% 
% ai = {1,2,3,4,5,6};
% parfor i=1:6
%     mnew = minew{inds(i)};
%     mold = miold{inds(i)};
%     ai{i} = (mnew == mold.'+(-1)^(i+1));
% end
    

a1p = (m1new == m1old.'+1);
% disp('a1p done');
a1m = (m1new == m1old.'-1);
a1s = (m1new == m1old');
clear m1new m1old;
% disp('a1m done');

a2p = (m2new == m2old.'+1);
% disp('a2p done');
a2m = (m2new == m2old.'-1);
a2s = (m2new == m2old');
clear m2new m2old;
% disp('a2m done');

a3p = (m3new == m3old.'+1);
% disp('a3p done');
a3m = (m3new == m3old.'-1);
a3s = (m3new == m3old');
clear m3new m3old;
% disp('a3m done');

% disp('Yayyyy')
v21 = (a1m & a2p & a3s);
% disp('Size of v21');
% disp(size(v21));
clear a1m;
% disp('1 done');
v22 = (a1p & a2m & a3s);
clear a1p;
% disp('2 done');
v23 = (a1s & a2m & a3p);
clear a2m a3p;
% disp('3 done');
v24 = (a2p & a3m);
clear a2p a3m;
% disp('4 done');

% v21 = (ai{2} & ai{3});
% v22 = (ai{1} & ai{4});
% v23 = (ai{4} & ai{5});
% v24 = (ai{3} & ai{6});
% disp('Yayyyyy')
v1 = v11'.*v21;
clear v11 v21;
% disp('v1')
v2 = v12'.*v22;
clear v12 v22;
% disp('v2')
v3 = v13'.*v23;
clear v13 v23;
% disp('v3')
v4 = v14'.*v24;
clear v14 v24;
% disp('v4')
A = v1+v2+v3+v4;
clear v1 v2 v3 v4;
% disp('all')
A = spdiags(v15,0,A);
nnz(A)
numel(A)
nnz(A)/numel(A)
% A = A';
% disp('Yayyyyyy')