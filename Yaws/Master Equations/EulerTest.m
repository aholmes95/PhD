clear;
clc;
houseMat=createHouseMat;
alpha = 0.1653;
numPeop = 5930;
N=1000;
t=zeros(1,N);
tmin = 0;
tmax = 500;
t(1)=tmin;
h = (tmax-tmin)/N;
numInf1 = 160;
for ii=1:1500
    n = houseMat(ii);
    matArray{ii} = zeros(0.5*(n+1)*(n+2), N);
    [Q, A, initCond] = master(n, (alpha*numInf1/numPeop), numInf1);
%     [Q, A, initCond] = master(n, 0.004, numInf1);
    matArray{ii}(:,1) = initCond;
end
yplot(1) = numInf1;
for n=1:N
    n
    for ii=1:1500
        if n==1
            numInfExt = numInf1;
        else
            numInfExt = yplot(n-1);
        end
        nn = houseMat(ii,1);
        [Q, A, initCond] = master(nn, (alpha*sum(numInfExt)/numPeop), numInf1);
%         [Q, A, initCond] = master(nn, 0.004, numInf1);
        k1 = h*Q*matArray{ii}(:,n);
        k2 = h*Q*(matArray{ii}(:,n) + k1/2);
        k3 = h*Q*(matArray{ii}(:,n) + k2/2);
        k4 = h*Q*(matArray{ii}(:,n) + k3);
        matArray{ii}(:,n+1) = matArray{ii}(:,n) + 1/6*(k1+2*k2+2*k3+k4);
        numInf(ii) = sum(A(:,2).*matArray{ii}(:,n));
    end
    yplot(n+1) = sum(numInf);
    t(n+1)=t(n)+h;
    
end
plot(t,yplot, '--g')
endInf = yplot(end);

        
    
