clear;
clc;
houseMat=createHouseMatMast;
numPeop = sum(houseMat(:,1));
% alpha = 0.1653;
% alpha = 0.17;
alpha = 0.164865777630342;
beta = 0.0516;
delta = 0.0513;
rho = 0.0165;
% rho = 0;
eps = 0.004;
N=200000;
t=zeros(1,N);
tmin = 0;
tmax = 700;
t(1)=tmin;
h = (tmax-tmin)/N;
numInf1 = 160;
for ii=1:13
    matArray{ii} = zeros(0.5*(ii+1)*(ii+2), N);
    [Q, A, initCond] = master(ii, (alpha*numInf1/numPeop), numInf1, beta, delta, rho);
%     [Q, A, initCond] = master(ii, 0.004, numInf1, beta, delta, rho);
    matArray{ii}(:,1) = initCond;
end
yplot(1) = numInf1;
for n=1:N
    n
    for ii=4
        if n==1
            numInfExt = numInf1;
        else
            numInfExt = newVal;
        end
        [Q, A, initCond] = master(ii, (alpha*sum(numInfExt)/numPeop), numInf1, beta, delta, rho);
%         [Q, A, initCond] = master(ii, 0.004, numInf1, beta, delta, rho);
        k1 = h*Q*matArray{ii}(:,n);
        k2 = h*Q*(matArray{ii}(:,n) + k1/2);
        k3 = h*Q*(matArray{ii}(:,n) + k2/2);
        k4 = h*Q*(matArray{ii}(:,n) + k3);
        matArray{ii}(:,n+1) = matArray{ii}(:,n) + 1/6*(k1+2*k2+2*k3+k4);
        numInf(ii) = sum(A(:,2).*matArray{ii}(:,n));
        numAsy(ii) = sum(A(:,3).*matArray{ii}(:,n));
    end
%     yplot(n+1) = 100*numInf(1)+250*numInf(2)+350*numInf(3)+330*numInf(4)+200*numInf(5)+100*numInf(6)+80*numInf(7)+50*numInf(8)+20*numInf(9)+10*numInf(10)+4*numInf(11)+2*numInf(12)+4*numInf(13);
% %     yplot(n+1) = sum(numInf);
    yplot(n+1) = 500*(numInf(4));
    newVal = 500*numInf(4);
    t(n+1)=t(n)+h;
    
end
% plot(t,yplot, '--b')
endInf = yplot(end);
        
    
