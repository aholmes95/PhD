clear;
clc;
tmax = 700;
numSims2 = 1600;
tq = 0:0.0005:tmax;
sumSimInf = zeros(1,length(tq));
delta1 = 0.0513;
delta2 = 0.0513;
% delta2 = 0;
beta = 0.0516;
rho = 0.0165;
% rho = 0;
lambda = 0.185;
% lambda = 0
initInfNum = 160;
% alpha = 0.1662;
% alpha = 0.164865777630342;
eps=0.004;
parfor ii=1:numSims2
    ii
    n=1;
    t = [];
    t(1) = 0;
    t1=0;
    I1 = initInfNum;

    houseMat = createHouseMat;
    sizeHouseMat = size(houseMat);
    numHouses = sizeHouseMat(1);
    initInf= randperm(numHouses,initInfNum);
    houseMat(initInf,2) = houseMat(initInf,2) - 1;
    houseMat(initInf,3) = houseMat(initInf,3) + 1;

    N = houseMat(:,1);
    S = houseMat(:,2);
    I = houseMat(:,3);
    A = houseMat(:,4);
    
    numPeop = sum(N);
    
    Sarray = [numPeop - initInfNum];
    Iarray1 = [initInfNum];
    Iarray2 = [0];
    Iarray3 = [initInfNum];
    

    while t1<tmax
%         if n==1
%             eps = alpha*initInfNum/numPeop;
%         else
%             eps = (alpha*(Iarray1(n))/numPeop);
%         end
%         eps=0.004;
        houseMat(:,5) = (eps+beta*I./(N-1)).*S;
        houseMat(N==1,5) = eps*S(N==1);
        houseMat(:,6) = delta1*I;
        houseMat(:,7) = delta2*A;
        houseMat(:,8) = rho*A;
        houseMat(:,9) = lambda*I;
        houseMat(:,10) = sum(houseMat(:,5:9),2);
        rand1=rand; rand2 = rand;
        Rtotal = sum(houseMat(:,10));
        dt=-log(rand1)/Rtotal;
        P=rand2*Rtotal; 
        cumTotMat = cumsum(houseMat(:,10));
        b = find(cumTotMat>=P);
        eventHouseNum = b(1);
        if eventHouseNum == 1
            withinHouseRate = cumTotMat(eventHouseNum);
            remainP = P;
        else 
            withinHouseRate = cumTotMat(eventHouseNum) - cumTotMat(eventHouseNum-1);
            remainP = P - cumTotMat(eventHouseNum-1);
        end
        withinCumTotMat = cumsum(houseMat(eventHouseNum, 5:9));
        c = find(withinCumTotMat>=remainP);
        eventNum = c(1);
        if eventNum == 1
            S(eventHouseNum) = S(eventHouseNum)-1;
            I(eventHouseNum) = I(eventHouseNum)+1;
        elseif eventNum == 2
            S(eventHouseNum) = S(eventHouseNum)+1;
            I(eventHouseNum) = I(eventHouseNum)-1;
        elseif eventNum == 3
            S(eventHouseNum) = S(eventHouseNum)+1;
            A(eventHouseNum) = A(eventHouseNum)-1;
        elseif eventNum == 4
            I(eventHouseNum) = I(eventHouseNum)+1;
            A(eventHouseNum) = A(eventHouseNum)-1;
        elseif eventNum == 5
            A(eventHouseNum) = A(eventHouseNum)+1;
            I(eventHouseNum) = I(eventHouseNum)-1;
        end
        n=n+1;
        if Rtotal==0
            t(n) = tmax
            Iarray1(n) = 0;
            Iarray2(n) = 0;
            Iarray3(n) = 0;
            t1=tmax;
        else
            t(n)=t(n-1)+dt;
            t1=t(n);
            numI = sum(I);
            numA = sum(A);
            Iarray1(n) = numI;
            Iarray2(n) = numA;
            Iarray3(n) = numI + numA;
        end
        I1 = Iarray3(n);
    end
    
    a1=0;
    a2=0;
    a3=0;
    a4=0;
    a5=0;
    a6=0;
    a7=0;
    a8=0;
    a9=0;
    a10=0;
    a11=0;
    a12=0;
    a13=0;
    a14=0;
    a15=0;

    for i=1:1500
        if S(i) == 0 && I(i) == 0
            a1 = a1+1;
        elseif S(i) == 0 && I(i) == 1
            a2 = a2+1;
        elseif S(i) == 0 && I(i) == 2
            a3 = a3+1;
        elseif S(i) == 0 && I(i) == 3
            a4 = a4+1;
        elseif S(i) == -4 && I(i) == 4
            a5 = a5+1;
        elseif S(i) == 1 && I(i) == 0
            a6 = a6+1;
        elseif S(i) == 1 && I(i) == 1
            a7 = a7+1;
        elseif S(i) == 1 && I(i) == 2
            a8 = a8+1;
        elseif S(i) == 1 && I(i) == 3
            a9 = a9+1;
        elseif S(i) == 2 && I(i) == 0
            a10 = a10+1;
        elseif S(i) == 2 && I(i) == 1
            a11 = a11+1;
        elseif S(i) == 2 && I(i) == 2
            a12 = a12+-1;
        elseif S(i) == 3 && I(i) == 0
            a13 = a13+1;
        elseif S(i) == 3 && I(i) == 1
            a14 = a14+1;
        elseif S(i) == 4 && I(i) == 0
            a15 = a15+1;
        end
    end
    dist(:,ii) = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14; a15]/1500;
%     hold on;

    newI1 = interp1(t, Iarray1, tq, 'previous');
%     newI2 = interp1(t, Iarray2, tq, 'previous');
%     newI3 = interp1(t, Iarray3, tq, 'previous');
    IMat1(ii,:) = newI1;
%     IMat2(ii,:) = newI2;
%     IMat3(ii,:) = newI3;
    sumSimInf = sumSimInf + newI1;
    
end
% GillDist = mean(dist,2);
% 
Iplot1 = mean(IMat1,1);
% Iplot2 = mean(IMat2,1);
% Iplot3 = mean(IMat3,1);
% Iplot1 = sumSimInf/numSims2;
% Iplot1(1) = initInfNum;
% stairs(tq, Iplot1, '--k');
% stairs(t, Iarray1)
% hold on;
% stairs(t, Iarray2)
% hold on;
stairs(tq, Iplot1, '--b')
% hold on;


    
