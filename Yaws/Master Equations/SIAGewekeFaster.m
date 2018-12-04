function [endSIA] = SIAGewekeFaster()
% tmax = 20000;
clear all;
tic
tmax = 50000;
initInfNum = 160;
alpha = 0.164865777630342;
beta = 0.0516;
delta = 0.0513;
% rho = 0.0165;
rho = 0;
delta1 = delta;
delta2 = delta;
lambda = 0.185;
n=1;
t(1) = 0;
t1=0;
% numHouses = 500;
% eps=0.004;

% [houseMat, numPeop] = hhSizeDist(numHouses);

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

rateMat(:,1) = houseMat(:,5);
rateMat(:,2) = houseMat(:,6);
rateMat(:,3) = houseMat(:,7);
rateMat(:,4) = houseMat(:,8);
rateMat(:,5) = houseMat(:,9);

numPeop = sum(N)-1;

Iarray1(1) = initInfNum;
numI = initInfNum;

eps = alpha*initInfNum/numPeop;

rateMat(:,1) = (eps+beta*I./(N-1)).*S;
% rateMat(1:100,1) = eps*S(1:100);
rateMat(:,2) = delta1*I;
rateMat(:,3) = delta2*A;
rateMat(:,4) = rho*A;
rateMat(:,5) = lambda*I;
sumMat = sum(rateMat(:,1:5),2);
cumTotMat = cumsum(sumMat);

while t1<tmax
    if n>1
        eps = (alpha*(Iarray1(n))/numPeop);
        
%         if eventHouseNum>100
%             newHouseRate(1) = (eps+beta*I(eventHouseNum)/(N(eventHouseNum)-1))*S(eventHouseNum);
%         else    
%             newHouseRate(1) = eps*S(eventHouseNum);
%         end
        newHouseRate(1) = (eps+beta*I(eventHouseNum)/(N(eventHouseNum)-1))*S(eventHouseNum);
        newHouseRate(2) = delta1*I(eventHouseNum);
        newHouseRate(3) = delta2*A(eventHouseNum);
        newHouseRate(4) = rho*A(eventHouseNum);
        newHouseRate(5) = lambda*I(eventHouseNum);


        newSum = sum(newHouseRate);
        oldSum = sum(rateMat(eventHouseNum,:));
        sumMat(eventHouseNum) = sumMat(eventHouseNum) - oldSum + newSum;


%         if eventHouseNum>100
%             rateMat(eventHouseNum,1) = newHouseRate(1);
%         else    
%             rateMat(eventHouseNum,1) = newHouseRate(1);
%         end
        rateMat(eventHouseNum,1) = newHouseRate(1);
        rateMat(eventHouseNum,2) = newHouseRate(2);
        rateMat(eventHouseNum,3) = newHouseRate(3);
        rateMat(eventHouseNum,4) = newHouseRate(4);
        rateMat(eventHouseNum,5) = newHouseRate(5);
        
%         cumTotMat(eventHouseNum:end) = cumTotMat(eventHouseNum:end)-oldSum + newSum;
        cumTotMat = cumsum(sumMat);
    end
    randNums = rand(2,1);
    Rtotal = cumTotMat(end);
    dt=-log(randNums(1))/Rtotal;
    P=randNums(2)*Rtotal; 
    eventHouseNum = find(cumTotMat>=P,1);
    
    if eventHouseNum == 1
        remainP = P;
    else 
        remainP = P - cumTotMat(eventHouseNum-1);
    end
    withinCumTotMat = cumsum(rateMat(eventHouseNum, 1:5));
    eventNum = find(withinCumTotMat>=remainP,1);
    
    if eventNum == 1
        S(eventHouseNum) = S(eventHouseNum)-1;
        I(eventHouseNum) = I(eventHouseNum)+1;
        numI = numI+1;
    elseif eventNum == 2
        S(eventHouseNum) = S(eventHouseNum)+1;
        I(eventHouseNum) = I(eventHouseNum)-1;
        numI = numI-1;
    elseif eventNum == 3
        S(eventHouseNum) = S(eventHouseNum)+1;
        A(eventHouseNum) = A(eventHouseNum)-1;
    elseif eventNum == 4
        I(eventHouseNum) = I(eventHouseNum)+1;
        A(eventHouseNum) = A(eventHouseNum)-1;
        numI = numI+1;
    elseif eventNum == 5
        A(eventHouseNum) = A(eventHouseNum)+1;
        I(eventHouseNum) = I(eventHouseNum)-1;
        numI = numI-1;
    end
    
    n=n+1;
    
    if Rtotal==0
        endSIA = SIAGewekeFaster(alpha, beta, delta, rho, tmax, initInfNum);
        return;
    else
        t(n)=t(n-1)+dt;
        t1=t(n);
        Iarray1(n) = numI;
    end
end
Iarray1(1) = initInfNum;
Iarray1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Geweke Diagnostic stuff%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% numBins = 20;
% for i=1:numBins   
%     discardAmount = (i-1)/numBins*0.4;
%     newData = Iarray1(floor(discardAmount*(length(Iarray1)))+1:end);
%     A = newData(1:floor(0.1*(length(newData))));
%     B = newData(floor(0.5*(length(newData))+1):end);
%     meanA = mean(A);
%     meanB = mean(B);
%     varA = var(A);
%     varB = var(B);
%     z(i) = (mean(A)-mean(B))/sqrt(varA+varB);
%     if abs(z(i)) < 2
%         flag=1;
%         break;
%     else 
%         flag=0;
%     end
% end
% 
% if flag==0
%     tMaxNew = tmax*1.1;
%     endSIA = SIAGewekeFaster(alpha, beta, delta, rho, tMaxNew, initInfNum);
%     return;
% end
    

% endSIA = meanB;
endSIA = mean(Iarray1(round(length(Iarray1/2)):end));
toc