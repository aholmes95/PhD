clear;
clc;
tmax = 1000;
numSims2 = 500;
tq = 0:0.05:tmax;
initInfNum = 160;
numI = initInfNum;
ageMat = csvread('ageData.csv');
ageDist = ageMat(:,2);
MDACostMat = zeros(numSims2, 1);
TTTCostMat = zeros(numSims2, 1);
timeToErad = zeros(numSims2, 1);

%%%%%%%%%%%%%%%%%
%Cost Parameters%
%%%%%%%%%%%%%%%%%

minAge = 5;
maxAge = 14;
minDoseCost = 0.27; %in US Dollars
maxDoseCost = 1.1;  %in US Dollars
finCost = 3.44; %in US Dollars (cost per person)
econCost = 4.76; % in US Dollars (cost per person)
fixedCost = 21400; 
costParams = [minAge, maxAge, minDoseCost, maxDoseCost, finCost, econCost, fixedCost];

%%%%%%%%%%%%%%%%%%%%%%
%Treatment parameters%
%%%%%%%%%%%%%%%%%%%%%%

contactCoverage = 1; %Probability that each contact will accept treatment. Per person chance, NOT per household.
massCoverage = 0.9; %Probability that a house receives treatment during MDA. Assumed the entire house would accept.
efficacy = 0.8; %This is the probability that each dose of the drug will work.
numExt = 0; %Number of times the disease is wiped out.

parfor ii=1:numSims2
    ii
    numI = initInfNum;
    numA = 0;
    delta1 = 0.0513;
    delta2 = 0.0513;
    beta = 0.0516;
    rho = 0.0165;
    lambda = 0.185;
    alpha = 0.1680;
    
    
    n=1;
    t1=0;
    t=[];
    t(1)=0;
    Iarray=[];
    Iarray2=[];
    Iarray3=[];
    Iarray(1) = initInfNum;
    Iarray2(1) = initInfNum;
    Iarray3(1) = 0;
    MDATreatmentCost = 0;
    TTTTreatmentCost = 0;
    
    totFinCost = 0;
    totEconCost = 0;
    totFixedCost = 0;

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

    while t1<tmax
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Calculate rates and time until mext event%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        flag = 0;
        numInf = sum(I);
        eps = alpha*numInf/5930;
        houseMat(:,5) = (eps+beta*I./(N-1)).*S;
        houseMat(N==1,5) = eps*S(N==1);
        houseMat(:,6) = delta1*I;
        houseMat(:,7) = delta2*A;
        houseMat(:,8) = rho*A;
        houseMat(:,9) = lambda*I;
        houseMat(:,10) = sum(houseMat(:,5:9),2);
        rand1=rand; rand2 = rand;
        Rtotal = sum(houseMat(:,10));
        
        if Rtotal == 0
            t(n+1) = t(n)+0.0001;
            t(n+2) = tmax;
            Iarray(n+1) = 0;
            Iarray(n+2) = 0;
            Iarray2(n+1) = 0;
            Iarray2(n+2) = 0;
            Iarray3(n+1) = 0;
            Iarray3(n+2) = 0;
            numExt = numExt + 1;
            timeToErad(ii) = t(n);
            break;
        end
            
        dt=-log(rand1)/Rtotal;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%    MDA    %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, MDATreatmentCost, totFinCost, totEconCost] = MDA(300, dt, numHouses, massCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, MDATreatmentCost, totFinCost, totEconCost, costParams);
%         if flag==1
%            t1=t(n);
%            continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, MDATreatmentCost, totFinCost, totEconCost] = MDA(312, dt, numHouses, massCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, MDATreatmentCost, totFinCost, totEconCost, costParams);
%         if flag==1
%            t1=t(n);
%            continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, MDATreatmentCost, totFinCost, totEconCost] = MDA(324, dt, numHouses, massCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, MDATreatmentCost, totFinCost, totEconCost, costParams);
%         if flag==1
%            t1=t(n);
%            continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, MDATreatmentCost, totFinCost, totEconCost] = MDA(336, dt, numHouses, massCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, MDATreatmentCost, totFinCost, totEconCost, costParams);
%         if flag==1
%            t1=t(n);
%            continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, MDATreatmentCost, totFinCost, totEconCost] = MDA(348, dt, numHouses, massCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, MDATreatmentCost, totFinCost, totEconCost, costParams);
%         if flag==1
%            t1=t(n);
%            continue;
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%    TTT    %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(300, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(312, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(324, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(336, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(348, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
        
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(360, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(372, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(384, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(396, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(408, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
        
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(420, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(432, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(444, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(456, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
        [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(468, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
        if flag==1
            t1=t(n);
            continue;
        end
% 
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(550, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
% 
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(600, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(650, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
% 
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(700, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
% 
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(750, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
% 
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(800, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
% 
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(850, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(900, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(950, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(1000, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(1050, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(1100, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(1150, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end
%         
%         [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(1200, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams);
%         if flag==1
%             t1=t(n);
%             continue;
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Determine which house the event takes place in%
        %%%%and which event it is that takes place.%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Update each array according to which event took place%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        t(n)=t(n-1)+dt;
        t1=t(n);
        numI=sum(I);
        numA=sum(A);
        Iarray(n)=numI;
        Iarray2(n) = numI+numA;
        Iarray3(n) = numA;
        
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Interpolate by taking previous value at fixed time steps%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     newI = interp1(t, Iarray, tq, 'previous');
%     newA = interp1(t, Iarray3, tq, 'previous');
%     IMat(ii,:) = newI;
%     AMat(ii,:) = newA;
    MDACost = MDATreatmentCost% + totEconCost + totFinCost;
    TTTCost = TTTTreatmentCost% + totFixedCost;
    MDACostMat(ii) = MDACost;
    TTTCostMat(ii) = TTTCost;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove simulations that resulted in the disease dying out
%%%%%%%%%        due to stochastcicity        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%zeroInd = [];
%for i=1:numSims2
%    if IMat(i,end) == 0
%        zeroInd = [zeroInd, i];
%    end
%end
%IMat2 = IMat(:,:);
%IMat2(zeroInd,:) = [];
%treatmentCostMat2(:) = treatmentCostMat;
%treatmentCostMat2(~zeroInd) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate mean of the simulations and plot results%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iplot = mean(IMat,1);
% Aplot = mean(AMat,1);
% Iplot(1) = initInfNum;
% Aplot(1) = 0;
% stairs(tq, Iplot+Aplot);
% simResults1 = {IMat, MDACostMat, TTTCostMat, timeToErad, numExt};
% save('simResults1','simResults1','-v7.3'); 
% save('Iplot4.mat','Iplot');
% save('Aplot4.mat','Aplot');

    
