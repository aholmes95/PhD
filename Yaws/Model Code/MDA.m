function [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, MDATreatmentCost, totFinCost, totEconCost] = MDA(eventTime, dt, numHouses, massCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, MDATreatmentCost, totFinCost, totEconCost, costParams)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    finCost = costParams(5);
    econCost = costParams(6);
    
    if t(n) < eventTime && t(n) + dt >= eventTime
        for j=1:numHouses
           treatedRand = rand;
           if treatedRand < massCoverage
               numInHouse = N(j);
               numInfInHouse = I(j);
               numAsInHouse = A(j);
               msize = numel(ageDist);
               agesOfPatients = ageDist(randperm(msize, N(j)));
               for jj=1:N(j)
                   totFinCost = totFinCost + finCost;
                   totEconCost = totEconCost + econCost;
                   MDATreatmentCost = MDATreatmentCost + cost(agesOfPatients(jj), costParams);
               end
         %     newI = interp1(t, Iarray, tq, 'previous');
%     newA = interp1(t, Iarray3, tq, 'previous');
      effInfRand = rand(numInfInHouse, 1);
               effAsRand = rand(numAsInHouse,1);
               numInfTreated = sum(effInfRand < efficacy);
               numAsTreated = sum(effAsRand < efficacy);

               I(j) = I(j) - numInfTreated;
               
               A(j) = A(j) - numAsTreated;
               
               S(j) = S(j) + numInfTreated + numAsTreated;
           end
        end
        n=n+1;
        t(n)=eventTime;
        numI=sum(I);
        numA=sum(A);
        Iarray(n)=numI;
        Iarray2(n)=numI+numA;
        Iarray3(n) = numA;
        flag = 1;
    else
        flag = 0;
    end
end

