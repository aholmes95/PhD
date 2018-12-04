function [S, I, A, t, n, Iarray, Iarray2, Iarray3, flag, TTTTreatmentCost, totFixedCost] = TTT(eventTime, dt, numHouses, contactCoverage, efficacy, S, I, A, N, t, Iarray, Iarray2, Iarray3, n, ageDist, TTTTreatmentCost, totFixedCost, costParams)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
        
    fixedCost = costParams(7);
    
    if t(n) < eventTime && t(n) + dt >= eventTime
        for j=1:numHouses
            if I(j) > 0
                infContCovRand = rand(I(j),1);
                numInfCov = sum(infContCovRand < contactCoverage);
                msize = numel(ageDist);
                agesOfPatients = ageDist(randperm(msize, numInfCov));
                for jj=1:numInfCov
                    totFixedCost = totFixedCost + fixedCost;
                    TTTTreatmentCost = TTTTreatmentCost + cost(agesOfPatients(jj), costParams);
                end
                effInfRand = rand(numInfCov,1);
                numInfTreated = sum(effInfRand < efficacy);
                I(j) = I(j) - numInfTreated;
                S(j) = S(j) + numInfTreated;

                if A(j) > 0
                    asyContCovRand = rand(A(j),1);
                    numAsyCov = sum(asyContCovRand < contactCoverage);
                    agesOfPatients = ageDist(randperm(msize, numAsyCov)); 
                    for jj=1:numAsyCov
                        totFixedCost = totFixedCost + fixedCost;
                        TTTTreatmentCost = TTTTreatmentCost + cost(agesOfPatients(jj), costParams);
                    end
                    effAsRand = rand(numAsyCov,1);
                    numAsTreated = sum(effAsRand < efficacy);
                    A(j) = A(j) - numAsTreated;
                    S(j) = S(j) + numAsTreated;
                end
                
                if S(j) > 0
                    susContCovRand = rand(S(j),1);
                    numSusCov = sum(susContCovRand < contactCoverage);
                    agesOfPatients = ageDist(randperm(msize, numSusCov)); 
                    for jj=1:numSusCov
                        totFixedCost = totFixedCost + fixedCost;
                        TTTTreatmentCost = TTTTreatmentCost + cost(agesOfPatients(jj), costParams);
                    end
                    
                end
            end
        end
        n=n+1;
        t(n)=eventTime;
        numI=sum(I);
        numA=sum(A);
        Iarray(n)=numI;
        Iarray2(n)=numI+numA;
        Iarray3(n)=numA;
        flag=1;
    else
        flag=0;
    end

end