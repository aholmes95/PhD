function [c] = cost(a, costParams)
    minAge = costParams(1); maxAge = costParams(2); minDoseCost = costParams(3); maxDoseCost = costParams(4);
    if a <= minAge
        c = minDoseCost;
    elseif a >= maxAge
        c = maxDoseCost;
    else
        c= ((maxDoseCost-minDoseCost)/(maxAge-minAge))*a + minDoseCost - (maxDoseCost-minDoseCost)/(maxAge-minAge)*minAge;
    end
end
