function [a] = createHouseMat
%This file creates a matrix that stores all the information about each of
%the households.
%Column1 = Number of people in each house
%Column2 = Number of susceptibles in each house
%Column3 = Number of infectious in each house
%Column4 = Number of asymptomatics  in each house.
%Column5 = S->S-1, I->I+1 rate
%Column6 = S->S+1, I->I-1 rate
%Column7 = A->A-1, S->S+1 rate
%Column8 = A->A-1, I->I+1 rate
%Column9 = I->I-1, A->A+1 rate
%Column10 = Total rate for household


houseMat = zeros(750,9);
houseMat(1:100,1:2) = 1;
houseMat(101:300,1:2) = 2;
houseMat(301:700,1:2) = 3;
houseMat(701:1030,1:2) = 4;
houseMat(1031:1230,1:2) = 5;
houseMat(1231:1330,1:2) = 6;
houseMat(1331:1410,1:2) = 7;
houseMat(1411:1460,1:2) = 8;
houseMat(1461:1480,1:2) = 9;
houseMat(1481:1490,1:2) = 10;
houseMat(1491:1494,1:2) = 11;
houseMat(1495:1496,1:2) = 12;
houseMat(1497:1500,1:2) = 13;

a = houseMat;
end

