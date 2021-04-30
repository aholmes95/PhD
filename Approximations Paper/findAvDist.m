function gillDist = findAvDist()
myDir = 'Distributions2';
addpath(myDir);
myPath1 = strcat(myDir,'/*.csv');
myPath2 = strcat(myDir,'/*.csv');
myPath3 = strcat(myDir,'/*.csv');
d1 = dir(myPath1);
d2 = dir(myPath2);
d3 = dir(myPath3);
numFiles = length(d1);
state1 = 0;
state2 = 0;
state3 = 0;
for i=1:numFiles
    b1 = csvread(d1(i).name);
    state1 = state1 + b1(1)/sum(b1);
    state2 = state2 + b1(2)/sum(b1);
    state3 = state3 + b1(3)/sum(b1);
end
state1 = state1/numFiles;
state2 = state2/numFiles;
state3 = state3/numFiles;
gillDist = [state1 state2 state3];