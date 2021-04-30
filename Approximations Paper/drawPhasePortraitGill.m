function [b1 b2] = drawPhasePortraitGill(startPt,myDir,N)
addpath(myDir);
myPath1 = strcat(myDir,'/*state1.csv');
myPath2 = strcat(myDir,'/*state2.csv');
myPath3 = strcat(myDir,'/*state3.csv');
d1 = dir(myPath1);
d2 = dir(myPath2);
d3 = dir(myPath3);
numFiles = size(d1,1);
b1 = {};
b2 = {};
b3 = {};
hold on;
% x = 2600:2800;
% y = 3000-x;
% plot(x,y)
for i=1
    b1{i} = csvread(d1(i).name);
    b2{i} = csvread(d2(i).name);
    b3{i} = csvread(d3(i).name);
    N = b1{1}(end)+b2{1}(end)+b3{1}(end);
%     scatter(b1{i}(200000:end),b2{i}(200000:end));
    
%     for j=1:mySize
%         scatter(b1{i}(2000000+20*(j-1):20:2000000+20*j),b2{i}(200000+20*(j-1):20:200000+20*j),'b.');
%         pause(1e-10);
%     end    
    scatter(b1{i}(startPt:end)/N,b2{i}(startPt:end)/N,'g.');
end
mean1 = mean(b1{1}(startPt:end)/N)
mean2 = mean(b2{1}(startPt:end)/N)