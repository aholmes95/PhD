folder = 'State Files';
myParams = [a; b; g; N];
command1='wolframscript -file test1.wls';
command2 = 'g++ -O3 -std=c++11 -o SIS2HousePhasePortrait SIS2HousePhasePortrait.cpp'
command3 = ['./SIS2HousePhasePortrait ', num2str(a), ' ', num2str(b), ' ', num2str(g), ' ', num2str(N)]
status2 = system(command2); 
status3 = system(command3);
csvwrite("myParams.csv",myParams);
status1 = system(command1);
evec1=csvread("myEvec1.csv");
evec1 = evec1/norm(evec1);
evec2=csvread("myEvec2.csv");
evec2 = evec2/norm(evec2);
mySS=csvread("mySS.csv");

mySSSHH=csvread("mySSSHH.csv");
[b1 b2] = drawPhasePortraitGill(startPt,folder,N)
% C = mnrnd(N,[mySS(1),mySS(2),1-mySS(1)-mySS(2)],10000);
% scatter(C(:,1)/N,C(:,2)/N,'r.')
drawPhasePortraitFP(a,b,g,N,mySS,evec1,evec2);
xlim([max(0,mySS(1)-0.12) min(1,mySS(1)+0.12)]);
ylim([max(0,mySS(2)-0.12) min(1,mySS(2)+0.12)]);
hold on;
xlabel('Household State (2,0)');
ylabel('Household State (1,1)');
myTitle = ['a = ' num2str(a) ', b = ' num2str(b) ', g = ' num2str(g) ', N = ' num2str(N)];
title(myTitle)
% saveas(gcf,['Phase Portraits/' myTitle '.png']);
