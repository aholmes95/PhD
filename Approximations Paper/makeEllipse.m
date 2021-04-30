figure;
% subplot(2,1,1)
% a=0.1;
% b=0.1;
% g=0.15;
N=500;
startPt = 50000;
drawPhasePortraitAll;
hold on
x = b1{1}(startPt:end)/N;
y = b2{1}(startPt:end)/N;
D = cov(x,y);
D = inv(D);
D11 = D(1,1);
D12 = D(1,2);
D21 = D(2,1);
D22 = D(2,2);
mu1 = mean(x);
mu2 =  mean(y);
plot(mu1,mu2,'k.','MarkerSize',25)
x1 = x;
y1 = y;
fp1 = mySS(1);
fp2 = mySS(2);
lightBlue = [91, 207, 244] / 255; 
scatter(fp1,fp2,25,lightBlue)
clear x y
f = @(x,y) D11*(x-mu1).^2 + (D12+D21)*(x-mu1).*(y-mu2) + D22*(y-mu2).^2 -5.991463;
% fimplicit(f,'LineWidth',2)
f = @(x,y) D11*(x-mu1).^2 + (D12+D21)*(x-mu1).*(y-mu2) + D22*(y-mu2).^2 - 9.21;
% fimplicit(f,'LineWidth',2)
a95 = D11*(x1-mu1).^2 + (D12+D21)*(x1-mu1).*(y1-mu2) + D22*(y1-mu2).^2 <= 5.991463;
a99 = D11*(x1-mu1).^2 + (D12+D21)*(x1-mu1).*(y1-mu2) + D22*(y1-mu2).^2 <= 9.21;
prop95 = sum(a95)/size(x1,1)
prop99 = sum(a99)/size(x1,1)


B = [evec1(1) evec2(1); evec1(2) evec2(2)];
invB = inv(B);
x = b1{1}(startPt:end);
y = b2{1}(startPt:end);
A = [x y];
A = [x-mean(x) y-mean(y)];
A = A';
C = invB*A/N;
hold on;
% subplot(3,1,2)
% scatter(C(1,:),C(2,:))
% mean(C(1,:))
% C = mnrnd(N,[mySS(1),mySS(2),1-mySS(1)-mySS(2)],10000);
% scatter(C(:,1)/N,C(:,2)/N,'r.')
% subplot(2,1,2)
% hh1 = binornd(N,mySSSHH(1),100000,1)/N;
% hh2 = binornd(N,mySSSHH(2),100000,1)/N;
% scatter(hh1,hh2)
% legend('','Gillespie Mean','','','','95%','')
% xlim([max(0,mySS(1)-0.12) min(1,mySS(1)+0.12)]);
% ylim([max(0,mySS(2)-0.12) min(1,mySS(2)+0.12)]);