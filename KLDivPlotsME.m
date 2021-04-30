clear all;
numParams = 75;
% KLMat = zeros(numParams,3);
% alphaVec = linspace(0.15,0.8,numParams);
% betaVec = linspace(0.1,0.8,numParams);
gammaVec = linspace(0.05,0.7,numParams);
% NVec = linspace(100,500,numParams);
% beta = linspace(0.1,1,75);
alpha = 0.5;
beta = 0.5;
% gamma = 0.15;
numHouses = 250;
MEDists = zeros(numParams,3);
for i=1:numParams %Calculate steady state distributions in parallel
%     alpha = alphaVec(i);
%     beta = betaVec(i);
    gamma = gammaVec(i);
    G = findSSDistME(alpha,beta,gamma,numHouses);
    MEDists(i,:) = G;
end


for j=1:numParams
    j
%     alpha = alphaVec(j);
%     beta = betaVec(j);
    gamma = gammaVec(j);
    myParams = [alpha; beta; gamma; numHouses];
    N = numHouses;
    command1='wolframscript -file test1.wls';
    csvwrite("myParams.csv",myParams);
    status1 = system(command1);
    mySS=csvread("mySS.csv");
    mySSSHH=csvread("mySSSHH.csv");
    mySSMC=csvread("mySSMC.csv");
    FP = [mySS(1) mySS(2) 1-mySS(1)-mySS(2)];
    SHH = [1-mySSSHH(1)-mySSSHH(2) mySSSHH(1) mySSSHH(2)];
%     MC = [mySSMC(1) mySSMC(2) N-mySSMC(1)-mySSMC(2)]/N;
%     G = findSSDistME(alpha,beta,gamma,N);
    G = MEDists(j,:);
    myDivFP = CalKLDiv(FP, G);
    myDivSHH = CalKLDiv(SHH, G);
%     myDivMC = CalKLDiv(MC,G);

    KLMatFP(j) = myDivFP;
    KLMatSHH(j) = myDivSHH;
%     KLMatMC(j) = myDivMC;
end

% a = KLMat(KLMat(:,2)==0.1 & KLMat(:,3)==0.15,:);
% plota = KLMat(:,1);
% fpkldiva = KLMatFP(:,1);
% shhkldiva = KLMatSHH(:,2);
% mckldiva = KLMatMC(:,3);

fpkldivb = KLMatFP;
shhkldivb = KLMatSHH;
% mckldivb = KLMatMC;

figure;
% heatmap(alphaVec',betaVec', KLMatFP);
% figure;
% heatmap(alphaVec',betaVec',KLMatSHH);
% figure;
% heatmap(alphaVec',betaVec',KLMatMC);
hold on;
plot(betaVec, fpkldivb,'r','LineWidth',2);
plot(betaVec, shhkldivb,'b','LineWidth',2);
% plot(alphaVec, mckldiva,'g','LineWidth',2);
% legend('Fokker Planck','Single Household','Moment Closure');
% scatter(plota, fpkldiva,'rx');
% scatter(plota, shhkldiva,'bx');
xlabel('Gamma');
ylabel('KL-Divergence');

% b = KLMat(KLMat(:,1)==0.1 & KLMat(:,3)==0.1,:);
% plotb = b(:,2);
% fpkldivb = b(:,5);
% shhkldivb = b(:,6);
% mckldivb = b(:,7);
% figure;
% hold on;
% plot(plotb, fpkldivb,'r','LineWidth',2);
% plot(plotb, shhkldivb,'b','LineWidth',2);
% plot(plotb, mckldivb,'g','LineWidth',2);
% legend('Fokker Planck','Single Household','Moment Closure');
% % scatter(plotb, fpkldivb,'rx');
% % scatter(plotb, shhkldivb,'bx');
% xlabel('Beta');
% ylabel('KL-Divergence');
% 
% g = KLMat(KLMat(:,1)==0.5 & KLMat(:,2)==0.5,:);
% plotg = g(:,3);
% fpkldivg = g(:,5);
% shhkldivg = g(:,6);
% mckldivg = g(:,7);
% figure;
% hold on;
% plot(plotg, fpkldivg,'r','LineWidth',2);
% plot(plotg, shhkldivg,'b','LineWidth',2);
% plot(plotg, mckldivg,'g','LineWidth',2);
% legend('Fokker Planck','Single Household','Moment Closure');
% % scatter(plotg, fpkldivg,'rx');
% % scatter(plotg, shhkldivg,'bx');
% xlabel('Gamma');
% ylabel('KL-Divergence');
