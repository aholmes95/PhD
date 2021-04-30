addpath('Distributions');
myFiles = dir('Distributions/*.txt');
numFiles = size(myFiles,1);
KLMat = zeros(53,7);
for i=1:numFiles
    fileID = fopen(myFiles(i).name);
    bI{i} = textscan(fileID, '%s %s %f');
    fclose(fileID);
end
for i=1:numFiles
    alpha = bI{i}{3}(1);
    beta = bI{i}{3}(2);
    gamma = bI{i}{3}(3);
    numHouses = bI{i}{3}(4);
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
    MC = [mySSMC(1) mySSMC(2) N-mySSMC(1)-mySSMC(2)]/N
    numState1 = bI{i}{3}(5);
    numState2 = bI{i}{3}(6);
    numState3 = bI{i}{3}(7);
    sumState = numState1+numState2+numState3;
    G = [numState1 numState2 numState3]/sumState;
    KLMat(i,1) = alpha;
    KLMat(i,2) = beta;
    KLMat(i,3) = gamma;
    KLMat(i,4) = numHouses;
    myDivFP = CalKLDiv(FP, G);
    myDivSHH = CalKLDiv(SHH, G);
    myDivMC = CalKLDiv(MC,G);
    KLMat(i,5) = myDivFP;
    KLMat(i,6) = myDivSHH;
    KLMat(i,7) = myDivMC;
end

a = KLMat(KLMat(:,2)==0.1 & KLMat(:,3)==0.15,:);
plota = a(:,1);
fpkldiva = a(:,5);
shhkldiva = a(:,6);
mckldiva = a(:,7);
figure;
hold on;
plot(plota, fpkldiva,'r','LineWidth',2);
plot(plota, shhkldiva,'b','LineWidth',2);
plot(plota, mckldiva,'g','LineWidth',2);
legend('Fokker Planck','Single Household','Moment Closure');
% scatter(plota, fpkldiva,'rx');
% scatter(plota, shhkldiva,'bx');
xlabel('Alpha');
ylabel('KL-Divergence');

b = KLMat(KLMat(:,1)==0.1 & KLMat(:,3)==0.1,:);
plotb = b(:,2);
fpkldivb = b(:,5);
shhkldivb = b(:,6);
mckldivb = b(:,7);
figure;
hold on;
plot(plotb, fpkldivb,'r','LineWidth',2);
plot(plotb, shhkldivb,'b','LineWidth',2);
plot(plotb, mckldivb,'g','LineWidth',2);
legend('Fokker Planck','Single Household','Moment Closure');
% scatter(plotb, fpkldivb,'rx');
% scatter(plotb, shhkldivb,'bx');
xlabel('Beta');
ylabel('KL-Divergence');

g = KLMat(KLMat(:,1)==0.5 & KLMat(:,2)==0.5,:);
plotg = g(:,3);
fpkldivg = g(:,5);
shhkldivg = g(:,6);
mckldivg = g(:,7);
figure;
hold on;
plot(plotg, fpkldivg,'r','LineWidth',2);
plot(plotg, shhkldivg,'b','LineWidth',2);
% plot(plotg, mckldivg,'g','LineWidth',2);
legend('Fokker Planck','Single Household','Moment Closure');
% scatter(plotg, fpkldivg,'rx');
% scatter(plotg, shhkldivg,'bx');
xlabel('Gamma');
ylabel('KL-Divergence');
