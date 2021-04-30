% clear all;
numParams = 50;
% KLMat = zeros(numParams,3);
% NVec = linspace(100,500,numParams);
% beta = linspace(0.1,1,75);
% alphaVec = linspace(0.4,0.8,numParams);
% alpha1 = alphaVec(1:20);
% alpha2 = alphaVec(21:40);
% alpha = alphaVec(41:50);
% betaVec = linspace(0.25,0.8,numParams);
% beta1 = betaVec(1:20);
% beta2 = betaVec(21:35);
% beta3 = betaVec(36:50);
alpha = 0.4;
% beta = 0.5;
gammaVec = linspace(0.1,0.5,numParams);
betaVec = linspace(0.25,0.8,numParams);

betaVec = betaVec(36:50);

% alpha = 0.4;
% beta = 0.25;
% gamma = 0.15;
numHouses = 250;
% MEDists = zeros(numParams,3);
% FPDists = zeros(numParams,3);
% SHHDists = zeros(numParams,3);
% MCDists = zeros(numParams,3);
% parfor i=1:numParams %Calculate steady state distributions in parallel
%     gamma = gammaVec(i)
%     G = findSSDistME(alpha,beta,gamma,numHouses);
%     MEDists(i,:) = G;
% end


% parfor j=1:length(beta3)
% parfor i=1:numParams
% %         beta = beta3(j)
% %         disp(beta)
%     gamma = gammaVec(i)
% %         disp(gamma)
% %         gamma = gammaVec(j);
%     N = numHouses;
%     command1=['wolframscript -file test1Par.wls', ' ', num2str(alpha), ' ', num2str(beta), ' ', num2str(gamma), ' ', num2str(N), ' ', num2str(i)];
%     status1 = system(command1);
%     SSfilename = "mySS"+num2str(i)+".csv";
%     mySS=csvread(SSfilename);
%     SSSHHfilename = "mySSSHH"+num2str(i)+".csv";
%     mySSSHH=csvread(SSSHHfilename);
%     SSMCfilename = "mySSMC"+num2str(i)+".csv";
%     mySSMC=csvread(SSMCfilename);
%     FP = [mySS(1) mySS(2) 1-mySS(1)-mySS(2)];
%     SHH = [1-mySSSHH(1)-mySSSHH(2) mySSSHH(1) mySSSHH(2)];
%     MC = [mySSMC(1) mySSMC(2) N-mySSMC(1)-mySSMC(2)]/N;
%     G = findSSDistME(alpha,beta,gamma,N);
% %         G = MEDists(j,:);
%     myDivFP = CalKLDiv(FP, G);
%     myDivSHH = CalKLDiv(SHH, G);
%     myDivMC = CalKLDiv(MC,G);
% 
%     KLMatFP(i) = myDivFP;
%     KLMatSHH(i) = myDivSHH;
%     KLMatMC(i) = myDivMC;
% end

parfor i=1:numParams
    i
    gamma = gammaVec(i);
    for j=1:15
        beta = betaVec(j);
        N = numHouses;
        command1=['wolframscript -file test1Par.wls', ' ', num2str(alpha), ' ', num2str(beta), ' ', num2str(gamma), ' ', num2str(N), ' ', num2str(i)];
        status1 = system(command1);
        SSfilename = "mySS"+num2str(i)+".csv";
        mySS=csvread(SSfilename);
        SSSHHfilename = "mySSSHH"+num2str(i)+".csv";
        mySSSHH=csvread(SSSHHfilename);
        SSMCfilename = "mySSMC"+num2str(i)+".csv";
        mySSMC=csvread(SSMCfilename);
        FP = [mySS(1) mySS(2) 1-mySS(1)-mySS(2)];
        SHH = [1-mySSSHH(1)-mySSSHH(2) mySSSHH(1) mySSSHH(2)];
        MC = [mySSMC(1) mySSMC(2) N-mySSMC(1)-mySSMC(2)]/N;
        G = findSSDistME(alpha,beta,gamma,N);
    %         G = MEDists(j,:);
        myDivFP = CalKLDiv(FP, G);
        myDivSHH = CalKLDiv(SHH, G);
        myDivMC = CalKLDiv(MC,G);

        KLMatFP(i,j) = myDivFP;
        KLMatSHH(i,j) = myDivSHH;
        KLMatMC(i,j) = myDivMC;
    end
end
% end

% a = KLMat(KLMat(:,2)==0.1 & KLMat(:,3)==0.15,:);
% plota = KLMat(:,1);
% fpkldiva = KLMatFP(:,1);
% shhkldiva = KLMatSHH(:,2);
% mckldiva = KLMatMC(:,3);

% fpkldiva = KLMatFP;
% shhkldiva = KLMatSHH;
% mckldiva = KLMatMC;

% figure;
% heatmap(alphaVec',betaVec', KLMatFP);
% figure;
% heatmap(alphaVec',betaVec',KLMatSHH);
% figure;
% heatmap(alphaVec',betaVec',KLMatMC);
% hold on;
% plot(gammaVec, fpkldiva,'r','LineWidth',2);
% plot(gammaVec, shhkldiva,'b','LineWidth',2);
% plot(gammaVec, mckldiva,'g','LineWidth',2);
% legend('Fokker Planck','Single Household','Moment Closure');
% % scatter(plota, fpkldiva,'rx');
% % scatter(plota, shhkldiva,'bx');
% xlabel('Gamma');
% ylabel('KL-Divergence');

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
% ratKLMat = KLMatFP./KLMatSHH;
% custommap = createCustomCMap(ratKLMat);
% h=heatmap(gammaVec,alphaVec,ratKLMat);
% h.Colormap = custommap;
