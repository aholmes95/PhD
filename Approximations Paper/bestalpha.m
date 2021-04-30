function output = bestalpha(e,b,g,N)

    epsSS = [g^2/((e+g)^2+b*e), 2*g*e/((e+g)^2+b*e), (e^2+b*e)/((e+g)^2+b*e)];
    command1 = ['wolframscript -file findAlphaWolfram.wls', ' ', num2str(e), ' ', num2str(b), ' ', num2str(g), ' ', num2str(N)];
    status1 = system(command1);
    fpalphafilename = "myFPalpha.csv";
    fpalpha = csvread(fpalphafilename);
    [minalpha, kldivmin] = solveforalpha(e,b,g,N);
    shhalpha = ((e + g)^2 + b*e)/(g + e + b);
    shhdist = findSSDistME(shhalpha,b,g,N);
    shhkldiv = CalKLDiv(epsSS,shhdist);
    fpdist = findSSDistME(fpalpha,b,g,N);
    fpkldiv = CalKLDiv(epsSS,fpdist);
    
    output = [kldivmin, shhkldiv, fpkldiv];
end

