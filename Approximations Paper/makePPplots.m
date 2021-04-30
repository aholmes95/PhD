n=30;
R0thresh = 1.1;
X = lhsdesign(n,3,'criterion','correlation','iterations',10)
X2=X((X(:,1)+sqrt(X(:,1).^2+4*X(:,1).*X(:,2)))./(2*X(:,3))>R0thresh,:);
n2 = length(X2)
for i=1:n2
    a = X2(i,1);
    b = X2(i,2);
    g = X2(i,3);
    makeEllipse;
    name='PP';
    returnDist = findSSDistME(a,b,g,500);
    lightBlue = [91, 207, 244] / 255; 
    scatter(returnDist(1),returnDist(2),25,'r','filled')
    output_file = ['Phase Portraits 2/output_file_' name '_value' num2str(a,'%.2f') '_' num2str(b,'%.2f') '_' num2str(g,'%.2f') '.png']
    saveas(gcf,output_file,'png')
    close all
end