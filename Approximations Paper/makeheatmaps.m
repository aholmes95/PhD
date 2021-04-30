alphaVec = linspace(0.4,0.8,50);
betaVec = linspace(0.25,0.8,50);
gammaVec = linspace(0.1,0.5,50);
h = heatmap(round(betaVec,3),round(alphaVec,3),ratmat)
h.ColorScaling = 'log'
custommap = createCustomCMap(log(ratmat));
h.Colormap = custommap
myvec = (mod(1:50,5)==1);
h.XDisplayLabels(~myvec) = {''}
h.YDisplayLabels(~myvec) = {''}
h.YDisplayLabels(50) = {0.8}
h.XDisplayLabels(50) = {0.5}
h.NodeChildren(3).YDir='normal'; 
xlabel('Gamma')
ylabel('Beta')
% set(gca,'fontsize', 18)
