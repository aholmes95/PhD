function customCMap = createCustomCMap(data)

L=nnz(data);
indexValue = 0;

topColor = [1 0 0];
indexColor = [1 1 1];
bottomcolor = [0 0 1];

largest = max(max(data));
smallest = min(min(data));
index = L*abs(indexValue-smallest)/(largest-smallest);

% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
            linspace(bottomcolor(2),indexColor(2),100*index)',...
            linspace(bottomcolor(3),indexColor(3),100*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
            linspace(indexColor(2),topColor(2),100*(L-index))',...
            linspace(indexColor(3),topColor(3),100*(L-index))'];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
% colormap(customCMap);
% psudo = pcolor(data);
% cb = colorbar;
% cb.Ruler.Scale = 'log';
% cb.Ruler.MinorTick = 'on';