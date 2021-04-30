addpath('IFiles');
addpath('tFiles');
dI = dir('IFiles/*myFileI.csv');
dt = dir('tFiles/*myFilet.csv');
numFiles = size(dI,1);
tmax = 700;
tq = 0:0.05:tmax;
IMat = zeros(numFiles, length(tq));
n = 1;
bI = {};
bt = {};
parfor i=1:numFiles
    i
    bI{i} = csvread(dI(i).name);
    bt{i} = csvread(dt(i).name);
    [bt{i}, index] = unique(bt{i});
    interpI = interp1(bt{i}, bI{i}(index), tq, 'previous');
    bI{i} = {};
    bt{i} = {};
    IMat(i,:) = interpI';
end
Iplot = mean(IMat,1);
%stairs(tq, Iplot, 'b')
% avArray = sumArray/numFiles;
mean(Iplot(15000:end))
