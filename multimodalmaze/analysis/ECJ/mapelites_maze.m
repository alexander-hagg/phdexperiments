set(0,'DefaultFigureWindowStyle','docked');
testConstraints = false;
testPenalty = false;
maxFitVal = 50;

runID=1
d = runData{runID}.d;
p = runData{runID}.p;
map = runData{runID}.acqMap;

fig(1) = figure(1);
[~,~,cRef] = viewMap(map.firstGoal,d);
caxis([-0.5 3.5]);colormap(parula(4));
cRef.Ticks = 0:4;
cRef.Label.String = 'Exit taken';
title(['Selval ' int2str(p.selectionValue) ', run ' int2str(runID)]);

fig(2) = figure(2);
firstGoal = runData{1}.acqMap.firstGoal(:);
firstGoal(isnan(firstGoal)) = [];
vTrajectories(runData{1}.acqMap.phenotype,firstGoal,d.maze,true,0);
%%
fig(3) = figure(3);
firstGoal = runData{1}.acqMap.firstGoal(:);
firstGoal(isnan(firstGoal)) = [];
selection = [100 200 300 400 500 600 700];
vTrajectories(runData{1}.acqMap.phenotype(selection,:,:),firstGoal(selection),d.maze,true,0);
grid on;
ax = gca;
ax.XTick = linspace(0,400,30);
ax.YTick = linspace(0,400,30);


%%
save_figures(fig, '.', 'map_64k', 12, [5 4])