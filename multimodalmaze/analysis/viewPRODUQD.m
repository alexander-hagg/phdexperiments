maxValFit = 12;
maxValFirst = 3;
d = runData{1}.d;
p = runData{1}.p;
set(0,'DefaultFigureWindowStyle','docked');

for runID=1:length(runData)
    map = runData{runID}.acqMap;
    
    figure((runID-1)*3+1);
    [~,~,cHandle] = viewMap(map.firstGoal,d);caxis([0 maxValFirst]);colormap(parula(4));
    cHandle.Ticks = 0:3;
    title(['Selected Exit: ' int2str(p.selectionValue)]);
    
    figure((runID-1)*3+2);
    [~,~,cHandle] = viewMap(log(map.fitness),d);
    caxis([0 maxValFit]);colormap(parula(32));
    title(['Log Fitness']);
    
    figure((runID-1)*3+3);
    curiousness = map.curiousness;
    curiousness(isnan(map.fitness(:))) = nan;
    [~,~,cHandle] = viewMap(curiousness,d);caxis([-30 -20]);colormap(parula(8));
    title(['Curiousness']);
            
end

%%
for runID = 1:3
    map = runData{runID}.acqMap;
    figure((runID-1)*2+1);
    [~,~,cHandle] = viewMap(map.fitness,d);caxis([0 maxValFit]);
    figure((runID-1)*2+2);
    [~,~,cHandle] = viewMap(map.firstGoal,d);caxis([0 maxValFirst]);
    title(['Selected Exit: ' int2str(p.selectionValue)]);
    colormap(parula(4));
    cHandle.Ticks = 0:3;
end

%%

%%
runID = 1;

maxValFit = 12;
maxValFirst = 3;
d = runData{1}.d;
p = runData{1}.p;
set(0,'DefaultFigureWindowStyle','docked');

map = runData_EXTRA{runID}.allMaps{end}
figure(1);
[~,~,cHandle] = viewMap(log(map.fitness),d);caxis([0 maxValFit]);
title('Log Fitness after 1st MAP-Elites');
figure(2);
[~,~,cHandle] = viewMap(map.firstGoal,d);caxis([0 maxValFirst]);
title(['Selected Exit: ' int2str(p.selectionValue)]);
colormap(parula(4));
cHandle.Ticks = 0:3;
% Show classes
map.classes = nan(size(map.fitness));
map.classes(~isnan(map.fitness)) = runData{runID}.constraints.classLabels;
figure(3);
[~,~,cHandle] = viewMap(map.classes,d);caxis([0 max(map.classes(:))]);
title(['Classes']);
colormap(hsv(max(map.classes(:))));
%cHandle.Ticks = 0:1;

% Show selected invidiuals
sel = ismember(runData{runID}.constraints.classLabels,runData{runID}.constraints.selectedClasses);
map.selected = nan(size(map.fitness));
map.selected(~isnan(map.fitness)) = sel;
figure(4);
[~,~,cHandle] = viewMap(map.selected,d);caxis([0 1]);
title(['Selection']);
colormap(parula(2));
cHandle.Ticks = 0:1;

%
map = runData_EXTRA{runID+1}.allMaps{1}
figure(5);
[~,~,cHandle] = viewMap(log(map.fitness),d);caxis([0 maxValFit]);
title('Log Fitness 1st map after selection');
figure(6);
[~,~,cHandle] = viewMap(map.firstGoal,d);
caxis([0 maxValFirst]);
title(['Selected Exit: ' int2str(p.selectionValue)]);
colormap(parula(4));
cHandle.Ticks = 0:3;

map = runData_EXTRA{runID+1}.allMaps{end}
figure(7);
[~,~,cHandle] = viewMap(log(map.fitness),d);caxis([0 maxValFit]);
title('Log Fitness last map after selection');
figure(8);
[~,~,cHandle] = viewMap(map.firstGoal,d);caxis([0 maxValFirst]);
title(['Selected Exit: ' int2str(p.selectionValue)]);
colormap(parula(4));
cHandle.Ticks = 0:3;

