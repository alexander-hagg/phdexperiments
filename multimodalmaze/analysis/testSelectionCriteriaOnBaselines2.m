set(0,'DefaultFigureWindowStyle','docked');


for runID=1:3
    d = runData{runID}.d;
    p = runData{runID}.p;
    map = runData{runID}.acqMap;
    
    figure(1);
    viewMap(map.firstGoal,d);
    caxis([-0.5 3.5]);colormap(parula(4));
    title(['Selval ' int2str(p.selectionValue) ', run ' int2str(runID)]);
    
    
    values = [];firstGoal = [];
    for i=1:size(map.phenotype,1)
        firstGoal(i) = evalFirstGoal(squeeze(map.phenotype(i,:,:)),d.goalLocations);
    end
    map.firstGoal(~isnan(map.firstGoal)) = firstGoal;
    
    %
    values{1} = map.alignment(~isnan(map.fitness(:)));
    values{2} = map.traject2Goal(~isnan(map.fitness(:)));
    values{3} = map.firstGoal(~isnan(map.fitness(:)));
    classification = runData{runID}.classification;
    minSelectedSamples = 100;
    thresholds = 1:-0.01:0.01;
    for selThresID = 1:length(thresholds)
        selThresh = thresholds(selThresID);
        [isSelected,selection] = mazeSelection(classification.labels, values, ...
            3, p.selectionValue, selThresh);
        disp(['Selected ' int2str(sum(isSelected)) ' individuals.']);
        if (selThresh < p.selectionThreshold) && sum(isSelected) > minSelectedSamples
            disp(['Threshold at: ' num2str(selThresh)]);
            break;
        end
    end
    
    % Show selected invidiuals
    sel = ismember(classification.labels,selection);
    selected = nan(size(map.fitness));
    selected(~isnan(map.fitness)) = sel;
    
    %map.selected = nan(size(map.fitness));
    %map.selected(~isnan(map.fitness)) = runData{runID}.isSelected;
    %selected = map.selected;
    figure(2);
    viewMap(selected,d);
    caxis([0 1]);colormap(parula(2));
    title(['Selected']);
    
    figure(3);
    curiousness = map.curiousness;
    curiousness(isnan(map.fitness)) = nan;
    viewMap(curiousness,d);
    caxis([-20 0]);colormap(parula(20));
    title(['Curiousness']);
    
    figure(4);
    viewMap(log(map.fitness),d);
    caxis([0 12]);colormap(parula(12));
    title(['Log(fitness)']);
    
    drawnow;
    pause(0.5);
end