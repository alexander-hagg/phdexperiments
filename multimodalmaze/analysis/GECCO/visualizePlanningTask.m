set(0,'DefaultFigureWindowStyle','docked');

for i=1
    uniqLbls = unique(runData{i}.classification.labels);
    clrs = hsv(length(uniqLbls));
    
    clrs = clrs(randperm(size(clrs,1)),:);
    clrs(1,:) = [1 1 1];
    
    figure(1);
    scatter(runData{i}.classification.simX(:,1),runData{i}.classification.simX(:,2),32, clrs(runData{i}.classification.labels,:), 'filled');
    
    figure(2);
    classMap = zeros(size(runData{i}.acqMap.fitness));
    classMap(~isnan(runData{i}.acqMap.fitness(:))) = runData{i}.classification.labels;
    viewMap(classMap,d)
    colormap(clrs);
    
    figure(3);
    viewMap(runData{i}.acqMap.ring1,d)
    figure(4);
    viewMap(runData{i}.acqMap.ring2,d)
    figure(5);
    viewMap(runData{i}.acqMap.ring3,d)
    drawnow;
end


