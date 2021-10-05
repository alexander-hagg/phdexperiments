set(0,'DefaultFigureWindowStyle','docked');
testConstraints = false;
testPenalty = false;
maxFitVal = 1000;

for runID=1
    d = runData{runID}.d;
    p = runData{runID}.p;
    map = runData{runID}.acqMap;
    
    
    %%
    values = [];firstGoal = [];
    for i=1:size(map.phenotype,1)
        firstGoal(i) = evalFirstGoal(squeeze(map.phenotype(i,:,:)),d.goalLocations);
    end
    map.firstGoal(~isnan(map.firstGoal)) = firstGoal;
    
    fig(1) = figure(1);
    [~,~,cRef] = viewMap(map.firstGoal,d);
    caxis([-0.5 3.5]);colormap(parula(4));
    cRef.Ticks = 0:4;
    cRef.Label.String = 'Exit taken';
    title(['Selval ' int2str(p.selectionValue) ', run ' int2str(runID)]);
    
    % Show selected invidiuals
    selection = ismember(runData{runID}.classification.labels,runData{runID}.constraints.selectedClasses);
    selected = nan(size(map.fitness));
    selected(~isnan(map.fitness)) = selection;
    
    fig(2) = figure(2);
    viewMap(selected.*map.firstGoal,d);
    caxis([0 3]);colormap(parula(4));
    title(['Selected']);
    
    fig(3) = figure(3);
    %curiousness = map.curiousness;
    %curiousness(isnan(map.fitness)) = nan;
    %viewMap(curiousness,d);
    %caxis([-50 0]);colormap(parula(50));
    %title(['Curiousness']);
    map2 = runData{runID+1}.acqMap;
    [~,~,cRef] = viewMap(map2.firstGoal,d);
    caxis([-0.5 3.5]);colormap(parula(4));
    cRef.Ticks = 0:4;
    cRef.Label.String = 'Exit taken';
    title(['Selval ' int2str(p.selectionValue) ', run ' int2str(runID+1)]);
    
    fig(4) = figure(4);
    viewMap(map2.fitness,d);
    caxis([0 maxFitVal]);colormap(parula(10));
    title(['Fitness']);
    
    if testPenalty
        fig(5) = figure(5);
        mapNew = runData{runID+1}.acqMap;
        genes = reshape(mapNew.genes,900,size(mapNew.genes,3)); genes(isnan(mapNew.fitness(:)),:) = [];
        pW = d.penaltyWeight;
        pW = 1;
        penalty = constraintPenalty(genes,runData{runID}.constraints,pW);
        mapNew.penalty = nan(size(mapNew.fitness));
        mapNew.penalty(~isnan(mapNew.fitness)) = penalty;
        viewMap(mapNew.penalty,d);
        caxis([0 1]);colormap(parula(25));
        title(['Penalty next map']);
    end
    
    %%
    if testConstraints
        fig(6) = figure(6);
        d.penaltyWeight = 5;
        maxFitVal = 200;
        d.fitfun  = @(X) objective(X, d.evalFcn, d.metricFitness, 'none', d.penaltyWeight);
        %t_noconstraint_fitness = d.fitfun(genes);
        noconstraint_fitness = nan(size(mapNew.fitness));
        noconstraint_fitness(~isnan(mapNew.fitness)) = t_noconstraint_fitness;
        viewMap(noconstraint_fitness,d);
        caxis([0 maxFitVal]);colormap(parula(25));
        title(['Fitness next map']);
        
        fig(7) = figure(7);
        d.fitfun  = @(X) objective(X, d.evalFcn, d.metricFitness, runData{runID}.constraints, d.penaltyWeight);
        %t_constraint_fitness = d.fitfun(genes);
        constraint_fitness = nan(size(mapNew.fitness));
        constraint_fitness(~isnan(mapNew.fitness)) = t_constraint_fitness;
        viewMap(constraint_fitness,d);
        caxis([0 maxFitVal]);colormap(parula(25));
        title(['Fitness next map']);
    end
    
    
    
    drawnow;
    pause(0.5);
end
%%
save_figures(fig, '.', 'selectionPenalty', 12, [5 4])
