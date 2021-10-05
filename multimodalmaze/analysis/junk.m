for runID=1:3
    disp(runID);
    d = runData{runID}.d;
    p = runData{runID}.p;
    for mapID=1:length(runData_EXTRA{runID}.allMaps)
        disp(mapID);
        map = runData_EXTRA{runID}.allMaps{mapID};
        genes = reshape(map.genes,900,size(map.genes,3));
        [~,values] = eval_maze(genes(~any(isnan(genes')),:),d.numHidden,d.maze,d.timesteps,d.ncores,d.goalLocations,d.useRNN,d.debug,false,true);
        runData_EXTRA{runID}.allMaps{mapID}.values = nan(size(map.fitness));
        runData_EXTRA{runID}.allMaps{mapID}.values(~any(isnan(genes'))) = values{2};
    end
end

%%
for runID=2:length(runData)
    d = runData{runID}.d;
    p = runData{runID}.p;
    for mapID=1:10:length(runData_EXTRA{runID}.allMaps)
        map = runData_EXTRA{runID}.allMaps{mapID};
        figure(1);
        viewMap(map.fitness,d)
        caxis([0 50]);
        colormap(parula(100));
        title(['SelVal: ' int2str(p.selectionValue)]);
        figure(2);
        viewMap(map.traject2Goal,d);
        caxis([1 3]);
        colormap(parula(3));
        title(['runID: ' int2str(runID) ', mapID: ' int2str(mapID)]);        
        drawnow;
        %incorrectSel((runID-1)*length(runData{runID}.allMaps)+mapID) = sum(values~=runData{1}.p.selectionValue);%/length(values);
        %pause(10.5)
    end
end
%%
firstGoal = [];runLengths = [0];
for runID=1:length(runData)
    d = runData{runID}.d;
    p = runData{runID}.p;
    runLengths(runID+1) = runLengths(runID) + length(runData_EXTRA{runID}.allMaps);
    for mapID=1:length(runData_EXTRA{runID}.allMaps)
        map = runData_EXTRA{runID}.allMaps{mapID};
        values = map.firstGoal; values = values(~isnan(values));
        firstGoal = [firstGoal;[sum(values==1),sum(values==2),sum(values==3)]];
    end
end
runLengths(1) = [];
fig(3+runData{1}.p.selectionValue) = figure(3+runData{1}.p.selectionValue); hold off;
plI = plot(firstGoal,'LineWidth',2);hold on;
plII = plot([runLengths;runLengths],[zeros(1,length(runLengths)); 200*ones(1,length(runLengths))],'--','Color',[0 0 0]);

exitID = int2str(runData{1}.p.selectionValue); if runData{1}.p.selectionValue == 9999; exitID='none';end
title(['Select elites that first pass exit ' exitID]);
grid on;axis([0 length(firstGoal) 0 200]);
xlabel('Generations'); ylabel('# Elites');
legend([plI(1) plI(2) plI(3) plII(4)],'Exit 1','Exit 2','Exit 3','Selection','Location','NorthWest');

%%
save_figures(fig, '.', 'correctElites2', 12, [7 7])


%% Show maps (fitness and values)
start =length(runData);
for runID=start:length(runData)
    fig(runData{1}.p.selectionValue) = figure(runData{1}.p.selectionValue);
    [figHandle, imageHandle, cHandle] = viewMap(runData{runID}.acqMap.traject2Goal,runData{runID}.d);
    caxis([1 3]);colormap(parula(3));
    exitID = int2str(runData{1}.p.selectionValue); if runData{1}.p.selectionValue == 9999; exitID='none';end
    title(['Select elites that first pass exit ' exitID]);
    cHandle.Ticks = [1 2 3];
    cHandle.Label.String = 'Goal #';
    drawnow;
    
    fig(runData{1}.p.selectionValue+3) = figure(runData{1}.p.selectionValue+3);
    [figHandle, imageHandle, cHandle] = viewMap(runData{runID}.acqMap.fitness,runData{runID}.d);
    caxis([0 50]);colormap(parula(32));
    title(['Select elites that first pass exit ' exitID]);
    cHandle.Label.String = 'Fitness';
    
    %%
    if runID > 1
        figID = runData{1}.p.selectionValue+6;
        numCells = runData{runID}.d.featureRes(1)*runData{runID}.d.featureRes(2);
        genes = runData{runID}.acqMap.genes;genes = reshape(genes,numCells,runData{runID}.d.dof);
        fitness = runData{runID}.acqMap.fitness;
        [outputs,confOutputs] = predictSimspace(genes(~isnan(fitness),:),runData{runID-1}.constraints.model);
        maxClassID = max(runData{runID-1}.classification.labels);
        colors = 0.7*ones(maxClassID,3);
        colors(runData{runID-1}.selection,:) = repmat([0 0 0],numel(runData{runID-1}.selection),1);
        [fig(runData{1}.p.selectionValue+6),pltHandle] = viewClasses(runData{runID-1}.classification,figID,colors);
        hold on;
        sPlt = scatter(outputs(:,1),outputs(:,2),16,[0 1 0],'s','filled');
        ax = gca;ax.XTick = [];ax.YTick = [];
        legend([pltHandle sPlt],'Classes (selected dark)','New Samples');
        hold off;
    end
    drawnow;
    %%
    pause(1);
    
end
%%
save_figures(fig, '.', 'maps', 12, [7 7])



%%
save_figures(figHandles, '.', 'class', 12, [4 4])

%% Simulate random solutions from map or sample vector
runID = 1;
d = runData{runID}.d;
%for iRobot=1:size(runData{runID}.samples,1)
%while true
timesteps = 3000;
iRobot = randi(size(runData{runID}.classification.X,1))
%
[trajectories,values] = eval_maze(runData{runID}.classification.X(iRobot,:,:),d.numHidden,d.maze,timesteps,d.ncores,d.goalLocations,d.useRNN,d.debug,true,true);

%
trajectory = squeeze(trajectories(1,:,:));
figure(1);hold off;
plotSim(trajectory,d);
title(iRobot);
drawnow;
%pause(2);
%end

%% Histograms
figure(98);clf

if runData{1}.p.selectionCriterionID==1
    edges = [0.0 0.1 0.9 1.0];
elseif runData{1}.p.selectionCriterionID==3
    edges = 1:1:4;
end

valueFcn = {'evalAlignment','evalTrajectoryClosestToGoal','evalFirstGoal'};
tic
clear h valueT;
for runID=1:length(runData)
    for ii=1:size(runData{runID}.acqMap.phenotype,1)
        valueT{runID}(ii) = feval(valueFcn{runData{1}.p.selectionCriterionID},squeeze(runData{runID}.acqMap.phenotype(ii,:,:)),runData{runID}.d.goalLocations);
        h{runID} = histcounts(valueT{runID},edges)';
    end
end
toc
%%
bar(edges(1:end-1),[h{:}])
xlabel('Selection Value');
ylabel('Count');
title(['Selected: ' num2str(runData{1}.p.selectionValue)]);
axis([0.5 3.5 0 350]);