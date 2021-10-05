
%%
set(0,'DefaultFigureWindowStyle','default');
d = domain_Maze;
angles = [30 90 150 210 270 330];
penaltyWeights = [0 0.5 1 2.5 5 10 25 50 75 100 150 200 250 500];
numA = length(angles);
exits = [1 2 3];

%load('GECCO2019/planner_parameterization/selectedMaps.mat');
load('GECCO2019/controller_parameterization/selectedMaps.mat');
data = data(:,angles);
disp('loaded');
%%

clear summedExits exitCorrectness exitInCorrectness medianVals prctileVals drift;
weightStrs = [];
weightOrIDs = [];
for weightID=1:length(penaltyWeights)
    weight = penaltyWeights(weightID)*1000+1;
    weightOrIDs(end+1) = weight;
    weightStrs(weight) = data(weight,1).description.penaltyweight;
    for angle=1:numA
        dataO = data(weight,angle).runData{end};
        if isempty(dataO.acqMap); continue; end
        exitCorrectness(weight,angle) = sum(dataO.acqMap.ring1(:)==dataO.p.selectionValue);
        exitInCorrectness(weight,angle) = sum(sum(dataO.acqMap.ring1(:)==setdiff(exits,dataO.p.selectionValue),2));
        drift(weight,angle) = nanmean(dataO.acqMap.penalty(:));
        
    end
end

exitCorrectness = exitCorrectness(weightOrIDs,:);
exitInCorrectness = exitInCorrectness(weightOrIDs,:);
drift = drift(weightOrIDs,:);


%% Correct
fig(1) = figure(1);hold off;
bpl = (100*(exitCorrectness-exitCorrectness(1,:))./900);
xx = repmat(penaltyWeights,1,6);
yy = reshape(bpl,6*size(penaltyWeights,2),1)';
pl(2) = semilogx(penaltyWeights,median(bpl'),'k:','LineWidth',2);hold on;
pl(3) = semilogx(penaltyWeights,mean(bpl'),'k-','LineWidth',2);hold on;
pl(1) = scatter(xx,yy,'filled');hold on;
grid on;
xlabel('Penalty Weight');
ylabel('Rel. Selected [%]');
axis([penaltyWeights(2) penaltyWeights(end) -5 50])
legend(pl,'Results','Median','Mean','Location','NorthWest');

% Incorrect
fig(2) = figure(2);hold off;
bpl = (100*(exitInCorrectness-exitInCorrectness(1,:))./900);
xx = repmat(penaltyWeights,1,6);
yy = reshape(bpl,6*size(penaltyWeights,2),1)';
pl(2) = semilogx(penaltyWeights,median(bpl'),'k:','LineWidth',2);hold on;
pl(3) = semilogx(penaltyWeights,mean(bpl'),'k-','LineWidth',2);hold on;
pl(1) = scatter(xx,yy,'filled');hold on;
grid on;
xlabel('Penalty Weight');
ylabel('Rel. Non-Selected [%]');

legend(pl,'Results','Median','Mean','Location','NorthWest');
axis([penaltyWeights(2) penaltyWeights(end) -50 20])

fig(3) = figure(3);hold off;
bpl = (drift)./900;
xx = repmat(penaltyWeights,1,6);
yy = reshape(bpl,6*size(penaltyWeights,2),1)';
pl(2) = semilogx(penaltyWeights,median(bpl'),'k:','LineWidth',2);hold on;
pl(3) = semilogx(penaltyWeights,mean(bpl'),'k-','LineWidth',2);hold on;
pl(1) = scatter(xx,yy,'filled');hold on;
grid on;
xlabel('Penalty Weight');
ylabel('User Selection Drift [%]');
legend(pl,'Results','Median','Mean','Location','NorthWest');
%axis([penaltyWeights(2) penaltyWeights(end) min(yy) max(yy)])
axis([penaltyWeights(2) penaltyWeights(end) min(yy) 6e-4])



%%
save_figures(fig, '.', 'parPenWeightController', 18, [6 5]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
grid on;
subplot(2,1,2);
boxplot((100*(exitCorrectness-exitCorrectness(1,:))./900)'); 
ylabel('Gain incorrect [%]');
ax = gca;
ax.XTickLabel = penaltyWeights;
grid on;
%axis([0.5 11.5 -50 10]);

%%
for weightID = 11:size(weightOrIDs,2)
angle = 3
weight = weightOrIDs(weightID)
runData = data(weight,angle).runData;

weights = reshape(runData{3}.acqMap.genes,900,14);
weights = weights(~any(isnan(weights')),:);
selection = 1:1:size(weights,1);
weights = weights(selection,:);

[~,values] = eval_maze_plan(weights,runData{3}.d);
selection = 1:length(values{2});
clrs = [0 0 0;1 0 0; 0 1 0; 0 0 1];
clrIDs = values{2}(selection)+1;
clrs = clrs(clrIDs,:);
vPlans(weights(selection,:),runData{3}.d,[],clrs);
title(['Weight: ' num2str(weightStrs(weight))])
drawnow
end
%%
for weightID=1:length(penaltyWeights)
%for weightID=[1]
    weight = penaltyWeights(weightID)*1000+1;
    for angle=5
        dataO = data(weight,angle).runData;
        figure(1);
        viewMap(dataO{3}.acqMap.ring1,dataO{3}.d);
        title(dataO{3}.d.penaltyWeight);
        colormap(parula(4))
        drawnow;
        pause(0.5)
    end
end


%%
xlabel('Penalty Weights');
ax = gca;
ax.XTick = penaltyWeights;
grid on;
axis([0.45 2.25 0 90]);ylabel('% Exits');
%legend([pl(1,1) pl(1,2)], 'Correct','Incorrect');

%%
save_figures(fig, '.', 'parPenWeight', 18, [10 6]);

