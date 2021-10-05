set(0,'DefaultFigureWindowStyle','default');
d = domain_Maze;
angles = [30 90 150 210 270 330];
numA =length(angles);
mutationRates = [0.001 0.005 0.01 0.015 0.02 0.05 0.08];
methods = {'UDHM','SEED','COMB'};
exits = [1 2 3];
%load('/scratch/ahagg2s/GECCO2019/controller_MUTATION_4096/selectedMaps.mat');
load('/scratch/ahagg2s/GECCO2019/planner_MUTATION_4096/selectedMaps.mat');
%data = data(:,1000*mutationRates,angles,:);
%fdata = data;
%fdata(cellfun(@isempty,{data(:,:,:,:).fname})) = [];
%fdata = reshape(fdata,size(data,1),size(data,2),size(data,3),[]);

fdataRunData = dataRunData;
fdataRunData(cellfun(@isempty,{dataRunData(:,:,:,:).runData})) = [];
fdataRunData = reshape(fdataRunData,size(dataRunData,1),size(dataRunData,2),size(dataRunData,3),[]);
disp('loaded');

%%
%for weightID = length(penaltyWeights)-1;%1:length(penaltyWeights)
angle = 6;
methodID = 2;
mutRateID = 1;
exitID = 3;
%runData = data(weight,angle).runData;
runData = fdataRunData(methodID,mutRateID,angle,exitID).runData;

weights = reshape(runData{1}.acqMap.genes,900,size(runData{1}.acqMap.genes,3));
%[~,values] = eval_maze_plan(weights,runData{1}.d);
[~,values] = eval_maze(weights,runData{1}.d);

map1 = runData{1}.acqMap.ring1;
map1(1:end) = values{2};
map1(isnan(runData{1}.acqMap.fitness)) = nan;

weights = reshape(runData{end}.acqMap.genes,900,size(runData{end}.acqMap.genes,3));
%[~,values] = eval_maze_plan(weights,runData{end}.d);
[~,values] = eval_maze(weights,runData{end}.d);
map2 = runData{3}.acqMap.ring1;
map2(1:end) = values{2};
map2(isnan(runData{end}.acqMap.fitness)) = nan;

disp('done');

%%
angle = 6;
methodID = 1;
mutRateID = 3;
exitID = 1;
runData = fdataRunData(methodID,mutRateID,angle,exitID).runData;

color = repmat([1 1 1],4,1);
color(2,:) = [0.8 0.8 0.8];
color(3,:) = [0.8 0.8 0.8];
color(4,:) = [0.8 0.8 0.8];
color(runData{3}.p.selectionValue+1,:) = [0.1 0.7 0.1];
%

map1 = runData{1}.acqMap.ring1;
map1(isnan(runData{1}.acqMap.fitness(:))) = nan;
fig(1) = figure(1);
[~,~,cHandle] = viewMap(map1,runData{1}.d);
colormap(color); 
%cHandle.Ticks = 0:2;
%cHandle.TickLabels = 0:2;
%cHandle.Limits = [0 2.25]
cHandle.Label.String = ['Exit' int2str(runData{3}.p.selectionValue)];

%
map2 = runData{3}.acqMap.ring1;
map2(isnan(runData{3}.acqMap.fitness(:))) = nan;
fig(2) = figure(2);
%cHandle.Limits = [0 2]
[~,~,cHandle] = viewMap(map2,runData{3}.d);
colormap(color); 
%cHandle.Ticks = 0:2;
%cHandle.TickLabels = 0:2;
%cHandle.Limits = [0 2.25]
cHandle.Label.String = ['Exit' int2str(runData{3}.p.selectionValue)];

drawnow;

%%
fig(3) = figure(3);
simX = runData{1}.classification.simX;
%clrs = parula(4);
mmm = map1(:);
scatter(simX(:,1),simX(:,2),16,color(map1(~isnan(mmm))+1,:),'filled');
ax = gca;
ax.XTick = -50:25:50;
ax.YTick = -50:25:50;
%axis([-50 50 -50 50]);
%save_figures(fig, '.', 'mapExamples', 18, [6 5]);

fig(4) = figure(4);
simX = runData{3}.classification.simX;
clrs = parula(4);
scatter(simX(:,1),simX(:,2),16,color(map2(~isnan(map2(:)))+1,:),'filled');
ax = gca;
ax.XTick = -50:25:50;
ax.YTick = -50:25:50;
%axis([-50 50 -50 50]);
save_figures(fig, '.', 'mapExamples_CONTROLLER', 18, [6 5]);

%% Trajectories before selection


samples = reshape(runData{1}.acqMap.genes,900,size(runData{1}.acqMap.genes,3));
nonans = any(~isnan(samples'));
samples = samples(nonans,:);
selection = 1:10:size(samples,1);
samples = samples(selection,:);
values = runData{1}.acqMap.ring1(nonans);
%clrs = parula(4);
clrIDs = values(selection)+1;
clrs = color(clrIDs,:);
fig(1) = vPlans(samples,runData{3}.d,1,clrs);
grid on;
ax = gca;
ax.XTick = 0:400/30:400;ax.XTickLabel = [];
ax.YTick = 0:400/30:400;ax.YTickLabel = [];

% Trajectories after selection
samples = reshape(runData{3}.acqMap.genes,900,size(runData{3}.acqMap.genes,3));
nonans = any(~isnan(samples'));
samples = samples(nonans,:);
selection = 1:10:size(samples,1);
samples = samples(selection,:);
values = runData{3}.acqMap.ring1(nonans);
%clrs = parula(4);
clrIDs = values(selection)+1;
clrs = color(clrIDs,:);
fig(2) = vPlans(samples,runData{3}.d,2,clrs);
grid on;
ax = gca;
ax.XTick = 0:400/30:400;ax.XTickLabel = [];
ax.YTick = 0:400/30:400;ax.YTickLabel = [];

%%
fig(1) = figure(1);
fig(2) = figure(2);
save_figures(fig, '.', 'pathplannerExamples_CONTROLLER', 18, [6 5]);


%end


