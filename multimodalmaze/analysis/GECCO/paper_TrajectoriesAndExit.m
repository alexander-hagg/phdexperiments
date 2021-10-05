runID = 1;
acqMap = runData{runID}.acqMap;
d = runData{runID}.d;

fig(1) = figure(1);
[~,~,cHandle] = viewMap(acqMap.ring1, runData{runID}.d, acqMap.edges,'flip');
colormap(parula(4));
caxis([0 4]); 
cHandle.Ticks = 0.5:4.5;
cHandle.TickLabels = [0:3];
cHandle.Label.String = 'Exit Taken';

fig(2) = figure(2);clf;
samples = reshape(acqMap.genes,900,size(acqMap.genes,3));
samples = samples(~any(isnan(samples)'),:);
subsamples = samples(1:9:size(samples,1),:);
vPlans(subsamples,d,2);
grid on;

save_figures(fig(1), '.', 'exampleExits2', 16, [4 4]);
save_figures(fig(2), '.', 'exampleTraj2', 16, [5 4]);


