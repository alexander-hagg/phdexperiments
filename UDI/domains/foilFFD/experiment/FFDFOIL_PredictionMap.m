clear predMap percImproved;
qdcfg.constraints.threshold = nan;
qdcfg.nGens = 30;
nChildren = [64 128 256 512];

for i=1:length(nChildren)
    qdcfg.nChildren = nChildren(i);
    tic;
    [predMap{i}, percImproved(i,:)] = createPredictionMap(output{1}.sail.modelPred,qdcfg,d2,'featureRes',p.qd.data.predMapRes);
    timings(i) = toc;
end

%%
fig(1) = figure(1);
smoothData = smoothdata(percImproved','gaussian',10);
plot(smoothData,'LineWidth',2);
l = legend('64','128','256','512');
l.Title.String = 'No. Children';
grid on;
ylabel('Percentage of map improved');
xlabel('Iterations');
axis([0 1000 0 0.05]);


save_figures(fig, '.', 'PredictionMapImprovement', 14, [7 5]);

