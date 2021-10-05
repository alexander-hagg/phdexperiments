clear predMap percImproved;
qdcfg.constraints.threshold = nan;
qdcfg.nGens = 200;
numModelSamples = [10 100 200 500 1000];
qdcfg.nChildren = 512;
for i=1:length(numModelSamples)
    models = output{1}.sail.modelPred;
    for m=1:2
        models{m}.trainInput = models{m}.trainInput(1:numModelSamples(i),:);
        models{m}.trainOutput = models{m}.trainOutput(1:numModelSamples(i),:);
    end
    
    tic;
    [predMap{i}, percImproved(i,:)] = createPredictionMap(models,qdcfg,d2,'featureRes',p.qd.data.predMapRes);
    timings(i) = toc;
end