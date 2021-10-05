
d = runData{2}.d;
acqMap = runData{2}.acqMap;
classification = runData{2}.classification;
samples = reshape(acqMap.genes,d.featureRes(1)*d.featureRes(2),d.dof); samples = samples(~any(isnan(samples)'),:); % Filter out NaNs
d.metricFitness               = @metricFitness;
d.fitfun                      = @(X) objective(X, d.evalFcn, d.metricFitness, [], d.penaltyWeight);

[fitness, phenotype, values] = d.fitfun(samples);


%%
figure(1);viewMap(acqMap.ring1,d); title('ring1');
figure(2);viewMap(acqMap.orgFitness,d); title('org. fitness'); caxis([0 200]);
figure(3);viewMap(acqMap.penalty,d); title('penalty');
figure(4);viewMap(acqMap.fitness,d); title('adj. fitness'); caxis([0 500]);

acqMap.newfitness = acqMap.fitness;
acqMap.newfitness(~isnan(acqMap.newfitness)) = values{end-1};
figure(5);viewMap(acqMap.newfitness,d); title('new fitness'); caxis([0 200]);

%%
figure(6);
ring1sc = acqMap.ring1(:);
ring1sc = ring1sc(~isnan(ring1sc));
clrs = parula(length(unique(ring1sc)));
scatter(classification.simX(:,1),...
    classification.simX(:,2),...
    32,...
    clrs(ring1sc+1,:),...
    'filled');

figure(7);
vPlans(samples,d,6)


