clear;clc;
DOF = 16;
DOMAIN = 'footprints';
addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));

d = domain(DOF);
p = defaultParamSet;
experimentName = 'SobolPredictionSymmetry'
numInitSamples = 64;
numReplicates = 64;
%% ----------------------------------------------------------------------------------
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,features] = d.fitfun(initSamples);
%%  
clear rmse;
cvIndices = crossvalind('Kfold',size(initSamples,1),numReplicates);
    
for rep=1:numReplicates
    disp([int2str(rep) '/' int2str(numReplicates)]);
    pModel = paramsGP(d.dof);
    fitModel = trainGP(initSamples(cvIndices~=rep,:),fitness(cvIndices~=rep),pModel);
    prediction = predictGP(fitModel, initSamples(cvIndices==rep,:));
    trueFit = fitness(cvIndices==rep);    
    rmse(1,rep) = sqrt(immse(prediction(:,1),trueFit))./(trueFit);
    
    featModel1 = trainGP(initSamples(cvIndices~=rep,:),features(cvIndices~=rep,1),pModel);
    prediction = predictGP(featModel1, initSamples(cvIndices==rep,:));
    trueFeat = features(cvIndices==rep,1);    
    rmse(2,rep) = sqrt(immse(prediction(:,1),trueFeat))./(trueFeat);

    featModel2 = trainGP(initSamples(cvIndices~=rep,:),features(cvIndices~=rep,2),pModel);
    prediction = predictGP(featModel2, initSamples(cvIndices==rep,:));
    trueFeat = features(cvIndices==rep,2);    
    rmse(3,rep) = sqrt(immse(prediction(:,1),trueFeat))./(trueFeat);    
end

save([experimentName '.mat']);

%%
fig(1) = figure(1);
boxplot(100*rmse');
grid on;
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'Symmetry', 'Area', 'Perimeter Length'};
ylabel('relative RMSE (%)');
title('Point Symmetry - leave one out crosscorr. ');
ax.YLim = [0 100];
save_figures(fig, '.', [experimentName '_errors'], 12, [5 4]);


