clear;clc;
DOF = 16;
DOMAIN = 'footprints';
addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));

d = domain(DOF);
d.fitfun = d.fitfunLettuce;
p = defaultParamSet;
load('initsamples64.mat');
experimentName = 'SobolPredictionLettuce'
numInitSamples = size(initSet.samples,1);
numReplicates = 10%numInitSamples;
%% ----------------------------------------------------------------------------------
clear mape;
cvIndices = crossvalind('Kfold',size(initSet.samples,1),numReplicates);
    
for rep=1:numReplicates
    disp([int2str(rep) '/' int2str(numReplicates)]);
    pModel = paramsGP(d.dof);
    
    featModel1 = trainGP(initSet.samples(cvIndices~=rep,:),initSet.features(cvIndices~=rep,1),pModel);
    prediction = predictGP(featModel1, initSet.samples(cvIndices==rep,:));
    trueFeat = initSet.features(cvIndices==rep,1);    
    mape(1,rep) = mean(abs(prediction(:,1)-trueFeat)./mean(trueFeat))*100;

    fitModel = trainGP(initSet.samples(cvIndices~=rep,:),initSet.fitness(cvIndices~=rep),pModel);
    prediction = predictGP(fitModel, initSet.samples(cvIndices==rep,:));
    trueFit = initSet.fitness(cvIndices==rep)';    
    mape(2,rep) = mean(abs(prediction(:,1)-trueFeat)./mean(trueFeat))*100;
    
end
%%
fig(1) = figure(1);hold off;
boxplot(mape');
grid on;
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'E', 'max U'};
ylabel('MAPE');
title('Lettuce - crossvalidation (10 replicates)');
ax.YLim = [0 100];
save_figures(fig, '.', [experimentName '_errors'], 12, [5 4]);




