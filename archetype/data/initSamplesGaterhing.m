clear;clc;
load('oldinitsamples64.mat');

initSet.samples = initSamples;
initSet.polygons = d.getPhenotype(initSamples);
initSet.shapeFeatures = d.categorize([],initSet.polygons,p,d);
initSet.flowFeatures = [movMeanE;movMeanMaxU]';
initSet.features = [initSet.shapeFeatures,initSet.flowFeatures];
initSet.fitness = movMeanMaxU';

save('initsamples64.mat','initSet');