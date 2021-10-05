clear;clc;
DOF = 16;
DOMAIN = 'footprints';
QD = 'grid';
SAQD = 'sail';

addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath(genpath('QD/grid')); rmpath(genpath('QD/voronoi')); addpath(genpath(['QD/' QD]));
rmpath(genpath('QD/sail')); rmpath(genpath('QD/sphen')); addpath(genpath(['QD/' SAQD]));

d = domain(DOF);
p = defaultParamSet;
p.intermediateMaps = true;

p.infill = infillParamSet;
p.infill.nTotalSamples = 1024;%p.nGens*p.nChildren;
p.infill.nAdditionalSamples = p.nChildren;

sailP = p;
sailP.nGens = 2^9;
sailP.nChildren = 2^6;

pMAPE = p;
pMAPE.nGens = 2^12;

d.fitfun = d.fitfunPointSymm; % Multimodal function: point symmetry (center); for testing purposes
experimentName = 'SAILvsMAPELITES_pointSymmetry'
numReplicates = 5;
%% ----------------------------------------------------------------------------------

for rep=1 : numReplicates
    sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = (rep-1)*p.numInitSamples + 1;
    initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
    [fitness,polygons] = d.fitfun(initSamples);
    features = d.categorize(initSamples, polygons, p, d);
    
    initmap{rep}                                            = createMap(d, p);
    [replaced, replacement, features]                       = nicheCompete(initSamples, fitness, initmap{rep}, d, p, features);
    initmap{rep}                                            = updateMap(replaced,replacement,initmap{rep},fitness,initSamples,features);
    
    %% MAP-Elites
    [mapMAPELITES{rep}, ~, ~, allMaps{rep}] = illuminate(initmap{rep},d.fitfun,pMAPE,d);
    
    %% SAIL (Surrogate-Assisted Illumination)
    SAQD = 'sail';
    rmpath(genpath('QD/sail')); rmpath(genpath('QD/SPHEN')); addpath(genpath(['QD/' SAQD]));
    
    [mapSAIL{rep},surrogate{rep},allMapsSAIL{rep}] = sail(initmap{rep},sailP,d);
    
    %% SPHEN (Surrogate-Assisted Phenotypic Niching
    SAQD = 'sphen';
    rmpath(genpath('QD/sail')); rmpath(genpath('QD/sphen')); addpath(genpath(['QD/' SAQD]));
    
    [mapSPHEN{rep},surrogateFitnessSPHEN{rep},surrogateFeaturesSPHEN{rep},allMapsSPHEN{rep}] = sail(initmap{rep},sailP,d);
    
end

save([experimentName '_4.mat']);
%% Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get true fitness for SAIL
for rep=1:numReplicates
tGenes = reshape(mapSAIL{rep}.genes,[],16);
tFitness = d.fitfun(tGenes);
mapSAIL{rep}.trueFitness = reshape(tFitness,p.resolution,p.resolution);
for i=1:length(allMapsSAIL{rep})
    disp([int2str(i) '/' int2str(length(allMapsSAIL{rep}))]);
    tGenes = reshape(allMapsSAIL{rep}{i}.genes,[],16);
    tFitness = d.fitfun(tGenes);
    allMapsSAIL{rep}{i}.trueFitness = reshape(tFitness,p.resolution,p.resolution);
end
end
%% SPHEN
for rep=1:numReplicates

tGenes = reshape(mapSPHEN{rep}.genes,[],16);
tFitness = d.fitfun(tGenes);
mapSPHEN{rep}.trueFitness = reshape(tFitness,p.resolution,p.resolution);

tGenes = reshape(mapSPHEN{rep}.genes,[],16); tGenes(any(isnan(tGenes')),:) = [];
[tFitness,tpolygons] = d.fitfun(tGenes);
trueFeatures = d.categorize(tGenes, tpolygons, p, d);
truemapSPHEN{rep}                                       = createMap(d, p);
[replaced, replacement, features]                       = nicheCompete(tGenes, tFitness, truemapSPHEN{rep}, d, p, trueFeatures);
truemapSPHEN{rep}                                       = updateMap(replaced,replacement,truemapSPHEN{rep},tFitness,tGenes,trueFeatures);
mapSPHEN_predFeatures{rep} = mapSPHEN{rep}.features(:); mapSPHEN_predFeatures{rep} = mapSPHEN_predFeatures{rep}(~isnan(mapSPHEN_predFeatures{rep}));
mapSPHEN_trueFeatures{rep} = trueFeatures(:);

for i=1:length(allMapsSPHEN{rep})
    disp([int2str(i) '/' int2str(length(allMapsSPHEN{rep}))]);
    tGenes = reshape(allMapsSPHEN{rep}{i}.genes,[],16); tGenes(any(isnan(tGenes')),:) = [];
    [tFitness,tpolygons] = d.fitfun(tGenes);
    trueFeatures = d.categorize(tGenes, tpolygons, p, d);
    allMapsSPHENTRUE{rep}{i}                            = createMap(d, p);
    [replaced, replacement, features]                   = nicheCompete(tGenes, tFitness, allMapsSPHENTRUE{rep}{i}, d, p, trueFeatures);
    allMapsSPHENTRUE{rep}{i}                            = updateMap(replaced,replacement,allMapsSPHENTRUE{rep}{i},tFitness,tGenes,trueFeatures);
    %allMapsSAIL{i}.trueFitness = reshape(tFitness,p.resolution,p.resolution);
end
end
%% Reporting RMSE
for rep=1:numReplicates
disp(rep)
for i=1:length(allMapsSPHEN{rep})
    disp('SPHEN')
    tGenes = reshape(allMapsSPHEN{rep}{i}.genes,[],16);
    predFeatures = reshape(allMapsSPHEN{rep}{i}.features,[],2); predFeatures(any(isnan(tGenes')),:) = [];
    predictedfits = reshape(allMapsSPHEN{rep}{i}.fitness,[],1);
    predictedfits(any(isnan(tGenes'))) = [];
    tGenes(any(isnan(tGenes')),:) = [];
    [tFitness,tpolygons] = d.fitfun(tGenes);
    trueFeatures = d.categorize(tGenes, tpolygons, p, d);
    
    rmseFEAT_SPHEN1(i,rep) = sqrt(immse(predFeatures(:,1),trueFeatures(:,1)));
    rmseFEAT_SPHEN2(i,rep) = sqrt(immse(predFeatures(:,2),trueFeatures(:,2)));
    rmseFIT_SPHEN(i,rep) = sqrt(immse(predictedfits,tFitness));
end

for i=1:length(allMapsSAIL{rep})
    disp('SAIL')
    tGenes = reshape(allMapsSAIL{rep}{i}.genes,[],16);
    predFitnessSAIL = reshape(allMapsSAIL{rep}{i}.fitness,[],1);
    predFitnessSAIL(any(isnan(tGenes'))) = [];
    tGenes(any(isnan(tGenes')),:) = [];
    [tFitness,~] = d.fitfun(tGenes);
    rmseFIT_SAIL(i,rep) = sqrt(immse(predFitnessSAIL,tFitness));
end
end

%%
fig(1) = figure(1);hold off;
plot(mean(rmseFIT_SAIL(1:end-1,:)'),'LineWidth',2);
hold on;
plot(mean(rmseFIT_SPHEN(1:end-1,:)'),'LineWidth',2);
plot(mean(rmseFEAT_SPHEN1(1:end-1,:)'),'LineWidth',2);
plot(mean(rmseFEAT_SPHEN2(1:end-1,:)'),'LineWidth',2);

legend('SAIL Fitness', 'SPHEN Fitness', 'SPHEN Area', 'SPHEN Perimeter');
grid on;
ylabel('RMSE');
ax = gca;
ax.XTick = [64/8 256/8 512/8 1024/8];
ax.XTickLabel = ax.XTick*p.nChildren;
xlabel('Precise Evaluations');
ax.YAxis.Limits = [0 0.3];
save_figures(fig, '.', [experimentName '_errors'], 12, [5 4]);


%% Visualization
set(0,'DefaultFigureWindowStyle','default')
rep = 1;

clear fig; close all;
fig(1) = figure(1);
[~,~,clabel] = viewMap(allMaps{rep}{128},d,'fitness');
clabel.Label.String = 'Fitness';
caxis([0 1]);
title(['MAP-Elites ' int2str(128*pMAPE.nChildren) ' PE']);


fig(2) = figure(2);
[~,~,clabel] = viewMap(allMaps{rep}{end},d,'fitness');
clabel.Label.String = 'Fitness';
caxis([0 1]);
title(['MAP-Elites ' int2str(pMAPE.nGens*pMAPE.nChildren) ' PE']);


fig(3) = figure(3);
[~,~,clabel] = viewMap(mapSAIL{rep},d,'fitness');
clabel.Label.String = 'Predicted Fitness';
caxis([0 1]);
title(['SAIL ' int2str(128*pMAPE.nChildren) ' PE']);


fig(4) = figure(4);
[~,~,clabel] = viewMap(mapSAIL{rep},d,'trueFitness');
clabel.Label.String = 'True Fitness';
caxis([0 1]);
title(['SAIL ' int2str(128*pMAPE.nChildren) ' PE']);

fig(5) = figure(5);
[~,~,clabel] = viewMap(mapSPHEN{rep},d,'fitness');
clabel.Label.String = 'Predicted Fitness';
caxis([0 1]);
title(['SPHEN ' int2str(128*pMAPE.nChildren) ' PE']);

fig(6) = figure(6);hold off;
[~,~,clabel] = viewMap(truemapSPHEN{rep},d,'fitness');
clabel.Label.String = 'True Fitness and Features';
caxis([0 1]);
title(['SPHEN ' int2str(128*pMAPE.nChildren) ' PE']);


save_figures(fig, '.', experimentName, 12, [5 4]);



%%
clear mape_medians SPHEN_mediansTrue sail_mediansTrue
for rep=1:numReplicates
    disp(rep);

for i=1:pMAPE.nGens
    mape_medians(i,rep) = nanmedian(allMaps{rep}{i}.fitness(:));
end
for i=1:length(allMapsSAIL{rep})
    %sail_medians(p.nGens-length(allMapsSAIL) +i) = nanmedian(allMapsSAIL{i}.fitness(:));
    sail_mediansTrue(p.nGens-length(allMapsSAIL{rep}) +i,rep) = nanmedian(allMapsSAIL{rep}{i}.trueFitness(:));
end
for i=1:length(allMapsSPHENTRUE{rep})
    SPHEN_mediansTrue(p.nGens-length(allMapsSPHENTRUE{rep}) +i,rep) = nanmedian(allMapsSPHENTRUE{rep}{i}.fitness(:));
end
end

clear mape_bins sail_bins SPHEN_bins
for rep=1:numReplicates
disp(rep);
for i=1:pMAPE.nGens
    mape_bins(i,rep) = 100*sum(~isnan(allMaps{rep}{i}.fitness(:)))./numel(allMaps{rep}{i}.fitness);
end
for i=1:length(allMapsSAIL{rep})
    sail_bins(p.nGens-length(allMapsSAIL{rep}) +i,rep) = 100*sum(~isnan(allMapsSAIL{rep}{i}.trueFitness(:)))./numel(allMapsSAIL{rep}{i}.trueFitness);
end
for i=1:length(allMapsSPHENTRUE{rep})
    SPHEN_bins(p.nGens-length(allMapsSPHENTRUE{rep}) +i,rep) = 100*sum(~isnan(allMapsSPHENTRUE{rep}{i}.fitness(:)))./numel(allMapsSPHENTRUE{rep}{i}.fitness);
end
end
%%
clear fig; close all;
fig(1) = figure(1); hold off;
semilogx(mean(mape_medians,2),'LineWidth',2);
hold on;
plot(mean(sail_mediansTrue,2),'LineWidth',2);
plot(mean(SPHEN_mediansTrue,2),'LineWidth',2);
ax = gca;
ax.YAxis.Limits = [0 1.0];
ax.XTick = [64/8 128 1024 pMAPE.nGens];
ax.XTickLabel = ax.XTick*p.nChildren;
xlabel('Precise Evaluations');
grid on;
title('Median True Fitness');
legend('MAP-Elites','SAIL','SPHEN', 'Location','SouthEast');



fig(2) = figure(2); hold off;
semilogx(mean(mape_bins,2),'LineWidth',2);
hold on;
plot(mean(sail_bins,2),'LineWidth',2);
plot(mean(SPHEN_bins,2),'LineWidth',2);
ax = gca;
ax.YAxis.Limits = [0 100]
ax.XTick = [64/8 128 1024 pMAPE.nGens];
ax.XTickLabel = ax.XTick*p.nChildren;
xlabel('Precise Evaluations');

grid on;
title('% Bins Filled');
legend('MAP-Elites','SAIL','SPHEN','Location','SouthEast');

save_figures(fig, '.', [experimentName 'fitANDfill'], 12, [5 4]);

