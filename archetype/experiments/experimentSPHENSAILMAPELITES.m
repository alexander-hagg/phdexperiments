clear;clc;
DOF = 16;
DOMAIN = 'footprints';
QD = 'grid';
SAQD = 'sphen'; 

addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath(genpath('QD/grid')); rmpath(genpath('QD/voronoi')); addpath(genpath(['QD/' QD]));
rmpath(genpath('QD/sail')); rmpath(genpath('QD/sphen')); addpath(genpath(['QD/' SAQD]));

d = domain(DOF);
p = defaultParamSet;

p.infill = infillParamSet;
p.infill.nTotalSamples = 1024;
p.infill.nAdditionalSamples = p.nChildren;
p.numInitSamples = 16;
p.display.illu = false;


pMAPE               = p;                pMAPE.nGens = 2^12;

pSAIL               = p;                pSAIL.nGens = 2^10; pSAIL.nChildren = 32;   
numSAILrounds       = p.infill.nTotalSamples/p.infill.nAdditionalSamples;

restrictedPESAIL    = (pMAPE.nGens*p.nChildren-p.infill.nTotalSamples)/numSAILrounds; % 1008 PE per SAIL iteration left
pSAILrestricted     = p;                pSAILrestricted.nGens = restrictedPESAIL/p.nChildren; % 48 generations per SAIL iteration left
sailTotalSamples = numSAILrounds*pSAIL.nGens*pSAIL.nChildren+p.infill.nTotalSamples;

pSPHEN              = p;                pSPHEN.nGens = 2^10; pSPHEN.nChildren = 32; 
disp('Set resolution to double to see whether we can get rid of the holes');
pSPHEN.resolution = 32;     

d.fitfun = d.fitfunPointSymm; % Multimodal function: point symmetry (center); for testing purposes
experimentName = 'SPHENSAILMAPELITES_pointSymmetry'
numReplicates = 5;
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');

%% ----------------------------------------------------------------------------------
set(0,'DefaultFigureWindowStyle','docked')

for rep=1 : numReplicates
    disp(['Replicate ' int2str(rep) '/' int2str(numReplicates)]);
    % Create initial sample set
    sobPoint = (rep-1)*p.numInitSamples + 1;
    initSet{rep}.samples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
    [initSet{rep}.fitness,initSet{rep}.features] = d.fitfun(initSet{rep}.samples);
    
    %% MAP-Elites
    initmap{rep}                                            = createMap(d, pMAPE);
    [replaced, replacement, features]                       = nicheCompete(initSet{rep}.samples, initSet{rep}.fitness, initmap{rep}, d, pMAPE, initSet{rep}.features);
    initmap{rep}                                            = updateMap(replaced,replacement,initmap{rep},initSet{rep}.fitness,initSet{rep}.samples,initSet{rep}.features);
    [mapMAPELITES{rep}, ~, ~, allMaps{rep}]                 = illuminate(initmap{rep},d.fitfun,pMAPE,d);    

    %% SAIL (Surrogate-Assisted Illumination)
    SAQD = 'sail';
    rmpath(genpath('QD/sail')); rmpath(genpath('QD/SPHEN')); addpath(genpath(['QD/' SAQD]));
    [mapSAIL{rep},surrogate{rep},allMapsSAIL{rep}] = sail(initSet{rep},pSAIL,d);

    %% SAIL (Surrogate-Assisted Illumination) with limited Feature Evaluations
    SAQD = 'sail';
    rmpath(genpath('QD/sail')); rmpath(genpath('QD/SPHEN')); addpath(genpath(['QD/' SAQD]));
    [mapSAIL2{rep},surrogate2{rep},allMapsSAIL2{rep}] = sail(initSet{rep},pSAILrestricted,d);
    
    %% SPHEN (Surrogate-Assisted Phenotypic Niching
    SAQD = 'sphen';
    rmpath(genpath('QD/sail')); rmpath(genpath('QD/SPHEN')); addpath(genpath(['QD/' SAQD]));
    pSPHEN.infill.acqFcn               = 'PureUCB'; %FeatureUnCertainty FeatureCertainty PureUCB
    [mapSPHEN{rep},surrogateFitnessSPHEN{rep},surrogateFeaturesSPHEN{rep},allMapsSPHEN{rep},trueFilledSPHEN(rep)] = sphen(initSet{rep},pSPHEN,d,1);
    
    % Save intermediate experiment state
    save([experimentName '_8.mat']);
    
end

%% Analysis

for rep=1:numReplicates
    [truemap{1,rep},~,              filled(1,rep),medianFitness(1,rep)] = analyzeMaps(mapMAPELITES{rep},d,p);
    [truemap{2,rep},errors(2,rep,:),filled(2,rep),medianFitness(2,rep)] = analyzeMaps(mapSAIL{rep},d,pSAIL,true);
    [truemap{3,rep},errors(3,rep,:),filled(3,rep),medianFitness(3,rep)] = analyzeMaps(mapSAIL2{rep},d,pSAILrestricted,true);
    [truemap{4,rep},errors(4,rep,:),filled(4,rep),medianFitness(4,rep)] = analyzeMaps(mapSPHEN{rep},d,pSPHEN,true);
    [truemap{5,rep},errors(5,rep,:),filled(5,rep),medianFitness(5,rep)] = analyzeMaps(mapSPHEN{rep},d,pSPHEN,true,pSPHEN.resolution/2);

    for a=1:length(allMaps{rep})
        disp(a)
        [alltruemapMAPELITES{a,rep},~,allfilledMAPELITES(a,rep),allmedianFitnessMAPELITES(a,rep)] = analyzeMaps(allMaps{rep}{a},d,p);        
    end
    for a=1:length(allMapsSPHEN{rep})
        disp(a)
        [alltruemapSAIL{a,rep},allerrorsSAIL(a,rep,:),allfilledSAIL(a,rep),allmedianFitnessSAIL(a,rep)] = analyzeMaps(allMapsSAIL{rep}{a},d,pSAIL,true);
        [alltruemapSAIL2{a,rep},allerrorsSAIL2(a,rep,:),allfilledSAIL2(a,rep),allmedianFitnessSAIL2(a,rep)] = analyzeMaps(allMapsSAIL2{rep}{a},d,pSAILrestricted,true);
        [alltruemapSPHEN{a,rep},allerrorsSPHEN(a,rep,:),allfilledSPHEN(a,rep),allmedianFitnessSPHEN(a,rep)] = analyzeMaps(allMapsSPHEN{rep}{a},d,pSPHEN,true);
        [alltruemapSPHEN2{a,rep},allerrorsSPHEN2(a,rep,:),allfilledSPHEN2(a,rep),allmedianFitnessSPHEN2(a,rep)] = analyzeMaps(allMapsSPHEN{rep}{a},d,pSPHEN,true,pSPHEN.resolution/2);
    end
end

%% Plot prediction maps

set(0,'DefaultFigureWindowStyle','default')


algs = {'MAP-Elites','SAIL', 'SAIL restricted', 'SPHEN training', 'SPHEN prediction'};
clear fig; close all;
for i=1:5
    fig(i) = figure(i);ax = gca;
    viewMap(truemap{i,rep},d,'fitness',ax);
    caxis([0.8 1.0]);title(algs{i});
end
save_figures(fig, '.', [experimentName 'QDmapcomparison'], 12, [5 4]);


%% Plot using Precise Fitness Evaluations
clear fig; close all;

fig(1) = figure(1);hold off;
plot([mean(allfilledMAPELITES,2);mean(filled(1,:))]);
hold on;
plot([mean(allfilledSAIL,2);mean(filled(2,:))]);
plot([mean(allfilledSAIL2,2);mean(filled(3,:))]);
plot([mean(allfilledSPHEN,2);mean(filled(4,:))]);
plot([mean(allfilledSPHEN2,2);mean(filled(5,:))]);
legend('MAP-Elites','SAIL','SAIL restricted','SPHEN training','SPHEN prediction','Location','SouthEast');
title('Filled %');grid on;
ax = gca; ax.XScale = 'log';
ax.XTick = [1 16 64 4096]
ax.XTickLabel = ax.XTick*p.nChildren
xlabel('Fitness Evaluations');

fig(2) = figure(2);hold off;
plot([mean(allmedianFitnessMAPELITES,2);mean(medianFitness(1,:))]);
hold on;
plot([mean(allmedianFitnessSAIL,2);mean(medianFitness(2,:))]);
plot([mean(allmedianFitnessSAIL2,2);mean(medianFitness(3,:))]);
plot([mean(allmedianFitnessSPHEN,2);mean(medianFitness(4,:))]);
plot([mean(allmedianFitnessSPHEN2,2);mean(medianFitness(5,:))]);
legend('MAP-Elites','SAIL','SAIL restricted','SPHEN training','SPHEN prediction','Location','SouthEast');
title('Fitness');grid on;
ax = gca; ax.XScale = 'log';
ax.XTick = [1 16 64 4096]
ax.XTickLabel = ax.XTick*p.nChildren
xlabel('Fitness Evaluations');

save_figures(fig, '.', [experimentName 'PE'], 12, [5 4]);
%% Plot using Precise Fitness/Feature Evaluations
clear fig; close all;
set(0,'DefaultFigureWindowStyle','default')

fig(1) = figure(1);hold off;
plot([mean(allfilledMAPELITES,2);mean(filled(1,:))]);
hold on;
sailGens = sailTotalSamples/p.nChildren;
plot(sailGens/64:sailGens/64:sailGens,[mean(allfilledSAIL,2);mean(filled(2,:))]);
initPtSAILrestricted = (pSAILrestricted.nGens*pSAILrestricted.nChildren+16)./p.nChildren;
%plot(initPtSAILrestricted:64:4096,mean(allfilledSAIL2,2));
plot(65:64:4096,mean(allfilledSAIL2,2));
plot([mean(allfilledSPHEN,2);mean(filled(4,:))]);
plot([mean(allfilledSPHEN2,2);mean(filled(5,:))]);
legend('MAP-Elites','SAIL','SAIL restricted','SPHEN training','SPHEN prediction','Location','SouthEast');
title('Filled %');grid on;
ax = gca; ax.XScale = 'log';
ax.XTick = [1 16 64 4096]
ax.XTickLabel = ax.XTick*p.nChildren
xlabel('Fitness or Feature Evaluations');


fig(2) = figure(2);hold off;
plot([mean(allmedianFitnessMAPELITES,2);mean(medianFitness(1,:))]);
hold on;
plot(sailGens/64:sailGens/64:sailGens,[mean(allmedianFitnessSAIL,2);mean(medianFitness(2,:))]);
plot(2:64:4096,[mean(allmedianFitnessSAIL2,2);mean(medianFitness(3,:))]);
plot([mean(allmedianFitnessSPHEN,2);mean(medianFitness(4,:))]);
plot([mean(allmedianFitnessSPHEN2,2);mean(medianFitness(5,:))]);
legend('MAP-Elites','SAIL','SAIL restricted','SPHEN training','SPHEN prediction','Location','SouthEast');
title('Fitness');grid on;
ax = gca; ax.XScale = 'log';
ax.XTick = [1 16 64 4096]
ax.XTickLabel = ax.XTick*p.nChildren
xlabel('Fitness or Feature Evaluations');

save_figures(fig, '.', [experimentName 'PFE'], 12, [5 4]);

%%
clear fig; close all;
fig(1) = figure(1);
titles = {'Fitness (Symmetry)','Area','Circumference'};
meanallerrorsSAIL = squeeze(mean(allerrorsSAIL,2));
meanallerrorsSAIL2 = squeeze(mean(allerrorsSAIL2,2));
meanallerrorsSPHEN = squeeze(mean(allerrorsSPHEN,2));
%meanallerrorsSPHEN2 = squeeze(mean(allerrorsSPHEN2,2));
for i=1:3
    subplot(3,1,i);
    %for j=2:3
        plot(meanallerrorsSAIL(:,i)); 
        hold on; 
        plot(meanallerrorsSAIL2(:,i)); 
        plot(meanallerrorsSPHEN(:,i)); 
        %plot(meanallerrorsSPHEN2(:,i)); 
        
    %end
    ax = gca; ax.YAxis.Limits = [0 0.12];
    title(titles{i});grid on;
    if i==1; legend('SAIL','SAIL reduced', 'SPHEN', 'Location','NorthEast'); end    
end


save_figures(fig, '.', [experimentName 'errors'], 12, [5 6]);

%% Show difference between SPHEN 32x32 and 16x16
SPHEN_original         = mapSPHEN{rep};
SPHEN_originaltrue     = analyzeMaps(mapSPHEN{rep},d,pSPHEN,true);
SPHEN_reduced          = analyzeMaps(mapSPHEN{rep},d,pSPHEN,false,pSPHEN.resolution/2);
SPHEN_reducedtrue      = analyzeMaps(mapSPHEN{rep},d,pSPHEN,true,pSPHEN.resolution/2);

fig(1) = figure(1);hold off;
viewMap(SPHEN_original,d);
caxis([0.8 1]);title('32x32');

fig(2) = figure(2);hold off;
viewMap(SPHEN_originaltrue,d);
caxis([0.8 1]);title('32x32 True');

fig(3) = figure(3);hold off;
viewMap(SPHEN_reduced,d);
caxis([0.8 1]); title('16x16');

fig(4) = figure(4);hold off;
viewMap(SPHEN_reducedtrue,d);
caxis([0.8 1]); title('16x16 True');

save_figures(fig, '.', [experimentName 'resolutioncomparison'], 12, [5 4]);



%% Significance tests
fig(1) = figure(1);hold off;
plot([(allfilledMAPELITES);(filled(1,:))]);
hold on;
sailGens = sailTotalSamples/p.nChildren;
plot(sailGens/64:sailGens/64:sailGens,[(allfilledSAIL);(filled(2,:))]);
initPtSAILrestricted = (pSAILrestricted.nGens*pSAILrestricted.nChildren+16)./p.nChildren;
%plot(initPtSAILrestricted:64:4096,mean(allfilledSAIL2,2));
plot(65:64:4096,(allfilledSAIL2));
%plot([(allfilledSPHEN);(filled(4,:))]);
plot([(allfilledSPHEN2);(filled(5,:))]);
legend('MAP-Elites','SAIL','SAIL restricted','SPHEN','Location','SouthEast');
title('Filled %');grid on;
ax = gca; ax.XScale = 'log';
ax.XTick = [1 16 64 4096]
ax.XTickLabel = ax.XTick*p.nChildren
xlabel('Fitness or Feature Evaluations');

fig(2) = figure(2);hold off;
plot([(allmedianFitnessMAPELITES);(medianFitness(1,:))]);
hold on;
plot(sailGens/64:sailGens/64:sailGens,[(allmedianFitnessSAIL);(medianFitness(2,:))]);
plot(2:64:4096,[(allmedianFitnessSAIL2);(medianFitness(3,:))]);
%plot([(allmedianFitnessSPHEN);(medianFitness(4,:))]);
plot([(allmedianFitnessSPHEN2);(medianFitness(5,:))]);
legend('MAP-Elites','SAIL','SAIL restricted','SPHEN','Location','SouthEast');
title('Fitness');grid on;
ax = gca; ax.XScale = 'log';
ax.XTick = [1 16 64 4096]
ax.XTickLabel = ax.XTick*p.nChildren
xlabel('Fitness or Feature Evaluations');

save_figures(fig, '.', [experimentName 'PFE'], 12, [5 4]);

%% Find PFE where filled/fitness values are the same (based on min filled/fitness value of all tested runs)
min1 = min([allfilledSPHEN2(end,:),allfilledMAPELITES(end,:)]);
min2 = min([allfilledSPHEN2(end,:),allfilledSAIL(end,:)]);
min3 = min([allfilledSPHEN2(end,:),allfilledSAIL2(end,:)]);


clear mfs
for i=1:5
    mfs(1,1,i) = min(find(allfilledMAPELITES(:,i) >= min1));
    mfs(1,2,i) = min(find(allfilledSPHEN2(:,i) >= min1));

    mfs(2,1,i) = min(find(allfilledSAIL(:,i) >= min2));
    mfs(2,2,i) = min(find(allfilledSPHEN2(:,i) >= min2));

    mfs(3,1,i) = min(find(allfilledSAIL2(:,i) >= min3));
    mfs(3,2,i) = min(find(allfilledSPHEN2(:,i) >= min3));
end

% MAP-Elites 4096 * p.nChildren = [1:4096]*16
% SAIL sailGens/64:sailGens/64:sailGens
%64+64*63

mfs(1,1,:) = mfs(1,1,:)*p.nChildren;
mfs(2,1,:) = mfs(2,1,:)*sailGens/64;
mfs(3,1,:) = mfs(3,1,:)*4096*16/63;
mfs(1,2,:) = mfs(1,2,:)*p.nChildren;
mfs(2,2,:) = mfs(2,2,:)*p.nChildren;
mfs(3,2,:) = mfs(3,2,:)*p.nChildren;

[~,pp] = ttest2(mfs(1,2,:),mfs(1,1,:))
[~,pp] = ttest2(mfs(2,2,:),mfs(2,1,:))
[~,pp] = ttest2(mfs(3,2,:),mfs(3,1,:))
%%

min1 = min([allmedianFitnessSPHEN2(end,:),allmedianFitnessMAPELITES(end,:)]);
min2 = min([allmedianFitnessSPHEN2(end,:),allmedianFitnessSAIL(end,:)]);
min3 = min([allmedianFitnessSPHEN2(end,:),allmedianFitnessSAIL2(end,:)]);


clear mfs
for i=1:5
    mfs(1,1,i) = min(find(allmedianFitnessMAPELITES(:,i) >= min1));
    mfs(1,2,i) = min(find(allmedianFitnessSPHEN2(:,i) >= min1));

    mfs(2,1,i) = min(find(allmedianFitnessSAIL(:,i) >= min2));
    mfs(2,2,i) = min(find(allmedianFitnessSPHEN2(:,i) >= min2));

    mfs(3,1,i) = min(find(allmedianFitnessSAIL2(:,i) >= min3));
    mfs(3,2,i) = min(find(allmedianFitnessSPHEN2(:,i) >= min3));
end

% MAP-Elites 4096 * p.nChildren = [1:4096]*16
% SAIL sailGens/64:sailGens/64:sailGens
%64+64*63

mfs(1,1,:) = mfs(1,1,:)*p.nChildren;
mfs(2,1,:) = mfs(2,1,:)*sailGens/64;
mfs(3,1,:) = mfs(3,1,:)*4096*16/63;
mfs(1,2,:) = mfs(1,2,:)*p.nChildren;
mfs(2,2,:) = mfs(2,2,:)*p.nChildren;
mfs(3,2,:) = mfs(3,2,:)*p.nChildren;

[~,pp] = ttest2(mfs(1,2,:),mfs(1,1,:))
[~,pp] = ttest2(mfs(2,2,:),mfs(2,1,:))
[~,pp] = ttest2(mfs(3,2,:),mfs(3,1,:))
