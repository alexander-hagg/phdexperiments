clear;clc;
DOF = 16;
workpath = ['/home/' getenv('USER') '/archetype/'];
addpath(genpath(workpath));
DOMAIN = 'footprints';  rmpath(genpath([workpath 'domain'])); addpath(genpath([workpath 'domain/' DOMAIN]));
QD = 'grid';            rmpath(genpath([workpath 'QD/grid'])); rmpath(genpath([workpath 'QD/voronoi'])); addpath(genpath([workpath 'QD/' QD]));
SAQD = 'sphen';          rmpath(genpath([workpath 'QD/sail'])); rmpath(genpath(['QD/sphen'])); addpath(genpath([workpath 'QD/' SAQD]));

load('LETTUCE_run_4.mat');
%%
%p = sphenP;
%p.nGens = 4096;
%[predMap, allMaps] = createPredictionMap([surrogateFitnessSPHEN{1},surrogateFeaturesSPHEN{1}{:}],sphenP,d)
%disp('TAKE THIS PREDMAP WITH SAMPLES INSTEAD');
%%

set(0,'DefaultFigureWindowStyle','default')

figHandle = figure(1); hold off;
fig(1) = figHandle;
genes = reshape(mapSPHEN{1}.genes,[],d.dof);
tRes = size(mapSPHEN{1}.fitness,1);
x = 1:tRes; y = 1:tRes;
[X,Y] = ndgrid(x,y); coordinates = [X(:),Y(:)];
coordinates(:,2) = (tRes+1)-coordinates(:,2);
nans = all(isnan(genes)');
genes(nans,:) = []; coordinates(nans,:) = [];
figHandle = showPhenotype(genes,d,figHandle,coordinates);
axisLimits = [-0.5 tRes+0.5 -0.5 tRes+0.5];
axis(axisLimits);
axis equal;
set(gca, 'Visible', 'off')

%%
set(0,'DefaultFigureWindowStyle','default')


fig(2) = figure(2);
viewMap(mapSPHEN{1},d);caxis([0 1]);colormap(parula(10));

%[reducedmap] = analyzeMaps(mapSPHEN{1},d,sphenP,false,3);
%fig(3) = figure(3);
%viewMap(reducedmap,d);caxis([0 1]);colormap(parula(10));

%% Find 9 examples
figHandle = figure(4); hold off;
fig(4) = figHandle;
%genes = reshape(reducedmap.genes,[],d.dof);
%features = reshape(reducedmap.features,[],2);
%fitness = reshape(reducedmap.fitness,[],1);
select = [1,4; ...      %A
          7,7; ...      %B
          1,21; ...    %C
          5,19; ...      %D
          4,31; ...      %E
          19,21; ...    %F
          18,27; ...    %G
          27,32; ...    %H
          29,28]        %I
      
mapLinIndx = sub2ind(mapSPHEN{1}.resolution,select(:,1),select(:,2));
genes = reshape(mapSPHEN{1}.genes,[],d.dof);
features = reshape(mapSPHEN{1}.features,[],2);
fitness = reshape(mapSPHEN{1}.fitness,[],1);

genes = genes(mapLinIndx,:);
features = features(mapLinIndx,:);
fitness = fitness(mapLinIndx);

tRes = ceil(sqrt(size(genes,1)))
x = 1:tRes; y = 1:tRes;
[X,Y] = ndgrid(x,y); coordinates = [X(:),Y(:)];
coordinates(:,2) = (tRes+1)-coordinates(:,2);
nans = all(isnan(genes)');
genes(nans,:) = []; coordinates(nans,:) = []; features(nans,:) = []; fitness(nans) = [];
figHandle = showPhenotype(genes,d,figHandle,coordinates);
axisLimits = [-0.5 tRes+0.5 -0.5 tRes+0.5];
axis(axisLimits);
axis equal;
set(gca, 'Visible', 'off')

%% Find locations of 9 examples in old map
locationmap                                         = createMap(d, sphenP);
[replaced, replacement, features]                   = nicheCompete(genes, fitness, locationmap, d, sphenP, features);
locationmap                                         = updateMap(replaced,replacement,locationmap,0.1:0.1:0.9,genes,features);
fig(5) = figure(5);
viewMap(locationmap,d)

save('finalResult.mat');

%% Run 9 examples through Lettuce (manual, but create the folders here)
d.nGPUs = 9
[TRUEfitness,TRUEfeatures,polyshapes,booleanMap,TRUEunnormalizedFeatures,rawData,TRUEallfeatures] = fitnessLettuce(genes,d)

save('finalResult_withSims.mat');

%% Errors and rank correlation
sqrt(immse(fitness,TRUEfitness))
sqrt(immse(features(:,1),TRUEfeatures(:,1)))
sqrt(immse(features(:,2),TRUEfeatures(:,2)))
[tau,pval] = corr(fitness,TRUEfitness,'Type','Kendall')
[tau,pval] = corr(features(:,1),TRUEfeatures(:,1),'Type','Kendall')
[tau,pval] = corr(features(:,2),TRUEfeatures(:,2),'Type','Kendall')

taufig(1) = figure(1);
subplot(3,1,1);hold off
scatter(1:length(fitness),fitness);hold on;scatter(1:length(fitness),TRUEfitness);
legend('Predicted','True');
subplot(3,1,2);hold off
scatter(1:length(fitness),features(:,1));hold on;scatter(1:length(fitness),TRUEfeatures(:,1));
%legend('Predicted','True');
subplot(3,1,3);hold off
scatter(1:length(fitness),features(:,2));hold on;scatter(1:length(fitness),TRUEfeatures(:,2));
%legend('Predicted','True');

save_figures(taufig, '.', ['lettuceSPHEN_shapes_taufig'], 12, [8 8]);


%% Create ground truth map
clear unnormalizedFlowFeatures e E
polyshapes = getPhenotypeFFD(genes,d.base);
winSize = 1000;
cd('/scratch/ahagg2s/ARCHETYPE/workfolderCLUSTER');
for i=1:9
    cd(int2str(i))
    e = csvread('AllEnstrophies'); ut = csvread('u0'); u = ut(:,1);
    E = movmean(e,min(size(e,1),winSize)); maxU = movmean(u,min(size(e,1),winSize));    
    unnormalizedFlowFeatures(i,:) = [E(end),maxU(end)];
    cd('..')
end

flowFeatures(:,1) = (unnormalizedFlowFeatures(:,1)-d.featureMin(3))./(d.featureMax(3)-d.featureMin(3));
flowFeatures(:,2) = (unnormalizedFlowFeatures(:,2)-d.featureMin(4))./(d.featureMax(4)-d.featureMin(4));
flowFeatures(flowFeatures>1) = 1; flowFeatures(flowFeatures<0) = 0;

% Shape features
[shapeFeatures,unnormalizedShapeFeatures] = categorize(polyshapes,d);
unnormalizedFeatures = [unnormalizedShapeFeatures,unnormalizedFlowFeatures];

allfeatures = [shapeFeatures,flowFeatures];
featuresGT = allfeatures(:,d.selectedFeatures);
fitnessGT = (1./(1+flowFeatures(:,2)))*2-1 ;
genesGT = genes;

%%

[reducedmap] = analyzeMaps(mapSPHEN{1},d,sphenP,false,3);
fig(3) = figure(3);
viewMap(reducedmap,d);caxis([0 0.5]);colormap(parula(10));

fig(6) = figure(6);
redp = sphenP;redp.resolution = 6;
gtmap                                         = createMap(d, redp);
[replaced, replacement, features]             = nicheCompete(genesGT, fitnessGT, gtmap, d, redp, featuresGT);
gtmap                                         = updateMap(replaced,replacement,gtmap,fitnessGT,genesGT,featuresGT);
viewMap(gtmap,d);caxis([0 0.5]);colormap(parula(10));
%%
fff =reshape(reducedmap.features,[],2);
figure(99);
subplot(2,1,1);hold off
scatter(1:length(fff(:,1)),fff(:,1),64,'k','filled');
hold on
scatter(1:length(featuresGT(:,1)),featuresGT(:,1),32,'r','filled');
title('Area');
subplot(2,1,2);hold off;
scatter(1:length(fff(:,2)),fff(:,2),64,'k','filled');
hold on
scatter(1:length(featuresGT(:,2)),featuresGT(:,2),32,'r','filled');
title('Enstrophy');
%%
err = sqrt(immse(fff(:,1),featuresGT(:,1)))
err = sqrt(immse(fff(:,2),featuresGT(:,2)))

%%
save_figures(fig, '.', ['lettuceSPHEN_shapes'], 12, [8 8]);

%% Create GIFs
set(0,'DefaultFigureWindowStyle','default')
nMaps = length(allMaps);
featMultiplier = 50;
filename = 'lettuceShapeEvolution.gif';
close all; clear figHandle;
for i=[1,nMaps]%1:nMaps
    figHandle(i) = figure(i);
    figHandle(i).Renderer = 'painters';
    figHandle(i).Position = [0 0 800 800];
    hold off;delete(findall(gcf,'type','annotation'));
    
    genes = reshape(allMaps{i}.genes,[],d.dof);
    fitness = reshape(allMaps{i}.fitness,[],1);
    tRes = size(allMaps{i}.fitness,1);
    x = 1:tRes; y = 1:tRes;
    [X,Y] = ndgrid(x,y); features = [X(:),Y(:)];
    features(:,2) = 33-features(:,2);
    nans = all(isnan(genes)');
    genes(nans,:) = []; features(nans,:) = []; fitness(nans) = [];
    clrs = fitness.*[0 2 0] + (1-fitness).*[1 0 0]; clrs(clrs>1) = 1;
    figHandle(i) = showPhenotype(genes,d,figHandle(i),features,clrs);
    set(gca, 'Visible', 'off')
    axisLimits = [-0.5 tRes+0.5 -0.5 tRes+0.5];
    axis(axisLimits);
    %axis tight manual % this ensures that getframe() returns a consistent size
    %axis manual
    %axis equal;
    dim = [.8 .65 .3 .3];
    str = ['Iteration ' int2str(i) '/' int2str(64)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
end
drawnow;

%%
for i=1:nMaps
    %% Capture the plot as an image
    frame = getframe(figHandle(i));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',0.3);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.3);
    end
end



