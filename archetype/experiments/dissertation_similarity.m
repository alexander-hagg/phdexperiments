clear;clc;
DOF = 16;DOMAIN = 'footprints';ALGORITHM = 'grid';
addpath(genpath('.'));
rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);
rmpath(genpath('generators/VAE'));
d = domain(DOF);
p = defaultParamSet(4);

%% Generate data
phenotypeRes = 64;
numShapes = 13;
scaling = [0.5 1.5];
rotation = [-pi pi];
load('matchedParams.mat');

scaleArray = [scaling(1):(scaling(2)-scaling(1))/(numShapes-1):scaling(2)]';
rotationArray = [rotation(1):(rotation(2)-rotation(1))/(numShapes-1):rotation(2)]';

genomes = []; iter = 1;
for i=1:size(matchedParams,1)
    for j=1:length(scaleArray)
        for k=1:length(rotationArray)
            genomes(iter,:) = [scaleArray(j).*matchedParams(i,1:DOF/2) rotationArray(k) + matchedParams(i,DOF/2+1:end)];
            iter = iter + 1;
        end
    end
end


%% Train genetic similarity space model 
rng default % for fair comparison
opts = statset('MaxIter',2500);
simspaceCoordinatesG = tsne(genomes,'Perplexity',50,'Options',opts);

% Train phenotypic similarity space model 
polyshapes = getPhenotypeFFD(genomes,d.base);
[flatbooleanMap,booleanMap] = getPhenotypeBoolean(polyshapes);
simspaceCoordinatesP = tsne(flatbooleanMap,'Perplexity',50,'Options',opts);

disp('t-SNE done');
%%
clrMat = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0];
clrs = [];
for c=1:size(matchedParams,1)
    clrs = [clrs; repmat(clrMat(c,:),size(genomes,1)./4,1)];
end
distMultiplier = 50;
fig(1) = figure(1); hold off;
figHandle = showPhenotype(genomes,d,fig(1),[],clrs);


fig(2) = figure(2); hold off;
gCoords= distMultiplier*mapminmax(simspaceCoordinatesG',0,1)';
figHandle = showPhenotype(genomes,d,fig(2),gCoords,clrs);

fig(3) = figure(3); hold off;
pCoords= distMultiplier*mapminmax(simspaceCoordinatesP',0,1)';
figHandle = showPhenotype(genomes,d,fig(3),pCoords,clrs);

%%
fig(1).Renderer='Painters';
fig(2).Renderer='Painters';
fig(3).Renderer='Painters';
save_figures(fig, '.', ['tsneGvsP'], 12, [6 6]);


