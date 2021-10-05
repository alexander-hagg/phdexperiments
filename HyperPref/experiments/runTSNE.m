clear;clc;
%% Configuration
addpath(genpath(pwd))                           % Set path to all modules
DOF = 16;DOMAIN = 'catmullRom';                 % Degrees of freedom, Catmull-Rom spline domain
axialBoundAdaptation = 0.1;
radialBoundAdaptation = 2;
d = domain(DOF,axialBoundAdaptation,radialBoundAdaptation);                                % Domain configuration
p = defaultParamSet;                            % Base Quality Diversity (QD) configuration (MAP-Elites)
m = cfgLatentModel('data/workdir',d.resolution);% VAE configuration
pm = poemParamSet(p,m);                    % Configure POEM ("Phenotypic niching based Optimization by Evolving a Manifold")
pm.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);  % Anonymous function ptr to phenotypic categorization function (= VAE)
hypID = 1;

%% Test
clear matchedParams;
fig = figure(1);ax=gca;
matchedParams(1,:) = [1 1 1 1 1 1 1 1 zeros(1,8)];
matchedParams(2,:) = [1 1 0 1 1 1 1 1 zeros(1,8)];
matchedParams(3,:) = [1 sin(0.25*pi) 1 sin(0.25*pi) 1 sin(0.25*pi) 1 sin(0.25*pi) zeros(1,8)];
matchedParams(4,:) = [1 0.3 1 0.3 1 0.3 1 0.3 zeros(1,8)];
showPhenotype(matchedParams,d,1.2,ax);
axis([-2 2 -2 2]); axis equal;



%% Generate data
phenotypeRes = 64;
numShapes = 13;
scaling = [0.5 1.5];
rotation = [-pi pi];
%load('matchedParams.mat');

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
fig(1) = figure(1); hold off; ax = gca;
showPhenotype(genomes,d,1.2,ax,[],clrs);

fig(2) = figure(2); hold off; ax = gca;
gCoords = distMultiplier*mapminmax(simspaceCoordinatesG',0,1)';
showPhenotype(genomes,d,1.2,ax,gCoords,clrs);

fig(3) = figure(3); hold off; ax = gca;
pCoords= distMultiplier*mapminmax(simspaceCoordinatesP',0,1)';
showPhenotype(genomes,d,1.2,ax,pCoords,clrs);

%%
fig(1).Renderer='Painters';
fig(2).Renderer='Painters';
fig(3).Renderer='Painters';
save_figures(fig, '.', ['tsneGvsP'], 12, [6 6]);



