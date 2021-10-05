%% TODO
%  - does point symmetry really work when we get NANs as well?

%% Configuration
clear;clc;
DOF = 16;
DOMAIN = 'npoly_ffd';
ALGORITHM = 'grid';

addpath(genpath('.'));
rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
d = domain(DOF);

rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);
p = defaultParamSet(4);


%% I) Multiobjective optimization with NSGA-II
numObjectives = evaluate_objective;
popsize = p.featureResolution(1)*p.featureResolution(2);
numgenerations = p.nGens;
tic
returnValuesNSGA = nsga_2(popsize,numgenerations, numObjectives, d.dof, d.ranges(1), d.ranges(2), d, p);
toc

genomes = returnValuesNSGA(:,1:d.dof);
fitness = 1-returnValuesNSGA(:,d.dof+1); % Turn back into maximized fitness
features = [1-returnValuesNSGA(:,d.dof+2),returnValuesNSGA(:,d.dof+3)];

phenotypes = getPhenotype(genomes,d);
NSGAmap = createMap(d, p);
[replaced, replacement, features]       = nicheCompete(genomes, -fitness, phenotypes, NSGAmap, d, p,features);
NSGAmap                                 = updateMap(replaced,replacement,NSGAmap,fitness,genomes,[],features,[]);

%% II) Illumination with QD
disp(['Illumination']);
p.numInitSamples = p.featureResolution(1)*p.featureResolution(2);
p.nChildren = p.featureResolution(1)*p.featureResolution(2);
[QDmap,fitnessFunction] = initialize([],d,p,false);
tic
QDmap = illuminate(QDmap,fitnessFunction,p,d);
toc

save('data/MOOvsQD/MOOvsQD.mat','returnValuesNSGA','NSGAmap','QDmap','d','p');

%% Show NSGA-II and MAP-Elites results in feature maps
d.featureLabels = d.featureLabels(d.selectedFeatures);


figure(1);hold off;figHandle = gca;
[figHandle, imageHandle, cHandle] = viewMap(NSGAmap,d,figHandle);caxis([0 1]);
title('MOO: NSGA-II');
cHandle.Label.String = 'point symmetry';

%
figure(2);hold off;figHandle = gca;
[figHandle, imageHandle, cHandle] = viewMap(QDmap,d,figHandle);caxis([0 1]);
title('QD: MAP-Elites');
cHandle.Label.String = 'point symmetry';

%% Comparison

% Compare genetic diversity

genesMAPE = reshape(QDmap.genes,[],d.dof);
genesMAPE(any(isnan(genesMAPE)'),:) = [];
fitMAPE = reshape(QDmap.fitness,[],1);
fitMAPE = fitMAPE(~isnan(fitMAPE));

X = [returnValuesNSGA(:,1:d.dof);genesMAPE];
mappedFit = [returnValuesNSGA(:,d.dof+1);fitMAPE];

[coeffs,score,latent,tsquared,explained,mu] = pca(X);
no_dims = 6;
score = score(:,1:no_dims);
explained = explained(1:no_dims);
%%
pcaDims = 3:4
mappedX = score(:,pcaDims);
disp(['PCA: % variance explained: ' num2str(sum(explained))]);

scaledFit = mappedFit;
fitThresh = 0.7;
scaledFit(mappedFit<fitThresh) = 8;
scaledFit(mappedFit>=fitThresh) = 32;

figure(3);hold off;
scatter(mappedX(1:size(returnValuesNSGA,1),1),mappedX(1:size(returnValuesNSGA,1),2),scaledFit(1:size(returnValuesNSGA,1)),[1 0 0],'filled');
hold on;
scatter(mappedX(size(returnValuesNSGA,1)+1:end,1),mappedX(size(returnValuesNSGA,1)+1:end,2),scaledFit(size(returnValuesNSGA,1)+1:end),[0 0 1],'filled');
xlabel(['PCA dim ' int2str(pcaDims(1))]);
ylabel(['PCA dim ' int2str(pcaDims(2))]);
title(['Genetic Diversity (' alg ' projection)']); legend('NSGA-II','MAP-Elites');

% Compare diversity in feature dimensions
disp(['NSGA-II % map filled: ' num2str(100*sum(~isnan(NSGAmap.fitness(:)))./numel(NSGAmap.fitness))]);
disp(['MAP-Elites % map filled: ' num2str(100*sum(~isnan(QDmap.fitness(:)))./numel(QDmap.fitness))]);

%% PCA on both sets separately
genesNSGA = returnValuesNSGA(:,1:d.dof);
alg='PCA';
[genesMAPE_PCA,score,latent,tsquared,explained,mu] = pca(genesMAPE);

disp(['PCA: % variance explained: ' num2str(sum(explained(1:9))) ' with 9 dims']);

[genesNSGA_PCA,score,latent,tsquared,explained,mu] = pca(genesNSGA);

disp(['PCA: % variance explained: ' num2str(sum(explained(1:2))) ' with 2 dims']);




%% Compare performance in NSGA-II space (compare Pareto fronts)
V = d.dof; M = 3;
frontNSGA = non_domination_sort_mod(returnValuesNSGA, M, V);

fitnessQD = reshape(QDmap.fitness,[],1);fitnessQD = fitnessQD(~isnan(fitnessQD));
featuresQD = reshape(QDmap.features,[],2);featuresQD = featuresQD(all(~isnan(featuresQD)'),:);
returnValuesMAPE = [genesMAPE,fitnessQD,featuresQD];

figure(4);hold off;
scatter3(1-returnValuesNSGA(:,d.dof + 1),1-returnValuesNSGA(:,d.dof + 2),returnValuesNSGA(:,d.dof + 3),16,[1 0 0],'filled');
hold on;
scatter3(returnValuesMAPE(:,d.dof + 1),returnValuesMAPE(:,d.dof + 2),returnValuesMAPE(:,d.dof + 3),16,[0 0 1],'filled');
legend('NSGA-II','MAP-Elites');
xlabel('point symmetry');
ylabel('surface area');
zlabel('perimeter length');
title('Objective Space (1 common objective, 1 feature)');grid on;
eps = 0.1;
axis([-eps 1+eps -eps 1+eps -eps 1+eps]);
view(90,270)


figure(5);hold off;
scatter3(1-returnValuesNSGA(:,d.dof + 1),1-returnValuesNSGA(:,d.dof + 2),returnValuesNSGA(:,d.dof + 3),16,[1 0 0],'filled');
hold on;
scatter3(returnValuesMAPE(:,d.dof + 1),returnValuesMAPE(:,d.dof + 2),returnValuesMAPE(:,d.dof + 3),16,[0 0 1],'filled');
legend('NSGA-II','MAP-Elites');
xlabel('point symmetry');
ylabel('surface area');
zlabel('perimeter length');
title('Objective Space (2 features)');grid on;
eps = 0.1;
axis([-eps 1+eps -eps 1+eps -eps 1+eps]);
view(90,0)
%% Compare optimality in common objective space

group = [    ones(size(returnValuesNSGA(:,d.dof+1)',2),1);
    2 * ones(size(fitnessQD,1),1)];
figure(6);
boxplot([1-returnValuesNSGA(:,d.dof+1); fitnessQD],group)
set(gca,'XTickLabel',{'NSGA-II','MAP-Elites'})
title('Performance in common objective');

%% Compare performance in common feature cells
commonCells = ~isnan(NSGAmap.fitness) & ~isnan(QDmap.fitness);

QDbetter = NSGAmap.fitness(commonCells)<QDmap.fitness(commonCells);
disp(['MAP-Elites performs better than NSGA-II in ' num2str(100*sum(QDbetter)/numel(QDbetter)) '% of cells they have in common (' int2str(sum(QDbetter)) ' cells)']);

%% Parameter plot
genesMAPE = reshape(QDmap.genes,[],d.dof);
genesMAPE(any(isnan(genesMAPE)'),:) = [];

figure(7);hold off;
plot(returnValuesNSGA(:,1:d.dof)','Color',[1 0 0]);
hold on;viewMap(~isnan(NSGAmap.fitness),d)
plot(genesMAPE','Color',[0 0 1]);

%% Clustering
genesMAPE = reshape(QDmap.genes,[],d.dof);
genesMAPE(any(isnan(genesMAPE)'),:) = [];

genesNSGA = returnValuesNSGA(:,1:d.dof);

allGenes = [genesNSGA;genesMAPE];


no_dims = 2;initial_dims=d.dof;perplexity=50;theta = 0.5;
simspaceCoordinates = fast_tsne(allGenes, no_dims, initial_dims, perplexity, theta);alg='tSNE';

%numClusters = 3;
%[labels,cen] = kmedoids(mappedX,numClusters);
coreneighbours      = max(2 * no_dims,3); %Rule of thumb
[~,t_distances]     = knnsearch(simspaceCoordinates,simspaceCoordinates,'K',coreneighbours+1);
t_distances(:,1)    = [];
t_distances         = sort(t_distances(:));
[maxVal ,~]         = getElbow(t_distances);
epsilon             = maxVal;
[~,labels,cen]      = dbscan(simspaceCoordinates', epsilon, coreneighbours);
%%
labels = labels  + 1;
uniqLabels = unique(labels);
colormap = hsv(max(uniqLabels));
figure(8);hold off;
scatter(simspaceCoordinates(:,1),simspaceCoordinates(:,2),64,colormap(labels,:),'filled');
hold on;
scatter(simspaceCoordinates(size(genesNSGA,1)+1:end,1),simspaceCoordinates(size(genesNSGA,1)+1:end,2),64,[1 1 1],'x');


%% Compare maps
figure(9); ax = gca;
occupanceMap.features = [];
occupanceMap.fitness = zeros(size(QDmap.fitness));
occupanceMap.fitness(~isnan(QDmap.fitness)) = 0.3;
occupanceMap.fitness(~isnan(NSGAmap.fitness)) = 1.0;
occupanceMap.fitness(~isnan(NSGAmap.fitness) & ~isnan(QDmap.fitness)) = 0.6;

[~,~,cbHandle] = viewMap(occupanceMap,d,ax);
cbHandle.Ticks = [0 0.3 0.6 1];
clear colormap;
%colormap(parula(4));
colormap([0 0 0; 0.3 0.3 0.3; 0.6 0.6 0.6;1 1 1]);
cbHandle.TickLabels = {'No solutions', 'MAP-Elites','Both','NSGA-II'}


%% Stats
fitNSGA = NSGAmap.fitness(~isnan(NSGAmap.fitness));
fitMAPE = QDmap.fitness(~isnan(QDmap.fitness));

disp('stddev common objective');
std(fitNSGA)
std(fitMAPE)

disp('Total Cells Filled');
numel(fitNSGA)
numel(fitMAPE)

disp('stddev genomes');
mean(std(genesNSGA))
mean(std(genesMAPE))


%% Create image dataset
genesNSGA = returnValuesNSGA(:,1:d.dof);
genesMAPE = reshape(QDmap.genes,[],d.dof);
genesMAPE(any(isnan(genesMAPE)'),:) = [];

resolution = 16;
mkdir([imgDir '/' 'mape' '/' int2str(resolution)]);
phenotypesMAPE = getPhenotype(genesMAPE,d);
for i=1:length(phenotypesMAPE)
    pixelCoordinates = ceil(phenotypesMAPE{i}.Vertices*resolution/2)+resolution/2;
    pixelCoordinates(all(isnan(pixelCoordinates)'),:) = [];
    bw = poly2mask(pixelCoordinates(:,1),pixelCoordinates(:,2),resolution,resolution);
    imwrite(bw,[imgDir '/' 'mape' '/' int2str(resolution) '/' int2str(i) '.png'],'png');
    masksMAPE{i} = bw;
end


mkdir([imgDir '/' 'nsga' '/' int2str(resolution)]);
phenotypesNSGA = getPhenotype(genesNSGA,d);
for i=1:length(phenotypesNSGA)
    pixelCoordinates = ceil(phenotypesNSGA{i}.Vertices*resolution/2)+resolution/2;
    pixelCoordinates(all(isnan(pixelCoordinates)'),:) = [];
    bw = poly2mask(pixelCoordinates(:,1),pixelCoordinates(:,2),resolution,resolution);
    imwrite(bw,[imgDir '/' 'nsga' '/' int2str(resolution) '/' int2str(i) '.png'],'png');
    masksMAPE{i} = bw;
end
