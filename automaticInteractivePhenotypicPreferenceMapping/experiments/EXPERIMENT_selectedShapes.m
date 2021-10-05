% EXPERIMENT_UserInteraction
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Apr 2020; Last revision: 01-May-2020
%
%------------- BEGIN CODE --------------

clear;clc;
%%
DOF = 32;DOMAIN = 'polygons';ALGORITHM = 'grid';LATENTMODEL = 'VAE';
addpath(genpath('/home/alex/IPPM')); rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN])); rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]); rmpath(genpath('latentmodels')); addpath(genpath(['latentmodels/' LATENTMODEL]));
d = domain(DOF);
p = defaultParamSet;
d.resolution = 128;

m = cfgLatentModel('data/workdir',d.resolution);
poemCfg = poemParamSet(p,m);
poemCfg.numInitSamples                = p.numInitSamples;
poemCfg.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);
poemCfg.numIterations = 1;

fname = ['asymmetry'];

%% Experiment
% User Interaction Demonstration
p.numInitSamples = 64;

% Initialize initial solution set
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,polygons] = fitfun(initSamples,d);

[map{1}, config{1}, stats{1}] = poem(initSamples,polygons,fitness,poemCfg,d,2);
close all; clear figs;save([fname 'step1.mat']);
disp('Optimization Step 1 Done');

%%
load([fname 'step1.mat']);
[genes,fitness,features,bins] = extractMap(map{1});

FITNESSFUNCTION = 'userANDasymmetry';  rmpath(genpath('domain/polygons/fitnessFunctions')); addpath(genpath(['domain/polygons/fitnessFunctions/' FITNESSFUNCTION]));
d.userModel = stats{1}.models{end};
selectionIDs = [17,279];

figs(2) = figure(4);
dSel = d;
dSel.phenoDistMult = 1;
showPhenotype(genes(selectionIDs,:),dSel,figs(2));title('Selected Shapes');
axis equal;

figs(1) = figure(5);
cmap = [0 0 0; 0 0 1];
colors = repmat(cmap(1,:),size(genes,1),1);
colors(selectionIDs,:) = repmat(cmap(2,:),numel(selectionIDs),1);
showPhenotype(genes,d,figs(1),bins,colors);


phenotypes = d.getPhenotype(genes);
features = predictFeatures(phenotypes,d.userModel);
d.selectedShapes = features(selectionIDs,:);
d.deselectedShapes = features; d.deselectedShapes(selectionIDs,:) = [];

% Initialize new initial solution set
% Inject selected solutions
newSamples = genes(selectionIDs,:); 

% Create mutated versions of selection
initMutationSigma = 0.1; % % % % zoomfactor ?

sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
newSampleMutations = 2*initMutationSigma*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)-initMutationSigma;
for i=1:length(selectionIDs)
    newSamples = [newSamples; genes(selectionIDs(i),:) + newSampleMutations];
end

% Recalculate fitness and phenotypes
[newSamplesfitness,newSamplespolygons] = fitfun(newSamples,d);


figs(3) = figure(6);
showPhenotype(newSamples,d,figs(3));title('Injected Shapes');
%%
[map{2}, config{2}, stats{2}] = poem(newSamples,newSamplespolygons,newSamplesfitness,poemCfg,d,2);

close all; clear figs; save([fname 'step2.mat']);
disp('Optimization Step 2 Done');

%%

load([fname 'step2.mat']);
[genes,fitness,features,bins] = extractMap(map{2});
figs(4) = figure(7);

distances = pdist2(genes,newSamples(1:numel(selectionIDs),:));
[row,col] = find(distances==0);
cmap = redgreencmap(256+1,'Interpolation','linear'); cmap = flipud(cmap);
fitBins = floor(256*fitness)+1;
colors = cmap(fitBins,:);
colors(row,:) = repmat([0 0 1],numel(row),1);
showPhenotype(genes,d,figs(4),bins,colors);
%%
save_figures(figs,'.','userselectionANDasymmetry',12,[4 4])








