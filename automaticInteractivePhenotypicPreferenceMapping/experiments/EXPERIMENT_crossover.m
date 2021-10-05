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
DOF = 16;DOMAIN = 'polygons';ALGORITHM = 'grid';LATENTMODEL = 'VAE';
addpath(genpath('/home/alex/IPPM'));
rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);
rmpath(genpath('latentmodels')); addpath(genpath(['latentmodels/' LATENTMODEL]));
fname = ['userInteraction' '.mat'];
d = domain(DOF);
p = defaultParamSet;

d.resolution = 32;

m = cfgLatentModel('data/workdir',d.resolution);
poemCfg = poemParamSet(p,m);
poemCfg.numInitSamples                = p.numInitSamples;
poemCfg.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);
poemCfg.numIterations = 2;

%% Experiment
% User Interaction Demonstration
p.numInitSamples = 64;

% Initialize initial solution set for MAP-Elites
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
initSamples(end-2,:) = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0];
initSamples(end-1,:) = [1 0.3 0.3 1 0.3 0.3 1 0.3 0 0 0 0 0 0 0 0];
initSamples(end,:)   = [1 1 0.4 0.4 1 1 0.4 0.4 -0 0 0 0 0 0 0 0];
[fitness,~,polygons,rawfitness] = fitfun(initSamples,d);
%showPhenotype(initSamples,d)


poemCfg.map.crossover = true;
[map{1}, config{1}, stats{1}] = poem(initSamples,polygons,fitness,poemCfg,d,2);

poemCfg.map.crossover = false;
[map{2}, config{2}, stats{2}] = poem(initSamples,polygons,fitness,poemCfg,d,3);
save(fname);

disp('Optimization Done');

%%
disp('...')
nanmean(map{1}.fitness(:))
sum(~isnan(map{1}.fitness(:)))
nanmean(map{2}.fitness(:))
sum(~isnan(map{2}.fitness(:)))

%%
for i=1:2
observations = reshape(map{i}.genes,[],d.dof);
fitness = reshape(map{i}.fitness,[],1);
features = reshape(map{i}.features,[],2);

nans = all(isnan(observations'));
observations(nans,:) = [];
fitness(nans) = [];
features(nans,:) = [];


figh = figure(4+i);
fitcolor = [0 1 0].*fitness + [1 0 0].*(1-fitness);
fitcolor = discretize(fitcolor,0:0.25:1)./5;
showPhenotype(observations,d,figh,features,fitcolor);
drawnow;
end




