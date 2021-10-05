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
showPhenotype(initSamples,d)
%%
m = cfgLatentModel('data/workdir',d.resolution);
poemCfg = poemParamSet(p,m);
poemCfg.numInitSamples                = p.numInitSamples;
poemCfg.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);

poemCfg.numIterations = 3;

poemCfg.selectionMethod               = 'all';
[map{1}, config{1}, stats{1}] = poem(initSamples,polygons,fitness,poemCfg,d,2);
save(fname);

disp('Optimization Done');

%%
mean(fitness(:))
nanmean(stats{1}.maps{1}.fitness(:))
nanmean(stats{1}.maps{2}.fitness(:))

observations = reshape(map{1}.genes,[],d.dof);
fitness = reshape(map{1}.fitness,[],1);
features = reshape(map{1}.features,[],2);

nans = all(isnan(observations'));
observations(nans,:) = [];
fitness(nans) = [];
features(nans,:) = [];


figh = figure(5);
fitcolor = [0 1 0].*fitness + [1 0 0].*(1-fitness);
fitcolor = discretize(fitcolor,0:0.25:1)./5;
showPhenotype(observations,d,figh,features,fitcolor);
drawnow;








