% EXPERIMENT_MAPandVAESELECTION
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Apr 2020; Last revision: 01-May-2020
%
%------------- BEGIN CODE --------------

clear;clc;
%%
DOF = 16;DOMAIN = 'polygons';ALGORITHM = 'voronoi';LATENTMODEL = 'VAE';
addpath(genpath('/home/alex/IPPM'));
rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);
rmpath(genpath('latentmodels')); addpath(genpath(['latentmodels/' LATENTMODEL]));
fname = ['mapInitWithSelection' '.mat'];
d = domain(DOF);
p = defaultParamSet;

d.resolution = 64;
%% Experiment
% User Interaction Demonstration
p.numInitSamples = 32;

% Initialize initial solution set for MAP-Elites
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,~,polygons] = fitfun(initSamples,d);

m = cfgLatentModel('data/workdir',d.resolution);
poemCfg = poemParamSet(p,m);
poemCfg.numInitSamples                = p.numInitSamples;
poemCfg.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);

poemCfg.numIterations = 2;

poemCfg.selectionMethod               = 'all';
[map{1}, config{1}, stats{1}] = poem(initSamples,polygons,fitness,poemCfg,d,2);

poemCfg.selectionMethod               = 'ascend';
[map{2}, config{2}, stats{2}] = poem(initSamples,polygons,fitness,poemCfg,d,3);

save(fname);

disp('Optimization Done');
