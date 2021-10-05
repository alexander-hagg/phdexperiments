% EXPERIMENT_UserInteraction
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Apr 2020; Last revision: 28-Apr-2020
%
%------------- BEGIN CODE --------------

clear;clc;
%%
DOF = 16;DOMAIN = 'polygons';ALGORITHM = 'voronoi';LATENTMODEL = 'VAE';
addpath(genpath('/home/alex/IPPM'));
rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);
rmpath(genpath('latentmodels')); addpath(genpath(['latentmodels/' LATENTMODEL]));
fname = ['ippm' '.mat'];
d = domain(DOF);
p = defaultParamSet;

d.resolution = 16;
%% Experiment
% User Interaction Demonstration
p.numInitSamples = 16;

% Initialize initial solution set for MAP-Elites
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,~,polygons] = fitfun(initSamples,d);

m = cfgLatentModel('data/workdir',d.resolution);
poemCfg = poemParamSet(p,m);
poemCfg.numInitSamples                = p.numInitSamples;
poemCfg.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);

poemCfg.numIterations = 3;
[map, config, stats] = poem(initSamples,polygons,fitness,poemCfg,d);

%save(fname);

disp('Optimization Done');
