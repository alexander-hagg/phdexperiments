clear;clc;
DOF = 16;
DOMAIN = 'footprints';
ALGORITHM = 'voronoi';

addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);

d = domain(DOF);
p = defaultParamSet;
d.fitfun = d.fitfunPointSymm; % Multimodal function


sc = 0.2;
experimentName = 'SimRes'
initSamples = [reshape([ones(1,d.dof/4);sc.*ones(1,d.dof/4)],1,d.dof/2), zeros(1,d.dof/2)];    

d.resolution = 64;
[fitness,polyshapes,booleanMap,movMeanE,movMeanMaxU,stdE,stdU] = getFitness(initSamples,d)
system(['mv ' d.workdir '1 ' d.workdir 'lowRes'])
d.resolution = 128;
[fitnessHR,polyshapesHR,booleanMapHR,movMeanEHR,movMeanMaxUHR,stdEHR,stdUHR] = getFitness(initSamples,d,true)
system(['mv ' d.workdir '1 ' d.workdir 'highRes'])

%%
