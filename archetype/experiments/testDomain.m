% test Domain - test domain
%
% Please include the following references in any publication using this code. For Bibtex please see the end of this file.
%
% Hagg, A., Asteroth, A. and B??ck, T., 2018, September. Prototype discovery using quality-diversity. In International Conference on Parallel Problem Solving from Nature (pp. 500-511). Springer, Cham.
% Hagg, A., Asteroth, A. and B??ck, T., 2019, July. Modeling user selection in quality diversity. In Proceedings of the Genetic and Evolutionary Computation Conference (pp. 116-124). ACM.
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Mar 2020; Last revision: 28-Mar-2020

%------------- BEGIN CODE --------------

clear;clc;
DOF = 16;
DOMAIN = 'footprints';
ALGORITHM = 'voronoi';

addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);

d = domain(DOF);
p = defaultParamSet;
d.fitfun = d.fitfunPointSymm; % Multimodal function

%%

numInitSamples = 10;
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+numInitSamples)-1,:)+d.ranges(:,1)';


[fitness,polyshapes,booleanMap,meanMovMeanE,meanMovMeanU,stdU] = getFitness(genomes,d);
%[fitness,polygons] = d.fitfun(initSamples);
%p.categorize = @(geno,pheno,p,d) categorize(pheno,d);












