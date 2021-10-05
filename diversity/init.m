DOF = 16;DOMAIN = 'npoly_ffd_imagePhenotype';
ALGORITHM = 'grid';
addpath(genpath('/home/alex/diversity'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);

d = domain(DOF);
p = defaultParamSet;
d.fitfun = d.fitfunPointSymm; % Multimodal function

p.nGens = 2048;
p.featureResolution = [20 20];
p.numInitSamples = p.featureResolution(1)*p.featureResolution(2);
p.nChildren = p.numInitSamples;

axialMultiplier = 0.5;
% radialMultiplier = 0.25;
% d.ranges(:,1) = [(1-axialMultiplier)*ones(d.dof/2,1);-radialMultiplier*pi*ones(d.dof/2,1)];
% d.ranges(:,2) = [(1+axialMultiplier)*ones(d.dof/2,1);radialMultiplier*pi*ones(d.dof/2,1)];

% Initial set
% sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
% initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
% [fitness,polygons] = d.fitfun(initSamples);
% features = categorize(polygons, d);


