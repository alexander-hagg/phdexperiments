%runUDI - Example usage script of User Driven Illumination
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jun 2019; Last revision: 18-Jun-2019


%------------- BEGIN CODE --------------

%% ------------- EXPERIMENT SETUP --------------
clear;clc;
domainname = 'FOILFFD';
systemInit;
d = ffd_Domain('/tmp/foilffd');
p = defaultParamSet(d);
p.nInitialSamples   = 32;
p.nAdditionalSamples= 8;
p.nTotalSamples     = 128;

%% ------------- EXPERIMENT STARTUP --------------
% Fix initial sample set to compare models
d.initialSampleSource   = 'initialSampleSource.mat';
if ~exist(d.initialSampleSource,'file')
    [observation, value]    = initialSampling(d,p.nInitialSamples);
    save(d.initialSampleSource,'observation','value');
end
d.loadInitialSamples    = true;

% Fix Sobol generator to compare models
d.commonSobolGen = scramble(sobolset(d.nDims,'Skip',1e3),'MatousekAffineOwen');

% Initialize surrogate models
for target=1:d.numSurrogateModels; d.params{target}     = feval(['params' d.surrogateName], d.dof); end

%% ------------- EXPERIMENT LOOP --------------
output = sail(p,d);


%% ------------- EXPERIMENT ANALYSIS --------------
X = reshape(output.predMapRecord(end).genes,size(output.predMapRecord(end).genes,1)*size(output.predMapRecord(end).genes,2),[]);
X = X(~any(isnan(X')),:);

[ estimatedLabels, reducedParSpace, t_distances, epsilon, coreneighbours, stats ] = dimReducedClustering( X, 'tSNE', 2 );
%net = selforgmap([5 5]);
%[net,tr] = train(net,X');
%y = net(X');
%estimatedLabels = vec2ind(y);

%% Visualization
set(0,'DefaultFigureWindowStyle','docked')
predMap = output.predMapRecord(end);
clusterMap = predMap.fitness; clusterMap(~isnan(predMap.fitness(:))) = estimatedLabels;
figure(7); viewMap(clusterMap,d);title('Genetic Clusters');

figure(4); clf;
[h(1), h(2)] = viewMap(-predMap.fitness, d, predMap.edges); title('Predicted Fitness')
caxis([0 6]);
colormap(h(1),parula(16));
drawnow;
figure(5); clf;
[h(3), h(4)] = viewMap(predMap.confidence, d, predMap.edges); title('Prediction Model Variance')
caxis([0 1]);
colormap(h(3),parula(32));
%%
diversityMap = getDiversityMap(predMap.genes,1);
figure(6); viewMap(diversityMap,d);title('Genetic Diversity');
colormap(parula(16));
caxis([0 0.5]);
%%
%mapViewer(output.predMap(end),d)



%------------- END CODE --------------
