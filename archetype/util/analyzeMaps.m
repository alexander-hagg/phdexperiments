function [truemap,errors,filled,medianFitness] = analyzeMaps(map,d,p,varargin)
%GETTRUEMAPS Summary of this function goes here
%   [truemap,errors,filled,medianFitness] = analyzeMaps(map,d,p,[opt: getTrueMap],[opt: getHalfRes])
errors = 0; truemap = map;
getTrueMap = false; if nargin > 3; getTrueMap = varargin{1}; end
resolution = p.resolution; if nargin > 4; resolution = varargin{2}; end; p.resolution = resolution;
disp(['Setting resolution of output true map to (' int2str(p.resolution) ')']);

tGenes = reshape(map.genes,[],16);
predFitness = reshape(map.fitness,[],1);
predFeatures = reshape(map.features,[],2);
predFitness(any(isnan(tGenes'))) = [];
predFeatures(any(isnan(tGenes')),:) = [];
nanVals = any(isnan(tGenes'));
tGenes(nanVals,:) = [];
if getTrueMap
    disp('Calculating true map (for surrogate-assisted methods like SAIL and SPHEN)');
    [tFitness,features] = d.fitfun(tGenes);
else
    tFitness = predFitness;
    features = predFeatures;
end
% Create true map
truemap                                       = createMap(d, p);
[replaced, replacement, features]             = nicheCompete(tGenes, tFitness, truemap, d, p, features);
truemap                                       = updateMap(replaced,replacement,truemap,tFitness,tGenes,features);

% Get prediction errors
% errors(1) = mean(abs(predFitness-tFitness)./tFitness)*100;
% errors(2) = mean(abs(predFeatures(:,1)-tFitness)./features(:,1))*100;
% errors(3) = mean(abs(predFeatures(:,2)-tFitness)./features(:,2))*100;

errors(1) = sqrt(immse(predFitness,tFitness));
errors(2) = sqrt(immse(predFeatures(:,1),features(:,1)));
errors(3) = sqrt(immse(predFeatures(:,2),features(:,2)));    

filled = 100*sum(~isnan(truemap.fitness(:)))/numel(truemap.fitness(:));
medianFitness = nanmedian(truemap.fitness(:));




end

