function [genes,fitness,features,bins] = extractMap(map)
%EXTRACTMAP Extract MAP-Elites map content into a usable form
%   [genes,fitness,features,bin] = extractMap(map)

genes = reshape(map.genes,[],size(map.genes,3));
fitness = reshape(map.fitness,[],size(map.fitness,3));
features = reshape(map.features,[],size(map.features,3));
nans = all(isnan(genes'));
genes(nans,:) = [];
fitness(nans,:) = [];
features(nans,:) = [];

squashedFeatures = features;
squashedFeatures(squashedFeatures<0) = 0;
squashedFeatures(squashedFeatures>1) = 1;
for iDim = 1:size(features,2)
    bins(:,iDim) = discretize(squashedFeatures(:,iDim),map.edges{iDim});
end
end

