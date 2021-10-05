function [ estimatedLabels, reducedParSpace, t_distances, epsilon, coreneighbours, stats ] = dimReducedClustering( X, dimreduxFcn, no_dims, varargin )
%%DIMREDUCEDCLUSTERING dimensionality reduction and (DBSCAN) clustering
%reduction. t-SNE parameters can be adjusted, DBSCAN parameters are
%automatically determined
%
% EXAMPLE ??

unscaled_X = X;
X = X - min(X(:));
X = X / max(X(:));

perplexity = min(50,floor(size(X,1)/2));
initial_dims = min(500,size(X,2));

if any(strcmpi(dimreduxFcn, {'none'}))
    reducedParSpace = unscaled_X;
else
    if any(strcmpi(dimreduxFcn, {'GPLVM', 'CFA'}))
        [reducedParSpace, mapping] = compute_mapping(unscaled_X, dimreduxFcn, no_dims);
    elseif any(strcmpi(dimreduxFcn, {'Isomap', 'LandmarkIsomap', 'LLE', 'Laplacian', 'MVU', 'CCA', 'FastMVU', 'LPP', 'NPE', 'LLTSA'}))
        [reducedParSpace, mapping] = compute_mapping(X, dimreduxFcn, no_dims, 'adaptive', 'k', perplexity);
    elseif any(strcmpi(dimreduxFcn, {'tSNE'}))
        disp(['t-SNE with perplexity: ' int2str(perplexity)]);        
        %[reducedParSpace] = fast_tsne(X, no_dims, size(X,2),perplexity);
        [reducedParSpace, mapping] = compute_mapping(X, dimreduxFcn, no_dims, initial_dims, perplexity);
    else
        [reducedParSpace, mapping] = compute_mapping(X, dimreduxFcn, no_dims);
    end
end

coreneighbours      = max(2 * no_dims,3); %Rule of thumb
[~,t_distances]     = knnsearch(reducedParSpace,reducedParSpace,'K',coreneighbours+1);
t_distances(:,1)    = [];
t_distances         = sort(t_distances(:));
[maxVal ,maxID]     = getElbow(t_distances);
epsilon             = maxVal;

tic; [~,estimatedLabels,cen] = dbscan(reducedParSpace', epsilon, coreneighbours); t = toc;
disp(['Elapsed time is ' num2str(t) ' seconds']);

y = sum(estimatedLabels==estimatedLabels');
[stats.uniqid,ia] = unique(estimatedLabels);
stats.sizes = y(ia);

end

