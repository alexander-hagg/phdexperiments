function [ estimatedLabels, reducedParSpace, t_distances, epsilon, coreneighbours, stats ] = tsne_dbscan( input, varargin )
%%TSNE_DBSCAN Clustering with DBSCAN, supported with t-SNE dimensionality
%reduction. t-SNE parameters can be adjusted, DBSCAN parameters are
%automatically determined
%
% EXAMPLE
% Create 100 Gaussian clusters (20 dimensions)
% numSets = 100;
% dims = 20;
% setPts = 3+randi(17,numSets,1);
% input = []; rlLab = [];
% for pp=1:length(setPts)
%     width = 3*rand(1);
%     input = [input; randi(10,1,dims)-5 + width*randn(setPts(pp),dims)];
%     rlLab = [rlLab; pp*ones(setPts(pp),1)];
% end
% realLabels = rlLab;
% [estimatedLabels, reducedParSpace, distances, epsilon, ~] = tsne_dbscan( input, 'perplexity', 12 );
% randIndex = rand_index(realLabels, estimatedLabels, 'adjusted');

parse = inputParser;
parse.addRequired('input');
parse.addOptional('no_dims', 2);
parse.addOptional('initial_dims', size(input,2));
parse.addOptional('theta', 0.5);
parse.addOptional('alg', 'svd');
parse.addOptional('max_iter', 1000);
parse.addOptional('perplexity', 10);
parse.addOptional('coreneighbours', NaN); % You can add domain knowledge here


parse.parse(input, varargin{:});
no_dims             = parse.Results.no_dims;
initial_dims        = parse.Results.initial_dims;
theta               = parse.Results.theta;
alg                 = parse.Results.alg;
max_iter            = parse.Results.max_iter;
perplexity          = parse.Results.perplexity;
coreneighbours      = parse.Results.coreneighbours;
if isnan(coreneighbours); coreneighbours = max(2 * no_dims,3);end %Rule of thumb
reducedParSpace     = fast_tsne(input, no_dims, initial_dims, perplexity, theta, alg, max_iter);


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

