function cfg = paramsHSM( inputs, cfg, varargin )
%% PARAMS_HSM - Set parameters for surrogate model
%------------- INPUT PARSING -----------
parse = inputParser;
parse.addRequired('inputs');
parse.addRequired('cfg');
parse.addOptional('minPCAVariance',NaN);
parse.addOptional('numSegments',4);
parse.parse(inputs, cfg, varargin{:});
cfg.pca_min_variance = parse.Results.minPCAVariance; % NaN to disable PCA
cfg.segmentation.kcenters = parse.Results.numSegments;

%% Modeling techniques
cfg.experts.modelName = cfg.basemodel;

%% Hierarchical Segmentation
cfg.hierarchy.num_layers = 2;
cfg.segmentation.min_samples = 5;
cfg.segmentation.method = 'kmeans';

end
