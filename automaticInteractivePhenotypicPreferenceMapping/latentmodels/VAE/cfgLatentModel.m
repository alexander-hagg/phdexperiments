function m = cfgLatentModel(workDir,resolution, varargin)
%CFGLATENTMODEL Configuration of Variational Autoencoder
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2019; Last revision: 13-Nov-2019
%
%------------- BEGIN CODE --------------

lossFunction                  = 'default';           % default binaryCrossEntropy
rmpath(genpath('latentmodels/VAE/lossFunctions')); addpath(genpath(['latentmodels/VAE/lossFunctions/loss_' lossFunction]));
cfg                         = 'convDefault'; %convDefault convSmall
latentDim                     = 5;
if nargin > 2
    latentDim = varargin{1};
    disp(['Latent dims: ' int2str(latentDim)]);
end

numFilters                    = 8;
trainPerc                     = 1.00;
numEpochs                     = 500;
maxBatchSize                  = 128;
learnRate                     = 1e-3;

filterSize                    = 3;
if nargin > 3
    filterSize = varargin{2};
    disp(['filterSize: ' int2str(filterSize)]);
end

stride                        = 1;
if nargin > 4
    stride = varargin{3};
    disp(['stride: ' int2str(stride)]);
end

m                             = feval(cfg,latentDim,resolution,numFilters,filterSize,stride);
m.predict                     = @(phenotypes,m) getLatent(phenotypes,m.model)';
m.uncertainty                 = @(genomes,latent,manifold,getPhenotypes) getReconstructionError(genomes,latent,manifold,getPhenotype);
m.train                       = @(phenotypes) trainLatentModel(m,getDataPoly(phenotypes,workDir,resolution,trainPerc),numEpochs,maxBatchSize,learnRate);


end

%------------- END CODE --------------
