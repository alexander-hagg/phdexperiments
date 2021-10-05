function m = getAEConfig(workDir,resolution, varargin)
%GETAECONFIG Configuration of Autoencoder
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2019; Last revision: 13-Nov-2019
%
%------------- BEGIN CODE --------------

lossFunction                  = 'default';           % default binaryCrossEntropy
rmpath(genpath('generators/conv_AE/lossFunctions')); addpath(genpath(['generators/conv_AE/lossFunctions/loss_' lossFunction]));
cfgAE                         = 'structureAE'; %structureAE structureAE_small
latentDim                     = 2; %10
if nargin > 2
    latentDim = varargin{1};
    disp(['Latent dims: ' int2str(latentDim)]);
end
numFilters                    = 8;
trainPerc                     = 0.95;
numEpochs                     = 350;
maxBatchSize                  = 128;
learnRate                     = 1e-3;

m                             = feval(cfgAE,latentDim,resolution,numFilters);
m.predict                     = @(phenotypes,m) getAELatent(phenotypes,m.model)';
m.uncertainty                 = @(genomes,latent,manifold,getPhenotypes) getReconstructionError(genomes,latent,manifold,getPhenotype);
m.train                       = @(phenotypes) trainAE(m,getDataPoly(phenotypes,workDir,resolution,trainPerc),numEpochs,maxBatchSize,learnRate);


end

%------------- END CODE --------------
