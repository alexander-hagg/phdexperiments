function m = cfgLatentModel(workDir,resolution, varargin)
%CFGLATENTMODEL Configuration of Autoencoder
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2019; Last revision: 13-Nov-2019
%
%------------- BEGIN CODE --------------

lossFunction                  = 'default';           % default binaryCrossEntropy
rmpath(genpath('lossFunctions')); addpath(genpath(['lossFunctions/loss_' lossFunction]));
cfg                         = 'convDefault'; %convDefault convSmall
latentDim                     = 5; %10
if nargin > 2
    latentDim = varargin{1};
    disp(['Latent dims: ' int2str(latentDim)]);
end
numFilters                    = 8;
trainPerc                     = 0.8;
numEpochs                     = 50;%350;
maxBatchSize                  = 128;
learnRate                     = 1e-3;

m                             = feval(cfg,latentDim,resolution,numFilters);
m.predict                     = @(phenotypes,m) getLatent(phenotypes,m.model)';
m.uncertainty                 = @(genomes,latent,manifold,getPhenotypes) getReconstructionError(genomes,latent,manifold,getPhenotype);
m.train                       = @(phenotypes) trainLatentModel(m,getDataPoly(phenotypes,workDir,resolution,trainPerc),numEpochs,maxBatchSize,learnRate);


end

%------------- END CODE --------------
