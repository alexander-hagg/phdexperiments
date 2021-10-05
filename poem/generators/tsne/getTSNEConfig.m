function m = getTSNEConfig(dof)
%GETTSNECONFIG Summary of this function goes here
%   Detailed explanation goes here

m.numDims_DR              = 2;
m.numDims_ORG             = dof;
m.perplexity              = 50;
m.speedQualitytradeOff    = 0.5;
%m.perplexity              = min(floor(size(samples,1)*0.33),50);

m.predict = @(phenotypes,m) predictSimspace(phenotypes,m.model)';
m.uncertainty = @(genomes,latent,manifold,getPhenotype) getVariance(genomes,latent,manifold,getPhenotype);
m.train = @(samples,m) trainConstraintModel(samples);

end

