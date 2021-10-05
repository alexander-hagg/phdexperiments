function uncertainty = getVariance(genomes,latent,manifold,getPhenotype)
%GETVARIANCE Summary of this function goes here
%   Detailed explanation goes here
phenotypes = getPhenotype(genomes);
[~,uncertainty] = predictSimspace(phenotypes,manifold.model);
uncertainty = uncertainty';
end

